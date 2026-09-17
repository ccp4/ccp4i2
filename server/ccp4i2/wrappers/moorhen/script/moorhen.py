"""The recorded Moorhen session.

Moorhen runs in a window of the app, not as a child process, so this
wrapper launches nothing. What it does is own the job that records the
session: the inputs say what the window loads (the load plan is derived
from them by ``lib/utils/jobs/interactive.py``), the window saves models
into ``MOORHEN_FILE_DROP`` through the drop endpoint, and this wrapper
harvests that directory when the session is finished.

Two ways the job reaches ``startProcess()``:

* From the app: Run opened a session without dispatching; Finish
  dispatched the job. The session is already finished, so the wrapper
  returns at once and harvests.
* From i2run: the job was dispatched immediately, so the wrapper waits
  here on the session row until a window finishes it, like ``i2run coot1``
  blocks on Coot.

See docs/moorhen-task-design.md.
"""

import logging
import os
import time
from pathlib import Path

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ModelData import CPdbDataFile

logger = logging.getLogger(f"ccp4i2:{__name__}")

DROP_DIR_NAME = "MOORHEN_FILE_DROP"
DESKTOP_LAUNCH_ENV = "CCP4I2_DESKTOP_LAUNCH"


def session_route(job):
    return f"/ccp4i2/moorhen-page/session/{job.id}"


def launch_session_window(route, environ=None):
    """Ask the desktop app to open ``route``; returns the spawned Popen or
    None if no app launch command is known or spawning failed."""
    import json
    import subprocess

    environ = os.environ if environ is None else environ
    raw = environ.get(DESKTOP_LAUNCH_ENV)
    if not raw:
        return None
    try:
        command = json.loads(raw)
        if not isinstance(command, list) or not command:
            raise ValueError("not a non-empty list")
    except (ValueError, TypeError) as err:
        logger.warning("%s is not a JSON list: %s", DESKTOP_LAUNCH_ENV, err)
        return None
    try:
        return subprocess.Popen(
            [str(part) for part in command] + ["--open-route", route],
            stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL, start_new_session=True)
    except OSError as err:
        logger.warning("Could not launch the desktop app for %s: %s", route, err)
        return None


class moorhen(CPluginScript):
    TASKNAME = "moorhen"
    TASKCOMMAND = None  # no external program: the window is the program
    WHATNEXT = ["servalcat_pipe", "prosmart_refmac", "moorhen"]
    POLL_SECONDS = 1.0

    ERROR_CODES = {
        201: {"description": "Harvesting the session's saved files failed"},
    }

    # -- session ------------------------------------------------------------

    def makeCommandAndScript(self):
        drop_dir = Path(self.getWorkDirectory()) / DROP_DIR_NAME
        drop_dir.mkdir(parents=True, exist_ok=True)
        return CPluginScript.SUCCEEDED

    def _job(self):
        """The Job row this plugin runs for, or None outside the database."""
        job_uuid = getattr(self, "_dbJobId", None)
        if not job_uuid:
            return None
        try:
            from ccp4i2.db import models

            return models.Job.objects.get(uuid=job_uuid)
        except Exception:
            logger.warning("moorhen: no Job row for %s", job_uuid)
            return None

    def startProcess(self):
        job = self._job()
        if job is None:
            # No database: nothing to wait for; harvest what is in the
            # drop directory (unit tests drive the wrapper this way).
            return CPluginScript.SUCCEEDED

        from ccp4i2.db import models
        from ccp4i2.lib.utils.jobs import interactive

        session = interactive.ensure_session_for_runner(job)
        if session.finished:
            return CPluginScript.SUCCEEDED

        self._announce_session(job)
        while True:
            time.sleep(self.POLL_SECONDS)
            session.refresh_from_db()
            if session.finished:
                return CPluginScript.SUCCEEDED
            job.refresh_from_db(fields=["status"])
            if job.status != models.Job.Status.RUNNING:
                # Cancelled from the job menu while we waited.
                return CPluginScript.INTERRUPTED

    def _announce_session(self, job):
        """Say where a window can attach, and open one if we can.

        The desktop app exports CCP4I2_DESKTOP_LAUNCH (a JSON list: its own
        executable, plus the app directory in development) into the Django
        child, so a job started from i2run can run it with --open-route;
        the app's single-instance lock forwards that to the running app,
        which opens the session window. Without it (a bare runserver, a web
        deployment) the route is logged and the app's job menu offers
        "Open session window".
        """
        route = session_route(job)
        print(f"Moorhen session open for job {job.number}: {route}")
        logger.info("moorhen session route: %s", route)
        launch_session_window(route)

    # -- harvesting ---------------------------------------------------------

    def processOutputFiles(self):
        from lxml import etree

        from ccp4i2.core import CCP4Utils
        from ccp4i2.cootbridge import harvest

        work_dir = Path(self.getWorkDirectory())
        n_models = n_dicts = 0
        try:
            n_models, n_dicts = harvest.harvest_drop_directory(
                work_dir, work_dir / DROP_DIR_NAME,
                self.container.outputData.XYZOUT,
                self.container.outputData.DICTOUT,
                self._annotate_model, self._annotate_dict)
        except Exception as err:
            self.appendErrorReport(201, str(err))
            return CPluginScript.FAILED

        for dict_file in self.container.outputData.DICTOUT[:n_dicts]:
            try:
                self.mergeDictToProjectLib(fileName=dict_file.__str__())
            except Exception:
                logger.warning("mergeDictToProjectLib failed for %s", dict_file)

        root = etree.Element("moorhen")
        etree.SubElement(root, "number_output_files").text = str(n_models)
        etree.SubElement(root, "number_output_dicts").text = str(n_dicts)
        CCP4Utils.saveEtreeToFile(root, self.makeFileName("PROGRAMXML"))

        if n_models + n_dicts > 0:
            return CPluginScript.SUCCEEDED
        # Nothing was saved: the job self-deletes rather than litter the
        # project (the Coot 0.9 convention).
        return CPluginScript.MARK_TO_DELETE

    def _annotate_model(self, item, path, meta=None):
        annotation = (meta or {}).get("annotation") or f"Moorhen output: {path.name}"
        item.annotation.set(annotation)
        item.subType.set(CPdbDataFile.SUBTYPE_MODEL)
        item.contentFlag.set(
            CPdbDataFile.CONTENT_FLAG_MMCIF if path.suffix == ".cif"
            else CPdbDataFile.CONTENT_FLAG_PDB)

    def _annotate_dict(self, item, path, meta=None):
        annotation = (meta or {}).get("annotation") or f"Moorhen ligand dictionary: {path.name}"
        item.annotation.set(annotation)
