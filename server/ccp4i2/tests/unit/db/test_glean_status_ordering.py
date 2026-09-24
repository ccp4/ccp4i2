"""Gleaning is part of the job, so the terminal status is recorded after it.

A run whose output files never reach the database has produced nothing a
later task can consume, so it is FAILED however cleanly the program exited.
That makes the *order* of the two writes load-bearing, in two ways:

* A job that gleaned must not be recorded FINISHED until it has, because
  ``get_job_report_xml`` writes any terminal-status report to disk and never
  regenerates it unasked. A report generated in the gap lists no output files
  and is then cached permanently -- and the report view refetches the moment
  it sees the status transition, so the gap is reachable, not theoretical.
* A job that failed to glean must end FAILED, and the failure must propagate,
  because ``run_job_async`` writes a belt-and-braces FINISHED after
  ``track_job`` exits and the Azure worker infers failure from the subprocess
  return code (``server/worker.py``, ``run_ccp4_analysis``).

These tests pin the ordering rather than the implementation: they assert on
the sequence of recorded calls, so a future refactor is free to move the code
around as long as gleaning still happens first.

They drive the coroutines with ``asyncio.run`` rather than
``@pytest.mark.asyncio``. pytest-asyncio is a dev dependency, so it is present
in the CCP4-free venv but *not* under ``ccp4-python`` -- where an
asyncio-marked test is silently skipped rather than failed. These guard an
invariant subtle enough to be worth running in both places.
"""

import asyncio
import uuid
from unittest.mock import AsyncMock, MagicMock, Mock, patch

import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.db import models
from ccp4i2.db.async_db_handler import (
    GLEAN_FAILED_ERROR_CODE,
    AsyncDatabaseHandler,
    record_glean_failure,
)


def _tracked_plugin(job_uuid):
    """A plugin that track_job will treat as an already-registered job."""
    plugin = Mock(spec=CPluginScript)
    plugin.get_db_job_id.return_value = job_uuid
    plugin.get_status.return_value = CPluginScript.SUCCEEDED
    plugin.container = MagicMock()
    plugin.container.outputData = MagicMock()
    plugin.errorReport = Mock()
    return plugin


@pytest.fixture
def handler():
    return AsyncDatabaseHandler(uuid.uuid4())


def test_terminal_status_is_recorded_after_gleaning(handler):
    """The FINISHED write must follow the glean, never precede it."""
    job_uuid = uuid.uuid4()
    plugin = _tracked_plugin(job_uuid)
    calls = []

    async def record_status(_uuid, status):
        calls.append(("status", status))

    async def record_glean(*_args, **_kwargs):
        calls.append(("glean", None))
        return []

    async def record_kpis(*_args, **_kwargs):
        calls.append(("kpis", None))
        return 0

    async def drive():
        async with handler.track_job(plugin):
            pass

    with patch.object(handler, "update_job_status", side_effect=record_status), \
         patch.object(handler, "glean_job_files", side_effect=record_glean), \
         patch.object(handler, "glean_performance_indicators", side_effect=record_kpis), \
         patch("ccp4i2.lib.utils.parameters.save_params.save_params_for_job"), \
         patch("ccp4i2.db.models.Job.objects.get", return_value=Mock()):
        asyncio.run(drive())

    kinds = [kind for kind, _ in calls]
    assert kinds.index("glean") < kinds.index("status", 1), (
        f"terminal status was written before gleaning: {calls}"
    )
    # RUNNING first, FINISHED last, and nothing terminal in between.
    assert calls[0] == ("status", models.Job.Status.RUNNING)
    assert calls[-1] == ("status", models.Job.Status.FINISHED)


def test_glean_failure_marks_job_failed_and_propagates(handler):
    """A job that cannot glean is FAILED, and the failure is not swallowed."""
    job_uuid = uuid.uuid4()
    plugin = _tracked_plugin(job_uuid)
    statuses = []

    async def record_status(_uuid, status):
        statuses.append(status)

    async def drive():
        async with handler.track_job(plugin):
            pass

    with patch.object(handler, "update_job_status", side_effect=record_status), \
         patch.object(handler, "glean_job_files", side_effect=OSError("disk gone")):
        with pytest.raises(OSError, match="disk gone"):
            asyncio.run(drive())

    assert models.Job.Status.FAILED in statuses
    assert models.Job.Status.FINISHED not in statuses, (
        "a job that failed to glean must never be recorded FINISHED"
    )


def test_glean_failure_reaches_the_error_report(handler):
    """The reason lands on errorReport, so diagnostic.xml can show it."""
    job_uuid = uuid.uuid4()
    plugin = _tracked_plugin(job_uuid)

    async def drive():
        async with handler.track_job(plugin):
            pass

    with patch.object(handler, "update_job_status", new_callable=AsyncMock), \
         patch.object(handler, "glean_job_files", side_effect=OSError("disk gone")):
        with pytest.raises(OSError):
            asyncio.run(drive())

    plugin.errorReport.append.assert_called_once()
    reported = plugin.errorReport.append.call_args.kwargs
    assert reported["code"] == GLEAN_FAILED_ERROR_CODE
    assert "disk gone" in reported["details"]
    assert reported["severity"] == 4


def test_record_glean_failure_never_masks_the_original_failure():
    """A plugin that cannot take an error report must not raise from here."""
    plugin = Mock()
    plugin.errorReport.append.side_effect = RuntimeError("no error report")

    record_glean_failure(plugin, OSError("disk gone"))  # must not raise


class TestSyncWrapper:
    """updateJobStatus is the legacy callback path, and used to disagree.

    It gleaned after writing FINISHED and downgraded a gleaning failure to a
    logger.warning, leaving the job FINISHED with silently missing outputs --
    the one place in the codebase that contradicted track_job and run_subjob.
    """

    def test_status_recorded_after_gleaning(self, handler):
        job_uuid = uuid.uuid4()
        calls = []

        async def record_status(_uuid, status):
            calls.append(("status", status))

        async def record_glean(*_args, **_kwargs):
            calls.append(("glean", None))
            return []

        with patch.object(handler, "update_job_status", side_effect=record_status), \
             patch.object(handler, "glean_job_files", side_effect=record_glean):
            handler.updateJobStatus(
                jobId=job_uuid,
                status=models.Job.Status.FINISHED,
                container=MagicMock(),
            )

        assert calls == [
            ("glean", None),
            ("status", models.Job.Status.FINISHED),
        ]

    def test_glean_failure_demotes_to_failed(self, handler):
        job_uuid = uuid.uuid4()
        statuses = []

        async def record_status(_uuid, status):
            statuses.append(status)

        with patch.object(handler, "update_job_status", side_effect=record_status), \
             patch.object(handler, "glean_job_files", side_effect=OSError("disk gone")):
            handler.updateJobStatus(
                jobId=job_uuid,
                status=models.Job.Status.FINISHED,
                container=MagicMock(),
            )

        assert statuses == [models.Job.Status.FAILED]
