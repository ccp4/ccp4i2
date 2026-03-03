"""
phasertng_picard - CCP4i2 wrapper for PhaserTNG Picard molecular replacement.

Picard is the master molecular replacement solution protocol in PhaserTNG.
It explores space groups, translational NCS, different cell contents, and
error estimations for search models.

This wrapper uses PhilPluginScript for native PHIL support:
- inputData/outputData use CCP4i2 rich file types (from .def.xml)
- controlParameters are populated at runtime from picard's master_phil
- At execution time, shims convert rich types to PHIL values and
  master_phil.fetch() assembles a validated working.phil
"""

import os
import shutil
import subprocess
import logging

from ccp4i2.core.PhilPluginScript import PhilPluginScript
from ccp4i2.utils.phil_shims import (
    MtzFileShim,
    PdbFileListShim,
    AsuContentShim,
    DictFileShim,
    DagFileShim,
)

logger = logging.getLogger(__name__)


class phasertng_picard(PhilPluginScript):
    TASKNAME = "phasertng_picard"
    TASKCOMMAND = "phasertng.picard"

    # PHIL scopes excluded from the GUI because they are handled by
    # CCP4i2 rich file types in inputData + shims at execution time
    PHIL_EXCLUDE_SCOPES = [
        # picard shortcut inputs (handled by shims)
        "picard.hklin",
        "picard.seqin",
        "picard.xyzin",
        # phasertng-level file inputs (also handled by shims)
        "phasertng.hklin",
        "phasertng.labin",
        "phasertng.biological_unit.sequence.filename",
        "phasertng.model.filename",
        "phasertng.cluster_compound.filename",
        "phasertng.trace.filename",
        "phasertng.file_discovery",
        # DAG input (handled by DagFileShim)
        "phasertng.put_solution.dag.filename",
        # Output settings (managed by CCP4i2)
        "output",
    ]

    def get_master_phil(self):
        """Import picard's master_phil at runtime via CCTBXParser.

        PhaserTNG uses custom PHIL converters (filesystem, mtzcol, scatterer,
        etc.) that must be registered before parsing. The CCTBXParser handles
        this automatically.
        """
        try:
            from phasertng.programs import picard
            from iotbx.cli_parser import CCTBXParser

            parser = CCTBXParser(
                program_class=picard.Program,
                logger=None,
                parse_phil=False,
            )
            return parser.master_phil
        except ImportError as e:
            logger.warning("Cannot import phasertng: %s", e)
            return None

    def get_shim_definitions(self):
        return [
            MtzFileShim("HKLIN", "picard.hklin"),
            PdbFileListShim("XYZIN", "picard.xyzin"),
            AsuContentShim("ASUIN", "picard.seqin"),
            DictFileShim("DICT", "phasertng.cluster_compound.filename"),
            DagFileShim("DAGIN", "phasertng.put_solution.dag.filename"),
        ]

    def get_command_target(self):
        """Return the phasertng.picard entry point."""
        # Try to find the command on PATH first
        cmd = shutil.which("phasertng.picard")
        if cmd:
            return cmd

        # Fallback: use the command_line module directly
        try:
            from phasertng.command_line import picard
            return picard.__file__
        except ImportError:
            return "phasertng.picard"

    def makeCommandAndScript(self):
        """Build working.phil and construct the command line.

        phasertng.picard can be invoked as:
            phasertng.picard working.phil
        or:
            ccp4-python -m phasertng.command_line.picard working.phil
        """
        phil_path = self.build_working_phil()

        # Use the command directly if on PATH
        cmd = shutil.which("phasertng.picard")
        if cmd:
            self.TASKCOMMAND = cmd
            self.appendCommandLine([phil_path])
        else:
            # Fallback to python module invocation
            self.TASKCOMMAND = "ccp4-python"
            self.appendCommandLine(["-m", "phasertng.command_line.picard", phil_path])

        return PhilPluginScript.SUCCEEDED

    def _elevate_to_job_dir(self, source_path):
        """Copy a file into the job root directory if it isn't there already.

        Output files must live in the job directory (getWorkDirectory()) for
        glean_job_files() to register them in the database.  Programs that
        write to subdirectories (like phasertng_picard/) need their outputs
        copied up before calling setFullPath().
        """
        work_dir = str(self.getWorkDirectory())
        dest = os.path.join(work_dir, os.path.basename(source_path))
        if os.path.normpath(source_path) != os.path.normpath(dest):
            shutil.copy2(source_path, dest)
        return dest

    def processOutputFiles(self):
        """Harvest phasertng output files into outputData.

        phasertng.picard writes output to a database subdirectory
        (phasertng_picard/) rather than directly in the work directory.
        The best solution files are:
          phasertng_picard/best.1.coordinates.pdb
          phasertng_picard/best.1.dag.cards

        Files found in subdirectories are copied into the job root
        before calling setFullPath(), because glean_job_files() only
        registers files that live directly in the job directory.
        """
        import glob

        work_dir = str(self.getWorkDirectory())
        db_dir = os.path.join(work_dir, "phasertng_picard")
        out = self.container.outputData

        # Look for best solution PDB in the database directory first,
        # then fall back to work directory
        search_dirs = [db_dir, work_dir] if os.path.isdir(db_dir) else [work_dir]

        for search_dir in search_dirs:
            for pattern in ("best.*.coordinates.pdb", "*.pdb", "*.cif"):
                matches = glob.glob(os.path.join(search_dir, pattern))
                if matches:
                    try:
                        dest = self._elevate_to_job_dir(matches[0])
                        out.XYZOUT.setFullPath(dest)
                        out.XYZOUT.annotation = "Positioned coordinates from molecular replacement"
                        logger.info("Set XYZOUT to %s", dest)
                    except Exception as e:
                        logger.debug("Could not set XYZOUT: %s", e)
                    break
            if out.XYZOUT.isSet():
                break

        # Look for output MTZ (in reflections_and_models subdirectory or work dir)
        refl_dir = os.path.join(db_dir, "reflections_and_models")
        mtz_search_dirs = (
            [refl_dir, db_dir, work_dir]
            if os.path.isdir(refl_dir) else search_dirs
        )
        for search_dir in mtz_search_dirs:
            mtz_files = glob.glob(os.path.join(search_dir, "*.mtz"))
            if mtz_files:
                try:
                    dest = self._elevate_to_job_dir(mtz_files[0])
                    out.HKLOUT.setFullPath(dest)
                    out.HKLOUT.annotation = "Reflections with MR phases"
                    logger.info("Set HKLOUT to %s", dest)
                except Exception as e:
                    logger.debug("Could not set HKLOUT: %s", e)
                break

        # Look for best DAG cards file (for chaining into subsequent runs)
        dag_pattern = os.path.join(db_dir, "best.*.dag.cards")
        dag_files = sorted(glob.glob(dag_pattern))
        if dag_files:
            try:
                dest = self._elevate_to_job_dir(dag_files[0])
                out.DAGOUT.setFullPath(dest)
                out.DAGOUT.annotation = "Solution state for incremental MR"
                logger.info("Set DAGOUT to %s", dest)
            except Exception as e:
                logger.debug("Could not set DAGOUT: %s", e)

        # Generate sigma-A weighted map coefficients from the MR solution
        if out.XYZOUT.isSet():
            self._run_sigmaa(work_dir, out)

        return PhilPluginScript.SUCCEEDED

    def _run_sigmaa(self, work_dir, out):
        """Run servalcat sigmaa to produce map coefficients from the MR solution.

        Uses the original HKLIN (observed data) and XYZOUT (positioned model)
        to calculate sigma-A weighted 2mFo-DFc and mFo-DFc map coefficients.
        The resulting MTZ is split into FPHIOUT and DIFFPHIOUT via splitHklout().
        """
        hklin_path = self.container.inputData.HKLIN.getFullPath()
        xyzout_path = out.XYZOUT.getFullPath()

        if not hklin_path or not xyzout_path:
            logger.info("Skipping sigmaa: HKLIN or XYZOUT not available")
            return

        sigmaa_prefix = os.path.join(work_dir, "sigmaa")
        sigmaa_mtz = sigmaa_prefix + ".mtz"

        cmd = [
            "servalcat", "sigmaa",
            "--hklin", str(hklin_path),
            "--model", str(xyzout_path),
            "--source", "xray",
            "-o", sigmaa_prefix,
        ]

        logger.info("Running servalcat sigmaa: %s", " ".join(cmd))
        try:
            result = subprocess.run(
                cmd, cwd=work_dir,
                capture_output=True, text=True, timeout=300,
            )
            if result.returncode != 0:
                logger.warning("servalcat sigmaa failed (rc=%d): %s",
                               result.returncode, result.stderr)
                return
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            logger.warning("servalcat sigmaa not available or timed out: %s", e)
            return

        if not os.path.exists(sigmaa_mtz):
            logger.warning("servalcat sigmaa produced no output MTZ")
            return

        # Split the combined MTZ into individual map coefficient files
        out.FPHIOUT.annotation = "Sigma-A weighted 2mFo-DFc map coefficients"
        out.FPHIOUT.subType = 1
        out.DIFFPHIOUT.annotation = "Sigma-A weighted mFo-DFc map coefficients"
        out.DIFFPHIOUT.subType = 2

        error = self.splitHklout(
            ["FPHIOUT", "DIFFPHIOUT"],
            ["FWT,PHWT", "DELFWT,PHDELWT"],
            inFile=sigmaa_mtz,
        )
        if error and error.maxSeverity() > 1:
            logger.warning("Failed to split sigmaa MTZ: %s", error)
        else:
            logger.info("Map coefficients generated from sigmaa")
