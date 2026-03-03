"""
phasertng_picard - CCP4i2 wrapper for PhaserTNG molecular replacement.

This is a unified wrapper that dispatches to either:
- **Picard** (default): The master MR solution protocol. Explores space groups,
  translational NCS, cell contents, and error estimations for search models.
- **Riker** (when FIXED is set): The "refine and continue" protocol. Takes
  coordinates already placed by a previous MR run, refines them, and optionally
  searches for additional components.

The dispatch is automatic: if the user provides a fixed model (FIXED input),
we run phasertng.riker; otherwise we run phasertng.picard.

This wrapper uses PhilPluginScript for native PHIL support:
- inputData/outputData use CCP4i2 rich file types (from .def.xml)
- controlParameters are populated at runtime from picard's master_phil
  (the phasertng.* scope is shared between picard and riker)
- At execution time, shims convert rich types to PHIL values and the
  appropriate master_phil.fetch() assembles a validated working.phil

Note: When running Riker, riker.control.xref_software is forced to 'refmac'
because CCP4 does not bundle Phenix.
"""

import os
import shutil
import subprocess
import logging

from ccp4i2.core.PhilPluginScript import PhilPluginScript
from ccp4i2.utils.phil_shims import (
    MtzFileShim,
    PdbFileShim,
    PdbFileListShim,
    AsuContentShim,
    DictFileShim,
    DagFileShim,
)

logger = logging.getLogger(__name__)


class phasertng_picard(PhilPluginScript):
    TASKNAME = "phasertng_picard"
    TASKCOMMAND = "phasertng.picard"

    # PHIL scopes excluded from the GUI.
    #
    # Two categories:
    # 1. Shim-handled — file I/O paths handled by CCP4i2 rich types + shims.
    #    These must never appear as text fields in controlParameters.
    # 2. Forced/internal — scopes that CCP4i2 manages directly (e.g., forcing
    #    xref_software=refmac) or pipeline bookkeeping with no expert_level
    #    marking in PHIL (mode, title, overwrite, suite, output).
    #
    # Scopes with high expert_level in PHIL (tncs, anisotropy, data, etc.)
    # are NOT excluded here — they are filtered by the PHIL_EXPERT_LEVEL
    # mechanism in PhilPluginScript, which controls both UI visibility and
    # serialization to working.phil.
    PHIL_EXCLUDE_SCOPES = [
        # picard shortcut inputs (handled by shims)
        "picard.hklin",
        "picard.seqin",
        "picard.xyzin",
        # riker shortcut inputs (handled by shims when FIXED is set)
        "riker.hklin",
        "riker.fixed",
        "riker.xyzin",
        "riker.seqin",
        # Refinement software choice (forced to refmac in CCP4)
        "riker.control.xref_software",
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
        # Pipeline bookkeeping — no expert_level in PHIL, not user params
        "phasertng.dag",
        "phasertng.notifications",
        "phasertng.mode",
        "phasertng.title",
        "phasertng.overwrite",
        "phasertng.suite",
        # Output settings (managed by CCP4i2)
        "output",
    ]

    def _prepare_model(self, file_obj, tag, work_dir):
        """Prepare a CPdbDataFile for PhaserTNG, returning a PDB path.

        Handles two concerns:
        1. Selection — if the user selected specific chains/residues, extract
           only those atoms via getSelectedAtomsPdbFile().
        2. Format — PhaserTNG (iotbx data_manager) converts CIF→PDB by
           writing next to the source file, which fails when the source
           directory is read-only or already contains a converted file from
           a previous job.  We always produce a PDB in *this* job's
           directory so iotbx never needs to convert.

        Args:
            file_obj: CPdbDataFile instance
            tag: short label for the output filename (e.g. "XYZIN_0")
            work_dir: job working directory

        Returns:
            Path to a PDB file in work_dir, or the original path if no
            conversion or selection is needed.

        Raises:
            RuntimeError: if conversion or selection extraction fails.
        """
        if file_obj is None or not file_obj.isSet():
            return None

        src = str(file_obj.getFullPath())
        is_cif = src.lower().endswith((".cif", ".mmcif"))
        has_sel = file_obj.isSelectionSet()

        if not is_cif and not has_sel:
            return src  # already PDB, no selection — use as-is

        pdb_path = os.path.join(work_dir, f"{tag}.pdb")

        if has_sel:
            # getSelectedAtomsPdbFile → writeSelection, which auto-detects
            # output format from the .pdb extension, so CIF→PDB is implicit.
            file_obj.loadFile()
            rc = file_obj.getSelectedAtomsPdbFile(pdb_path)
            if rc != 0:
                raise RuntimeError(
                    f"Failed to extract selected atoms from {src}"
                )
            logger.info("Wrote selected atoms: %s → %s", src, pdb_path)
        else:
            # CIF with no selection — convert with gemmi
            import gemmi
            logger.info("Converting CIF→PDB: %s → %s", src, pdb_path)
            st = gemmi.read_structure(src)
            st.write_pdb(pdb_path)

        return pdb_path

    def processInputFiles(self):
        """Prepare model inputs: apply selections and convert CIF→PDB.

        Stores a mapping from original path → prepared PDB path in
        self._model_path_map, which build_working_phil() uses to
        rewrite the PHIL values emitted by the shims.
        """
        from ccp4i2.core.CCP4ErrorHandling import CErrorReport

        error = CErrorReport()
        work_dir = str(self.getWorkDirectory())
        path_map = {}

        # XYZIN is a CList of CPdbDataFile
        for i, item in enumerate(self.container.inputData.XYZIN):
            src = str(item.getFullPath()) if item and item.isSet() else None
            if src:
                try:
                    prepared = self._prepare_model(item, f"XYZIN_{i}", work_dir)
                    if prepared and prepared != src:
                        path_map[src] = prepared
                except Exception as e:
                    error.append(
                        klass=self.__class__.__name__, code=201,
                        details=f"Failed to prepare search model {src}: {e}",
                        name="processInputFiles", severity=4,
                    )
                    return error

        # FIXED (riker mode) is a single CPdbDataFile
        try:
            fixed = self.container.inputData.FIXED
            if fixed is not None and fixed.isSet():
                src = str(fixed.getFullPath())
                try:
                    prepared = self._prepare_model(fixed, "FIXED", work_dir)
                    if prepared and prepared != src:
                        path_map[src] = prepared
                except Exception as e:
                    error.append(
                        klass=self.__class__.__name__, code=201,
                        details=f"Failed to prepare fixed model {src}: {e}",
                        name="processInputFiles", severity=4,
                    )
                    return error
        except AttributeError:
            pass

        self._model_path_map = path_map
        return error

    def _use_riker(self):
        """Return True if the user has provided a fixed model (FIXED input)."""
        try:
            fixed = self.container.inputData.FIXED
            return fixed is not None and fixed.isSet()
        except AttributeError:
            return False

    def get_master_phil(self):
        """Import picard's master_phil at runtime via CCTBXParser.

        We always use picard's master_phil for GUI population because the
        phasertng.* scope is shared between picard and riker. When running
        riker, we fetch against riker's master_phil at execution time instead.
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

    def _get_riker_master_phil(self):
        """Import riker's master_phil for execution-time validation."""
        try:
            from phasertng.programs import riker
            from iotbx.cli_parser import CCTBXParser

            parser = CCTBXParser(
                program_class=riker.Program,
                logger=None,
                parse_phil=False,
            )
            return parser.master_phil
        except ImportError as e:
            logger.warning("Cannot import phasertng riker: %s", e)
            return None

    def get_shim_definitions(self):
        """Return shims appropriate for the current dispatch mode."""
        if self._use_riker():
            return [
                MtzFileShim("HKLIN", "riker.hklin"),
                PdbFileShim("FIXED", "riker.fixed"),
                PdbFileListShim("XYZIN", "riker.xyzin"),
                AsuContentShim("ASUIN", "riker.seqin"),
                DictFileShim("DICT", "phasertng.cluster_compound.filename"),
                DagFileShim("DAGIN", "phasertng.put_solution.dag.filename"),
            ]
        else:
            return [
                MtzFileShim("HKLIN", "picard.hklin"),
                PdbFileListShim("XYZIN", "picard.xyzin"),
                AsuContentShim("ASUIN", "picard.seqin"),
                DictFileShim("DICT", "phasertng.cluster_compound.filename"),
                DagFileShim("DAGIN", "phasertng.put_solution.dag.filename"),
            ]

    def get_command_target(self):
        """Return the appropriate phasertng entry point."""
        if self._use_riker():
            cmd = shutil.which("phasertng.riker")
            if cmd:
                return cmd
            try:
                from phasertng.command_line import riker
                return riker.__file__
            except ImportError:
                return "phasertng.riker"
        else:
            cmd = shutil.which("phasertng.picard")
            if cmd:
                return cmd
            try:
                from phasertng.command_line import picard
                return picard.__file__
            except ImportError:
                return "phasertng.picard"

    def build_working_phil(self):
        """Assemble working.phil using the appropriate master_phil.

        When running riker, we use riker's master_phil for validation and
        inject xref_software=refmac (CCP4 does not bundle Phenix).
        """
        from libtbx.phil import parse

        use_riker = self._use_riker()

        if use_riker:
            master_phil = self._get_riker_master_phil()
            logger.info("Dispatching to phasertng.riker (FIXED model provided)")
        else:
            master_phil = self.get_master_phil()
            logger.info("Dispatching to phasertng.picard")

        if master_phil is None:
            self.appendErrorReport(
                203,
                "Cannot import PhaserTNG master_phil — is phasertng installed?",
            )
            phil_path = os.path.join(str(self.getWorkDirectory()), "working.phil")
            with open(phil_path, "w") as f:
                f.write("")
            return phil_path

        # Collect user-set PHIL parameters from controlParameters
        user_params = self.extract_phil_parameters()

        # Run shims to convert rich CCP4i2 types to PHIL values
        work_dir = str(self.getWorkDirectory())
        for shim in self.get_shim_definitions():
            user_params.extend(shim.convert(self.container, work_dir))

        # Substitute any model paths that were prepared in processInputFiles
        # (CIF→PDB conversion and/or atom selection extraction).
        path_map = getattr(self, "_model_path_map", {})
        if path_map:
            user_params = [
                (name, path_map.get(str(val), val))
                for name, val in user_params
            ]

        # Inject forced parameters for riker
        if use_riker:
            user_params.append(("riker.control.xref_software", "refmac"))

        # Build PHIL string
        if user_params:
            user_lines = [f"{name}={val}" for name, val in user_params]
            user_phil = parse("\n".join(user_lines))
            working_phil = master_phil.fetch(sources=[user_phil])
        else:
            working_phil = master_phil

        # Write to working.phil in job directory
        phil_path = os.path.join(work_dir, "working.phil")
        from io import StringIO
        buf = StringIO()
        working_phil.show(out=buf)
        with open(phil_path, "w") as f:
            f.write(buf.getvalue())

        return phil_path

    def makeCommandAndScript(self):
        """Build working.phil and construct the command line.

        Dispatches to either phasertng.picard or phasertng.riker based
        on whether the FIXED input is set.
        """
        phil_path = self.build_working_phil()

        if self._use_riker():
            program = "phasertng.riker"
            module = "phasertng.command_line.riker"
        else:
            program = "phasertng.picard"
            module = "phasertng.command_line.picard"

        cmd = shutil.which(program)
        if cmd:
            self.TASKCOMMAND = cmd
            self.appendCommandLine([phil_path])
        else:
            self.TASKCOMMAND = "ccp4-python"
            self.appendCommandLine(["-m", module, phil_path])

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

    def _get_db_dir(self, work_dir):
        """Return the phasertng database subdirectory.

        Picard uses phasertng_picard/, Riker uses phasertng_riker/.
        Falls back to whichever exists.
        """
        for name in ("phasertng_picard", "phasertng_riker"):
            candidate = os.path.join(work_dir, name)
            if os.path.isdir(candidate):
                return candidate
        return os.path.join(work_dir, "phasertng_picard")

    def _extract_log_errors(self, max_lines=50):
        """Scan the PhaserTNG log file for error/failure lines.

        Returns a summary string suitable for inclusion in an error report,
        or None if no notable errors found.
        """
        try:
            log_text = self.logFileText()
        except Exception:
            return None
        if not log_text:
            return None

        markers = ("Sorry:", "Error:", "FATAL", "No solution",
                   "failed", "Traceback (most recent call last)")
        hits = []
        for line in log_text.splitlines()[-max_lines:]:
            stripped = line.strip()
            if any(m in stripped for m in markers):
                hits.append(stripped)
        return "\n".join(hits) if hits else None

    def processOutputFiles(self):
        """Harvest phasertng output files into outputData.

        Both picard and riker write output to a database subdirectory
        (phasertng_picard/ or phasertng_riker/). The best solution files are:
          best.1.coordinates.pdb
          best.1.dag.cards

        Files found in subdirectories are copied into the job root
        before calling setFullPath(), because glean_job_files() only
        registers files that live directly in the job directory.

        Only best.*.coordinates.pdb is accepted as solution output —
        intermediate files (.monostructure.pdb, .trace.pdb, etc.) from
        earlier pipeline steps must NOT be harvested.
        """
        import glob

        work_dir = str(self.getWorkDirectory())
        db_dir = self._get_db_dir(work_dir)
        out = self.container.outputData

        # Look for best solution PDB — only the specific best.N.coordinates.pdb
        # pattern produced by a completed pipeline. Do NOT fall back to *.pdb
        # which would match intermediate files (models, monostructure, trace).
        search_dirs = [db_dir, work_dir] if os.path.isdir(db_dir) else [work_dir]

        annotation = (
            "Refined coordinates from molecular replacement"
            if self._use_riker()
            else "Positioned coordinates from molecular replacement"
        )

        for search_dir in search_dirs:
            matches = sorted(glob.glob(
                os.path.join(search_dir, "best.*.coordinates.pdb")
            ))
            if matches:
                try:
                    dest = self._elevate_to_job_dir(matches[0])
                    out.XYZOUT.setFullPath(dest)
                    out.XYZOUT.annotation = annotation
                    logger.info("Set XYZOUT to %s", dest)
                except Exception as e:
                    logger.debug("Could not set XYZOUT: %s", e)
                break

        if not out.XYZOUT.isSet():
            detail = "No MR solution coordinates found (best.*.coordinates.pdb)"
            log_errors = self._extract_log_errors()
            if log_errors:
                detail += f"\n\nPhaserTNG log:\n{log_errors}"
            self.appendErrorReport(201, detail)
            return PhilPluginScript.FAILED

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
                self.appendErrorReport(
                    202,
                    f"servalcat sigmaa failed (rc={result.returncode}): "
                    f"{result.stderr[:500]}",
                    severity=2,  # WARNING — non-fatal
                )
                return
        except (subprocess.TimeoutExpired, FileNotFoundError) as e:
            self.appendErrorReport(
                202, f"servalcat sigmaa not available or timed out: {e}",
                severity=2,
            )
            return

        if not os.path.exists(sigmaa_mtz):
            self.appendErrorReport(
                202, "servalcat sigmaa produced no output MTZ",
                severity=2,
            )
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
