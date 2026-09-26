"""``pandda_campaign`` -- run PanDDA 2 over a declared list of datasets.

Design note: docs/pandda-campaign-design.md, sections 3, 4 and 6.

The task is a pure function of its ``DATASETS`` list: it never asks the
database which datasets a campaign has (campaign-awareness lives in job
construction, 3.1). It stages the invocation contract's input tree inside
its own job directory (3.2), records what produced the tree (4.4), runs
``pandda2.analyse`` with the contract's argv and environment (4.1-4.2), and
verifies the output tree (4.3). Failures are classified by stderr against
the contract's catalogue, never by exit-code value (4.5).

``RUN_MODE=stage_only`` is the escape hatch of 5.5 for a run that will not
fit on this machine: the job stages and finishes, the tree is shipped, and
fan-out later takes the resulting ``pandda2_out/`` with this job's manifest
(8.4), so the receipts still trace to this job.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import threading
import time
from pathlib import Path
from xml.etree import ElementTree as ET

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ErrorHandling import CErrorReport, SEVERITY_WARNING, SEVERITY_ERROR

from ccp4i2.lib.utils.jobs import dispatch_record

from . import pandda_invocation as contract
from . import pandda_run_summary as analysis
from .pandda_staging import DatasetSpec, stage_datasets

import logging
logger = logging.getLogger(f"ccp4i2:{__name__}")


class pandda_campaign(CPluginScript):

    TASKNAME = 'pandda_campaign'
    TASKVERSION = 1.0
    TASKMODULE = 'ligands'
    TASKCOMMAND = contract.PROGRAM
    ASYNCHRONOUS = False
    PERFORMANCECLASS = 'CPanddaRunPerformance'
    WHATNEXT = []

    ERROR_CODES = {
        201: {'description': 'No datasets: the DATASETS list is empty'},
        202: {'description': 'A dataset is missing its model or reflections'},
        203: {'description': 'Two datasets carry the same label'},
        204: {'description': 'pandda2.analyse could not be found'},
        205: {'severity': SEVERITY_WARNING,
              'description': f'Fewer than {contract.MIN_DATASETS} datasets: PanDDA will not characterise a ground state'},
        206: {'severity': SEVERITY_WARNING, 'description': 'The scratch directory has little free space'},
        207: {'severity': SEVERITY_WARNING, 'description': 'Estimated peak memory exceeds this machine'},
        208: {'description': 'Staging the input tree failed'},
        210: {'description': 'PanDDA failed: out of memory'},
        211: {'description': 'PanDDA failed: free-R column not recognised'},
        212: {'description': 'PanDDA failed: ligand dictionary block not read'},
        213: {'description': 'PanDDA failed: no dataset passed its range filter'},
        214: {'description': 'PanDDA failed: scratch disk full'},
        215: {'description': 'PanDDA failed: CCP4 not visible to its environment'},
        216: {'description': 'PanDDA failed with an unrecognised error'},
        217: {'description': 'PanDDA failed: scratch directory path too long for its socket'},
        209: {'description': 'The scratch directory path is too long for Ray to open its socket there'},
        220: {'description': 'PanDDA exited cleanly but wrote no output tree'},
        221: {'severity': SEVERITY_WARNING,
              'description': 'PanDDA wrote processed datasets but no events table: a partial run'},
        222: {'severity': SEVERITY_WARNING,
              'description': 'Stage-only run: the tree is staged; run PanDDA elsewhere, then fan out from its output'},
        223: {'severity': SEVERITY_WARNING,
              'description': 'PanDDA ran to the end but analysed no dataset'},
        224: {'severity': SEVERITY_WARNING,
              'description': 'PanDDA left some datasets unanalysed'},
        225: {'description': 'this deployment registers no program run target'},
        226: {'description': 'the run target to dispatch to is not chosen or not registered'},
        227: {'description': 'the run target refused the submission'},
        228: {'description': 'PanDDA dispatched to a run target; the job waits for the reconcile'},
        229: {'description': 'the dispatched run was cancelled'},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._staging_root = None
        self._manifest = None
        self._resolved = None
        self._probe = None
        self._started = None
        self._tailer_stop = None
        self._progress = None
        self._summary = None
        self._analysis = None

    # -- interface-time methods, through the generic object_method endpoint --
    #
    # Campaign-awareness lives in job construction, not in the run (3.1):
    # these read the database to *propose* a list; the run reads only the
    # list. The reading itself lives in lib/utils/jobs/pandda_fanin.py.

    def campaignCandidates(self):
        """What "fill from campaign" would add, and why anything is left out."""
        from ccp4i2.lib.utils.jobs.pandda_fanin import campaign_candidates
        job = self._job_row()
        if job is None:
            return {"campaigns": [], "candidates": [], "skipped": [],
                    "reason": "this job is not in the database"}
        return campaign_candidates(job)

    def fillDatasetsFromCampaign(self):
        """Append every member of the campaign that has a finished dimple job
        and is not already listed, and save the parameters."""
        from ccp4i2.lib.utils.jobs.pandda_fanin import fill_datasets
        job = self._job_row()
        if job is None:
            return {"success": False, "error": "this job is not in the database"}
        result = fill_datasets(self, job)
        return result.to_dict() if hasattr(result, "to_dict") else {
            "success": result.success, "data": result.data, "error": result.error}

    def _job_row(self):
        """This job's database row: from the dbHandler context when there is
        one, else from the jobId the params header records."""
        from ccp4i2.db import models
        job_id = self.get_db_job_id() if hasattr(self, 'get_db_job_id') else None
        if not job_id:
            for name in ('input_params.xml', 'params.xml'):
                path = Path(self.workDirectory) / name
                if path.is_file():
                    try:
                        job_id = ET.parse(path).getroot().findtext('ccp4i2_header/jobId')
                    except ET.ParseError:
                        job_id = None
                    if job_id:
                        break
        if not job_id:
            return None
        return models.Job.objects.filter(uuid=job_id).first()

    # -- validation (two tiers) ------------------------------------------

    def validity(self) -> CErrorReport:
        error = super().validity()
        datasets = self.container.inputData.DATASETS
        name = f'{self.TASKNAME}.container.inputData.DATASETS'
        if len(datasets) == 0:
            error.append(klass=self.TASKNAME, code=201, details='Add at least one dataset',
                         name=name, severity=SEVERITY_ERROR)
            return error
        seen = set()
        for i, item in enumerate(datasets):
            label = str(item.DTAG) if item.DTAG.isSet() else ''
            if not item.XYZIN.isSet() or not item.HKLIN.isSet():
                error.append(klass=self.TASKNAME, code=202,
                             details=f'Dataset {i + 1} ({label or "unlabelled"}) needs a model and reflections',
                             name=f'{name}[{i}]', severity=SEVERITY_ERROR)
            if label:
                if label in seen:
                    error.append(klass=self.TASKNAME, code=203, details=f'Label {label!r} is used twice',
                                 name=f'{name}[{i}].DTAG', severity=SEVERITY_ERROR)
                seen.add(label)
        if self._mode() == 'dispatch':
            self._check_dispatch_target(error)
        return error

    def _check_dispatch_target(self, error: CErrorReport) -> None:
        """Dispatch mode needs a registered program run target. The registry
        is the deployment's (docs/run-target-dispatch.md); the desktop
        registers none, and this says so instead of failing at submit."""
        from ccp4i2.lib.dispatch import available_targets
        targets = [t for t in available_targets() if t.get('runs_programs')]
        name = self._target_name()
        field = f'{self.TASKNAME}.container.controlParameters.DISPATCH_TARGET'
        if not targets:
            error.append(klass=self.TASKNAME, code=225,
                         details='This deployment registers no run target that runs programs, so PanDDA '
                                 'cannot be dispatched from here. Run locally, or stage only and run elsewhere',
                         name=field, severity=SEVERITY_ERROR)
        elif name is None:
            error.append(klass=self.TASKNAME, code=226,
                         details='Choose the run target: ' + ', '.join(t['name'] for t in targets),
                         name=field, severity=SEVERITY_ERROR)
        elif name not in {t['name'] for t in targets}:
            error.append(klass=self.TASKNAME, code=226,
                         details=f"'{name}' is not a registered program run target; registered: "
                                 + ', '.join(t['name'] for t in targets),
                         name=field, severity=SEVERITY_ERROR)

    def runTimeValidity(self) -> CErrorReport:
        error = super().runTimeValidity()
        if error.maxSeverity() >= SEVERITY_ERROR:
            return error
        par = self.container.controlParameters
        datasets = self.container.inputData.DATASETS
        n = len(datasets)
        if self._mode() == 'dispatch':
            self._check_dispatch_target(error)
            if error.maxSeverity() >= SEVERITY_ERROR:
                return error

        if self._mode() == 'local':
            resolved = self._resolve_executable()
            if resolved is None:
                error.append(klass=self.TASKNAME, code=204,
                             details='Set PANDDA_EXECUTABLE, the PANDDA2_EXECUTABLE preference, '
                                     'or install a CCP4 that ships pandda2.analyse; '
                                     'or choose the stage-only run mode',
                             name=f'{self.TASKNAME}.container.controlParameters.PANDDA_EXECUTABLE',
                             severity=SEVERITY_ERROR)

        configured = self._configured_min_datasets()
        minimum = self._min_datasets()
        if minimum < contract.MIN_DATASETS:
            lowered = (f'lowered from {configured} to {minimum}, the number of datasets, '
                       if minimum < configured else f'set to {minimum} (PanDDA default {contract.MIN_DATASETS}) ')
            error.append(klass=self.TASKNAME, code=205,
                         details=(f'minimum datasets to characterise a ground state {lowered}so the run '
                                  'proceeds; expect the ground state to be poorly characterised and the '
                                  'events to be noisy. Fine for a test; add datasets for a real screen'),
                         name=f'{self.TASKNAME}.container.controlParameters.MIN_CHARACTERISATION_DATASETS',
                         severity=SEVERITY_WARNING)

        if self._mode() == 'local':
            scratch = self._scratch_dir()
            if not contract.scratch_fits(scratch):
                error.append(klass=self.TASKNAME, code=209,
                             details=(f'{scratch} is {len(str(scratch))} characters; Ray needs its socket path '
                                      f'under {contract.AF_UNIX_PATH_LIMIT}. Leave SCRATCH_DIR empty or '
                                      'choose a short path such as /tmp/ray'),
                             name=f'{self.TASKNAME}.container.controlParameters.SCRATCH_DIR',
                             severity=SEVERITY_ERROR)
                return error
            try:
                scratch.mkdir(parents=True, exist_ok=True)
                free_gib = shutil.disk_usage(scratch).free / 2 ** 30
                if free_gib < 20:
                    error.append(klass=self.TASKNAME, code=206,
                                 details=f'{free_gib:.0f} GB free at {scratch}; Ray needs tens of GB',
                                 name=f'{self.TASKNAME}.container.controlParameters.SCRATCH_DIR',
                                 severity=SEVERITY_WARNING)
            except OSError as e:
                error.append(klass=self.TASKNAME, code=206, details=f'{scratch}: {e}',
                             name=f'{self.TASKNAME}.container.controlParameters.SCRATCH_DIR',
                             severity=SEVERITY_WARNING)

            hint = self._sizing_hint()
            cpus = int(par.LOCAL_CPUS) if par.LOCAL_CPUS.isSet() else 4
            estimate = contract.estimate_peak_gib(n, hint['cell_volume_class'], cpus)
            total = self._physical_gib()
            if total is not None and estimate > total:
                error.append(klass=self.TASKNAME, code=207,
                             details=(f'~{estimate:.0f} GB estimated for {n} {hint["cell_volume_class"]}-cell '
                                      f'datasets at {cpus} CPUs; this machine has {total:.0f} GB. '
                                      'Fewer datasets per run is the lever; or stage only and run elsewhere'),
                             name=f'{self.TASKNAME}.container.inputData.DATASETS',
                             severity=SEVERITY_WARNING)
        return error

    # -- lifecycle ---------------------------------------------------------

    def processInputFiles(self):
        out = self.container.outputData
        self._staging_root = Path(self.workDirectory) / 'staging'
        self._resolved = self._resolve_executable() if self._mode() == 'local' else None
        self._probe = contract.probe_executable(self._resolved) if self._resolved else {}
        if self._is_harvest():
            # The run happened elsewhere; the tree staged at submit is the
            # record of what it was given. Re-staging would re-import every
            # input into this project and hand out new uuids.
            manifest_path = self._staging_root / 'manifest.json'
            if not manifest_path.is_file():
                self.appendErrorReport(208, f'no staged tree to complete from: {manifest_path} is missing')
                return CPluginScript.FAILED
            self._manifest = json.loads(manifest_path.read_text())
            out.MANIFEST.setFullPath(str(Path(self.workDirectory) / 'manifest.json'))
            out.STAGING_DIR.set(str(self._staging_root))
            out.CONTRACT_VERSION.set(contract.CONTRACT_VERSION)
            out.PERFORMANCE.nDatasets.set(len(self._manifest['datasets']))
            self._record_dispatch_output()
            return None
        job_uuid = self.get_db_job_id() if hasattr(self, 'get_db_job_id') else None
        provenance = {
            'task': self.TASKNAME,
            'run_job_uuid': str(job_uuid) if job_uuid else None,
            'run_mode': self._mode(),
            'contract': contract.CONTRACT_VERSION,
            'executable': self._resolved,
            'probe': self._probe,
            'sizing_hint': self._sizing_hint(),
        }
        try:
            self._manifest = stage_datasets(self._specs(), self._staging_root, provenance=provenance)
        except (OSError, ValueError) as e:
            self.appendErrorReport(208, str(e))
            return CPluginScript.FAILED
        # The manifest is an output in its own right (fan-out takes it with the
        # tree), so it lives at the job root where the gleaner expects one.
        manifest_path = Path(self.workDirectory) / 'manifest.json'
        shutil.copy2(self._staging_root / 'manifest.json', manifest_path)
        out.MANIFEST.setFullPath(str(manifest_path))
        out.STAGING_DIR.set(str(self._staging_root))
        out.CONTRACT_VERSION.set(contract.CONTRACT_VERSION)
        if self._resolved:
            out.PROVENANCE_EXECUTABLE.set(self._resolved)
            out.PROVENANCE_PROBE.set(json.dumps(self._probe, sort_keys=True))
        out.PERFORMANCE.nDatasets.set(len(self._manifest['datasets']))
        self.container.controlParameters.MIN_CHARACTERISATION_DATASETS.set(self._min_datasets())
        self._write_program_xml(state='staged')
        return None

    def makeCommandAndScript(self):
        argv = contract.build_argv(self._staging_root / 'datasets', self._out_dir(),
                                   self._local_cpus(), self._min_datasets())
        self.commandLine = list(argv)
        self.container.outputData.PROVENANCE_ARGV.set(' '.join([contract.PROGRAM] + argv))
        if self._resolved:
            # _prepareProcessExecution resolves TASKCOMMAND through the
            # preferences and PATH; an explicit path short-circuits that.
            self.TASKCOMMAND = self._resolved
        return None

    def startProcess(self):
        if self._mode() == 'stage_only':
            self.appendErrorReport(
                222, f'Input tree staged at {self._staging_root}. Run: {contract.PROGRAM} '
                     f'{" ".join(self.commandLine)}', stack=False)
            self._write_program_xml(state='stage_only')
            return CPluginScript.SUCCEEDED
        if self._mode() == 'dispatch':
            return self._dispatch_or_harvest()
        self._started = time.time()
        self._write_program_xml(state='running')
        stop = threading.Event()
        self._tailer_stop = stop
        tailer = threading.Thread(target=self._tail_progress, args=(stop,), daemon=True)
        tailer.start()
        try:
            result = super().startProcess()
        finally:
            stop.set()
            tailer.join(timeout=5)
        # A non-zero exit comes back as an error report, and process() stops
        # there without postProcess: classify here, or the run ends with the
        # raw exit code and a report still saying "running".
        if result:
            self._classify_failure()
        return result

    def _classify_failure(self):
        text = self._read(self.makeFileName('STDERR')) + '\n' + self._read(self.makeFileName('LOG'))
        name, code, prompt = contract.classify_failure(text)
        self.appendErrorReport(code, f'{name}: {prompt}', stack=False)
        self._record_tree()
        self._write_program_xml(state='failed', failure=name)

    def _prepareProcessExecution(self):
        prep = super()._prepareProcessExecution()
        prep['env'] = contract.build_env(prep['env'], self._scratch_dir())
        return prep

    def postProcessCheck(self, processId=None):
        status, exit_status, exit_code = super().postProcessCheck(processId)
        if status != CPluginScript.SUCCEEDED and self._mode() == 'local':
            self._classify_failure()
        return status, exit_status, exit_code

    def processOutputFiles(self):
        if self._mode() == 'stage_only':
            return CPluginScript.SUCCEEDED
        summary = self._record_tree()
        if summary['processed'] == 0:
            self.appendErrorReport(220, f'nothing under {self._out_dir()}')
            self._write_program_xml(state='failed', failure='no_output')
            return CPluginScript.FAILED
        if not summary['complete']:
            self.appendErrorReport(221, f"{summary['processed']} processed dataset(s), no events table", stack=False)
            self._write_program_xml(state='partial')
            return CPluginScript.UNSATISFACTORY
        reasons = '; '.join(summary['reasons']) or 'PanDDA gave no reason in its log'
        if summary['analysed'] == 0:
            # Exit 0, a header-only events table, empty dataset directories:
            # the run got to the end and did nothing. Say so.
            self.appendErrorReport(223, f"{summary['processed']} loaded, none analysed. {reasons}", stack=False)
            self._write_program_xml(state='empty')
            return CPluginScript.UNSATISFACTORY
        if summary['unanalysed']:
            self.appendErrorReport(
                224, f"{len(summary['unanalysed'])} of {summary['processed']} unanalysed "
                     f"({', '.join(summary['unanalysed'][:10])}). {reasons}", stack=False)
        self._write_program_xml(state='finished')
        return CPluginScript.SUCCEEDED

    # -- helpers ------------------------------------------------------------

    def _mode(self) -> str:
        par = self.container.controlParameters
        return str(par.RUN_MODE) if par.RUN_MODE.isSet() else 'local'

    # -- dispatch (axis B, docs/run-target-dispatch.md) ---------------------

    def _target_name(self):
        """DISPATCH_TARGET, else the one program target the deployment registers."""
        par = self.container.controlParameters
        if par.DISPATCH_TARGET.isSet() and str(par.DISPATCH_TARGET).strip():
            return str(par.DISPATCH_TARGET).strip().lower()
        from ccp4i2.lib.dispatch import available_targets
        programs = [t['name'] for t in available_targets() if t.get('runs_programs')]
        return programs[0] if len(programs) == 1 else None

    def _is_harvest(self) -> bool:
        return self._mode() == 'dispatch' and dispatch_record.is_harvest(self.workDirectory)

    def _record_dispatch_output(self):
        """Copy dispatch.json into the typed DISPATCH output."""
        record = dispatch_record.read_record(self.workDirectory) or {}
        out = self.container.outputData.DISPATCH
        for field, key in (('TARGET', 'target'), ('HANDLE', 'handle'), ('SUBMITTED_AT', 'submitted_at'),
                           ('STATE', 'state'), ('STDERR', 'stderr')):
            if record.get(key):
                getattr(out, field).set(str(record[key]))

    def _dispatch_or_harvest(self):
        record = dispatch_record.read_record(self.workDirectory)
        if record and record.get('state') in dispatch_record.TERMINAL_STATES:
            return self._harvest(record)
        return self._submit()

    def _submit(self):
        from ccp4i2.lib.dispatch import RunTargetError, UnknownRunTarget, get_target
        name = self._target_name()
        try:
            target = get_target(name) if name else None
        except (UnknownRunTarget, RunTargetError) as err:
            self.appendErrorReport(226, str(err), stack=False)
            return CPluginScript.FAILED
        if target is None:
            self.appendErrorReport(226, 'no program run target to dispatch to', stack=False)
            return CPluginScript.FAILED
        self._write_program_xml(state='submitting')
        try:
            handle = target.submit(self._staging_root / 'datasets', list(self.commandLine),
                                   self._out_dir(), self._sizing_hint())
        except Exception as err:  # noqa: BLE001 -- the target's failure is this job's failure
            self.appendErrorReport(227, f"run target '{name}' refused the submission: "
                                        f"{type(err).__name__}: {err}", stack=False)
            self._write_program_xml(state='failed', failure='submit_refused')
            return CPluginScript.FAILED
        dispatch_record.new_record(self.workDirectory, target=name, handle=str(handle),
                                   out_dir=str(self._out_dir()),
                                   sizing_hint=self._sizing_hint())
        self._record_dispatch_output()
        self._write_program_xml(state='dispatched')
        self.appendErrorReport(228, f"PanDDA dispatched to run target '{name}' (handle {handle}); "
                                    "the job waits until the run is reconciled", stack=False)
        return CPluginScript.DISPATCHED

    def _harvest(self, record):
        """The run ended elsewhere: complete this job from what it left."""
        self._record_dispatch_output()
        state = record.get('state')
        if state == 'succeeded':
            self._write_program_xml(state='running')
            return CPluginScript.SUCCEEDED      # processOutputFiles verifies the tree
        stderr = record.get('stderr')
        text = self._read(stderr) if stderr and Path(stderr).is_file() else ''
        if state == 'cancelled':
            self.appendErrorReport(229, f"the run on '{record.get('target')}' was cancelled", stack=False)
            self._write_program_xml(state='failed', failure='cancelled')
            return CPluginScript.INTERRUPTED
        name, code, prompt = contract.classify_failure(text)
        self.appendErrorReport(code, f'{name}: {prompt}', stack=False)
        self._record_tree()
        self._write_program_xml(state='failed', failure=name)
        return CPluginScript.FAILED

    def reconcileDispatch(self):
        """Plugin method (object_method endpoint): ask the target how the
        dispatched run is doing and act on it. Idempotent."""
        from ccp4i2.db import models
        job_uuid = self.get_db_job_id() if hasattr(self, 'get_db_job_id') else None
        if not job_uuid:
            return {"action": "error", "reason": "this job is not in the database"}
        job = models.Job.objects.get(uuid=job_uuid)
        return dispatch_record.reconcile(job)

    def _configured_min_datasets(self) -> int:
        par = self.container.controlParameters
        return int(par.MIN_CHARACTERISATION_DATASETS) if par.MIN_CHARACTERISATION_DATASETS.isSet() \
            else contract.MIN_DATASETS

    def _min_datasets(self) -> int:
        """The minimum PanDDA is given: the parameter, but never more than
        the datasets there are. A run with fewer datasets than PanDDA's
        default of 25 would otherwise skip every one of them and finish
        empty; lowering it lets a small set run, and the warning says what
        that costs. The parameter is set to the value used, so params.xml
        records what ran."""
        configured = self._configured_min_datasets()
        n = len(self.container.inputData.DATASETS)
        return min(configured, n) if n else configured

    def _local_cpus(self) -> int:
        par = self.container.controlParameters
        return int(par.LOCAL_CPUS) if par.LOCAL_CPUS.isSet() else 4

    def _out_dir(self) -> Path:
        return Path(self.workDirectory) / contract.OUT_DIR_NAME

    def _scratch_dir(self) -> Path:
        par = self.container.controlParameters
        if par.SCRATCH_DIR.isSet() and str(par.SCRATCH_DIR).strip():
            return Path(str(par.SCRATCH_DIR)).expanduser()
        token = (self.get_db_job_id() if hasattr(self, 'get_db_job_id') else None) or ''
        token = str(token).replace('-', '') or Path(self.workDirectory).name
        return contract.default_scratch_dir(self.workDirectory, token)

    def _resolve_executable(self):
        par = self.container.controlParameters
        if par.PANDDA_EXECUTABLE.isSet() and str(par.PANDDA_EXECUTABLE).strip():
            candidate = Path(str(par.PANDDA_EXECUTABLE)).expanduser()
            return str(candidate) if candidate.is_file() and os.access(candidate, os.X_OK) else None
        from ccp4i2.config.program_discovery import resolve_program
        found = resolve_program(contract.PROGRAM)
        if found:
            return found
        ccp4 = os.environ.get('CCP4')
        if ccp4:
            shipped = Path(ccp4) / 'bin' / contract.PROGRAM
            if shipped.is_file() and os.access(shipped, os.X_OK):
                return str(shipped)
        return None

    def _specs(self):
        specs = []
        for item in self.container.inputData.DATASETS:
            uuids = {}
            for role, member in (('xyzin', item.XYZIN), ('hklin', item.HKLIN), ('dict', item.DICT)):
                if member.isSet() and member.dbFileId.isSet():
                    uuids[role] = str(member.dbFileId)
            if item.PROJECT_UUID.isSet() and str(item.PROJECT_UUID).strip():
                project_uuid = str(item.PROJECT_UUID).strip()
            else:
                project_uuid = str(item.XYZIN.project) if item.XYZIN.project.isSet() else None
            source_job = (str(item.SOURCE_JOB_UUID).strip()
                          if item.SOURCE_JOB_UUID.isSet() and str(item.SOURCE_JOB_UUID).strip() else None)
            specs.append(DatasetSpec(
                label=str(item.DTAG) if item.DTAG.isSet() else f'dataset-{len(specs) + 1}',
                xyzin=Path(str(item.XYZIN.fullPath)),
                hklin=Path(str(item.HKLIN.fullPath)),
                dictionary=Path(str(item.DICT.fullPath)) if item.DICT.isSet() else None,
                project_uuid=project_uuid,
                source_job_uuid=source_job,
                source_file_uuids=uuids,
            ))
        return specs

    def _sizing_hint(self):
        import gemmi
        cells = []
        for item in self.container.inputData.DATASETS:
            if not item.HKLIN.isSet():
                continue
            try:
                mtz = gemmi.read_mtz_file(str(item.HKLIN.fullPath))
                c = mtz.cell
                cells.append((c.a, c.b, c.c, c.alpha, c.beta, c.gamma))
            except Exception:
                continue
        return contract.sizing_hint(len(self.container.inputData.DATASETS), cells)

    @staticmethod
    def _physical_gib():
        try:
            import psutil
            return psutil.virtual_memory().total / 2 ** 30
        except Exception:
            return None

    @staticmethod
    def _read(path) -> str:
        try:
            return Path(path).read_text(errors='replace')
        except OSError:
            return ''

    def _record_tree(self):
        """Count what PanDDA wrote (see ``summarise_output_tree``); set the
        output-tree fields and KPIs. Returns the summary."""
        out = self.container.outputData
        tree = self._out_dir()
        summary = contract.summarise_output_tree(tree, self._read(self.makeFileName('LOG')))
        self._summary = summary
        if summary['processed']:
            out.PANDDA_OUT_DIR.set(str(tree))
        out.PERFORMANCE.nDatasetsProcessed.set(summary['processed'])
        out.PERFORMANCE.nDatasetsAnalysed.set(summary['analysed'])
        out.PERFORMANCE.nEvents.set(summary['events'])
        # The run's own analysis, recapitulated from its tables (the bundled
        # PanDDA 2 writes no HTML summary). Never worth failing the job over.
        try:
            self._analysis = analysis.summarise_run(tree)
            out.PERFORMANCE.nSites.set(self._analysis['stats']['n_sites'])
        except Exception as e:      # noqa: BLE001
            logger.warning('PanDDA analysis tables not summarised: %s', e)
            self._analysis = None
        if self._started is not None:
            out.PERFORMANCE.wallSeconds.set(round(time.time() - self._started, 1))
        return summary

    def _tail_progress(self, stop: threading.Event):
        log = self.makeFileName('LOG')
        last = None
        while not stop.wait(2.0):
            progress = contract.parse_progress(self._read(log))
            if progress != last:
                last = progress
                self._progress = progress
                try:
                    self._write_program_xml(state='running')
                except Exception as e:      # noqa: BLE001 - a report is never worth a job
                    logger.debug('running report not written: %s', e)

    def _write_program_xml(self, state: str, failure: str = ''):
        root = ET.Element('pandda_campaign')
        ET.SubElement(root, 'state').text = state
        ET.SubElement(root, 'run_mode').text = self._mode()
        if failure:
            ET.SubElement(root, 'failure').text = failure
        if self._staging_root is not None:
            ET.SubElement(root, 'staging_dir').text = str(self._staging_root)
        ET.SubElement(root, 'out_dir').text = str(self._out_dir())
        ET.SubElement(root, 'executable').text = self._resolved or ''
        ET.SubElement(root, 'contract').text = contract.CONTRACT_VERSION
        if self._probe:
            ET.SubElement(root, 'probe').text = json.dumps(self._probe, sort_keys=True)
        record = dispatch_record.read_record(self.workDirectory)
        if record:
            node = ET.SubElement(root, 'dispatch')
            for key in ('target', 'handle', 'state', 'submitted_at', 'polled_at', 'stderr'):
                if record.get(key):
                    ET.SubElement(node, key).text = str(record[key])
        n = len(self._manifest['datasets']) if self._manifest else len(self.container.inputData.DATASETS)
        ET.SubElement(root, 'n_datasets').text = str(n)
        if self._progress:
            ET.SubElement(root, 'progress', done=str(self._progress[0]), total=str(self._progress[1]))
        perf = self.container.outputData.PERFORMANCE
        for tag, field in (('n_processed', perf.nDatasetsProcessed), ('n_analysed', perf.nDatasetsAnalysed),
                           ('n_events', perf.nEvents), ('wall_seconds', perf.wallSeconds)):
            if field.isSet():
                ET.SubElement(root, tag).text = str(field)
        if self._summary and self._summary.get('reasons'):
            reasons = ET.SubElement(root, 'reasons')
            for text in self._summary['reasons']:
                ET.SubElement(reasons, 'reason').text = text
        if self._manifest:
            datasets = ET.SubElement(root, 'datasets')
            for entry in self._manifest['datasets']:
                ET.SubElement(datasets, 'dataset', xtal=entry['xtal'], label=entry['label'],
                              dict='yes' if 'dict.cif' in entry['files'] else 'no')
        if self._analysis:
            analysis.analysis_to_xml(root, self._analysis)
        tree = ET.ElementTree(root)
        ET.indent(tree, space='  ')
        tree.write(self.makeFileName('PROGRAMXML'), encoding='utf-8', xml_declaration=True)
