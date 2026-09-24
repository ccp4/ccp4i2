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

from . import pandda_invocation as contract
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
        220: {'description': 'PanDDA exited cleanly but wrote no output tree'},
        221: {'severity': SEVERITY_WARNING,
              'description': 'PanDDA wrote processed datasets but no events table: a partial run'},
        222: {'severity': SEVERITY_WARNING,
              'description': 'Stage-only run: the tree is staged; run PanDDA elsewhere, then fan out from its output'},
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
        return error

    def runTimeValidity(self) -> CErrorReport:
        error = super().runTimeValidity()
        if error.maxSeverity() >= SEVERITY_ERROR:
            return error
        par = self.container.controlParameters
        datasets = self.container.inputData.DATASETS
        n = len(datasets)

        if self._mode() == 'local':
            resolved = self._resolve_executable()
            if resolved is None:
                error.append(klass=self.TASKNAME, code=204,
                             details='Set PANDDA_EXECUTABLE, the PANDDA2_EXECUTABLE preference, '
                                     'or install a CCP4 that ships pandda2.analyse; '
                                     'or choose the stage-only run mode',
                             name=f'{self.TASKNAME}.container.controlParameters.PANDDA_EXECUTABLE',
                             severity=SEVERITY_ERROR)

        if n < contract.MIN_DATASETS:
            error.append(klass=self.TASKNAME, code=205,
                         details=f'{n} dataset(s); PanDDA needs {contract.MIN_DATASETS} to run',
                         name=f'{self.TASKNAME}.container.inputData.DATASETS',
                         severity=SEVERITY_WARNING)

        if self._mode() == 'local':
            scratch = self._scratch_dir()
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
        provenance = {
            'task': self.TASKNAME,
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
        self._write_program_xml(state='staged')
        return None

    def makeCommandAndScript(self):
        argv = contract.build_argv(self._staging_root / 'datasets', self._out_dir(),
                                   self._local_cpus())
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
        self._started = time.time()
        self._write_program_xml(state='running')
        stop = threading.Event()
        self._tailer_stop = stop
        tailer = threading.Thread(target=self._tail_progress, args=(stop,), daemon=True)
        tailer.start()
        try:
            return super().startProcess()
        finally:
            stop.set()
            tailer.join(timeout=5)

    def _prepareProcessExecution(self):
        prep = super()._prepareProcessExecution()
        prep['env'] = contract.build_env(prep['env'], self._scratch_dir())
        return prep

    def postProcessCheck(self, processId=None):
        status, exit_status, exit_code = super().postProcessCheck(processId)
        if status != CPluginScript.SUCCEEDED and self._mode() == 'local':
            text = self._read(self.makeFileName('STDERR')) + '\n' + self._read(self.makeFileName('LOG'))
            name, code, prompt = contract.classify_failure(text)
            self.appendErrorReport(code, f'{name}: {prompt}', stack=False)
            self._record_tree()
            self._write_program_xml(state='failed', failure=name)
        return status, exit_status, exit_code

    def processOutputFiles(self):
        if self._mode() == 'stage_only':
            return CPluginScript.SUCCEEDED
        processed, events, complete = self._record_tree()
        if processed == 0:
            self.appendErrorReport(220, f'nothing under {self._out_dir()}')
            self._write_program_xml(state='failed', failure='no_output')
            return CPluginScript.FAILED
        if not complete:
            self.appendErrorReport(221, f'{processed} processed dataset(s), no events table', stack=False)
            self._write_program_xml(state='partial')
            return CPluginScript.UNSATISFACTORY
        self._write_program_xml(state='finished')
        return CPluginScript.SUCCEEDED

    # -- helpers ------------------------------------------------------------

    def _mode(self) -> str:
        par = self.container.controlParameters
        return str(par.RUN_MODE) if par.RUN_MODE.isSet() else 'local'

    def _local_cpus(self) -> int:
        par = self.container.controlParameters
        return int(par.LOCAL_CPUS) if par.LOCAL_CPUS.isSet() else 4

    def _out_dir(self) -> Path:
        return Path(self.workDirectory) / contract.OUT_DIR_NAME

    def _scratch_dir(self) -> Path:
        par = self.container.controlParameters
        if par.SCRATCH_DIR.isSet() and str(par.SCRATCH_DIR).strip():
            return Path(str(par.SCRATCH_DIR)).expanduser()
        return Path(self.workDirectory) / 'ray_scratch'

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
            project_uuid = str(item.XYZIN.project) if item.XYZIN.project.isSet() else None
            specs.append(DatasetSpec(
                label=str(item.DTAG) if item.DTAG.isSet() else f'dataset-{len(specs) + 1}',
                xyzin=Path(str(item.XYZIN.fullPath)),
                hklin=Path(str(item.HKLIN.fullPath)),
                dictionary=Path(str(item.DICT.fullPath)) if item.DICT.isSet() else None,
                project_uuid=project_uuid,
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
        """Count what PanDDA wrote; set the output-tree fields. Returns
        ``(processed, events, complete)``."""
        out = self.container.outputData
        tree = self._out_dir()
        processed_dir = tree / 'processed_datasets'
        processed = len([d for d in processed_dir.iterdir() if d.is_dir()]) if processed_dir.is_dir() else 0
        table = tree / 'analyses' / 'pandda_analyse_events.csv'
        complete = table.is_file()
        events = 0
        if complete:
            events = max(0, sum(1 for _ in open(table)) - 1)
        if processed:
            out.PANDDA_OUT_DIR.set(str(tree))
        out.PERFORMANCE.nDatasetsProcessed.set(processed)
        out.PERFORMANCE.nEvents.set(events)
        if self._started is not None:
            out.PERFORMANCE.wallSeconds.set(round(time.time() - self._started, 1))
        return processed, events, complete

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
        n = len(self._manifest['datasets']) if self._manifest else len(self.container.inputData.DATASETS)
        ET.SubElement(root, 'n_datasets').text = str(n)
        if self._progress:
            ET.SubElement(root, 'progress', done=str(self._progress[0]), total=str(self._progress[1]))
        perf = self.container.outputData.PERFORMANCE
        for tag, field in (('n_processed', perf.nDatasetsProcessed), ('n_events', perf.nEvents),
                           ('wall_seconds', perf.wallSeconds)):
            if field.isSet():
                ET.SubElement(root, tag).text = str(field)
        if self._manifest:
            datasets = ET.SubElement(root, 'datasets')
            for entry in self._manifest['datasets']:
                ET.SubElement(datasets, 'dataset', xtal=entry['xtal'], label=entry['label'],
                              dict='yes' if 'dict.cif' in entry['files'] else 'no')
        tree = ET.ElementTree(root)
        ET.indent(tree, space='  ')
        tree.write(self.makeFileName('PROGRAMXML'), encoding='utf-8', xml_declaration=True)
