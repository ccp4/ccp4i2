"""``pandda_fanout`` -- fan a PanDDA run out into per-dataset receipts.

Design note: docs/pandda-campaign-design.md, section 8 (and decision 1 as
amended on 2026-09-24: fan-out is a task after all).

Takes the manifest a PanDDA job wrote when it staged its run, finds the
output tree (that job's, or one named for a tree produced elsewhere), and
creates one ``pandda_events`` receipt in each dataset's own project, running
each to completion as it is created. The job is the record of what a
fan-out did: which tree, which manifest, and per dataset what happened. It
is not a parent of the receipts, and it is not a claim about how they fare
afterwards: a receipt's own job carries that.

Needs no CCP4: it copies files and creates jobs, which is why it is
``ccp4_free`` and can run wherever the API runs.
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from xml.etree import ElementTree as ET

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ErrorHandling import CErrorReport, SEVERITY_WARNING, SEVERITY_ERROR

logger = logging.getLogger(f"ccp4i2:{__name__}")


class pandda_fanout(CPluginScript):

    TASKNAME = 'pandda_fanout'
    TASKVERSION = 1.0
    TASKMODULE = 'ligands'
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS = False
    PERFORMANCECLASS = 'CPanddaFanoutPerformance'

    ERROR_CODES = {
        201: {'description': 'The manifest could not be read'},
        202: {'description': 'No PanDDA output tree to fan out'},
        203: {'severity': SEVERITY_WARNING, 'description': 'Some receipts could not be created or run'},
        204: {'severity': SEVERITY_WARNING, 'description': 'The run wrote nothing for some datasets'},
        205: {'severity': SEVERITY_WARNING, 'description': 'Some datasets belong to projects not in this database'},
        206: {'severity': SEVERITY_WARNING, 'description': 'Preview only: nothing was created'},
        207: {'severity': SEVERITY_WARNING, 'description': 'Nothing to do: every dataset already has a receipt, or none is in the tree'},
    }

    def runTimeValidity(self) -> CErrorReport:
        error = super().runTimeValidity()
        if error.maxSeverity() >= SEVERITY_ERROR:
            return error
        from ccp4i2.lib.utils.jobs.pandda_fanout import read_manifest
        try:
            manifest = read_manifest(str(self.container.inputData.MANIFEST.fullPath))
        except (OSError, ValueError, json.JSONDecodeError) as e:
            error.append(klass=self.TASKNAME, code=201, details=str(e),
                         name=f'{self.TASKNAME}.container.inputData.MANIFEST', severity=SEVERITY_ERROR)
            return error
        tree = self._tree(manifest)
        if tree is None or not (tree / 'processed_datasets').is_dir():
            error.append(klass=self.TASKNAME, code=202,
                         details=(f'{tree or "(no tree)"} has no processed_datasets/. The manifest names no run '
                                  'job in this database; set PANDDA_OUT_DIR to the tree'),
                         name=f'{self.TASKNAME}.container.controlParameters.PANDDA_OUT_DIR',
                         severity=SEVERITY_ERROR)
        return error

    def makeCommandAndScript(self):
        return None

    def startProcess(self):
        from ccp4i2.lib.utils.jobs.pandda_fanout import execute_fanout, plan_fanout, read_manifest
        par = self.container.controlParameters
        try:
            manifest = read_manifest(str(self.container.inputData.MANIFEST.fullPath))
        except (OSError, ValueError, json.JSONDecodeError) as e:
            self.appendErrorReport(201, str(e))
            return CPluginScript.FAILED
        tree = self._tree(manifest)
        if tree is None or not (tree / 'processed_datasets').is_dir():
            self.appendErrorReport(202, str(tree))
            return CPluginScript.FAILED
        run_job_uuid = (manifest.get('provenance') or {}).get('run_job_uuid')
        dry_run = bool(par.DRY_RUN)
        run_receipts = bool(par.RUN_RECEIPTS) if par.RUN_RECEIPTS.isSet() else True

        plan = plan_fanout(tree, manifest, run_job_uuid)
        if not dry_run:
            execute_fanout(plan, run=run_receipts)

        counts = {a: plan.count(a) for a in ('created', 'skipped', 'absent', 'failed', 'no_project')}
        perf = self.container.outputData.PERFORMANCE
        perf.nCreated.set(counts['created'])
        perf.nSkipped.set(counts['skipped'])
        perf.nAbsent.set(counts['absent'])
        perf.nFailed.set(counts['failed'])
        perf.nNoProject.set(counts['no_project'])
        self._write_program_xml(plan, dry_run, run_receipts)

        names = lambda action: ', '.join(o.xtal for o in plan.outcomes if o.action == action)
        if dry_run:
            self.appendErrorReport(206, f"would create {counts['created']}", stack=False)
        if counts['absent']:
            self.appendErrorReport(204, names('absent'), stack=False)
        if counts['no_project']:
            self.appendErrorReport(205, names('no_project'), stack=False)
        if counts['failed']:
            self.appendErrorReport(203, '; '.join(f'{o.xtal}: {o.reason}' for o in plan.outcomes
                                                   if o.action == 'failed'), stack=False)
            return CPluginScript.UNSATISFACTORY
        if not dry_run and counts['created'] == 0 and counts['skipped'] == 0:
            self.appendErrorReport(207, '', stack=False)
            return CPluginScript.UNSATISFACTORY
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        return CPluginScript.SUCCEEDED

    def _tree(self, manifest):
        par = self.container.controlParameters
        if par.PANDDA_OUT_DIR.isSet() and str(par.PANDDA_OUT_DIR).strip():
            return Path(str(par.PANDDA_OUT_DIR)).expanduser()
        run_job_uuid = (manifest.get('provenance') or {}).get('run_job_uuid')
        if not run_job_uuid:
            return None
        from ccp4i2.db import models
        job = models.Job.objects.filter(uuid=run_job_uuid).first()
        if job is None:
            return None
        return Path(job.directory) / 'pandda2_out'

    def _write_program_xml(self, plan, dry_run, run_receipts):
        root = ET.Element('pandda_fanout')
        ET.SubElement(root, 'tree').text = str(plan.tree)
        ET.SubElement(root, 'run_job_uuid').text = plan.run_job_uuid or ''
        ET.SubElement(root, 'incomplete').text = str(plan.incomplete)
        ET.SubElement(root, 'dry_run').text = str(dry_run)
        ET.SubElement(root, 'run_receipts').text = str(run_receipts)
        for action in ('created', 'skipped', 'absent', 'failed', 'no_project'):
            ET.SubElement(root, f'n_{action}').text = str(plan.count(action))
        datasets = ET.SubElement(root, 'datasets')
        for o in plan.outcomes:
            ET.SubElement(datasets, 'dataset', xtal=o.xtal, label=o.label, action=o.action,
                          project=o.project or '', job_uuid=o.job_uuid or '').text = o.reason
        tree = ET.ElementTree(root)
        ET.indent(tree, space='  ')
        tree.write(self.makeFileName('PROGRAMXML'), encoding='utf-8', xml_declaration=True)
