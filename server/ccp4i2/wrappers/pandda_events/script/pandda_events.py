"""``pandda_events`` -- a receipt for one dataset's share of one PanDDA run.

Design note: docs/pandda-campaign-design.md, section 7.

One job per dataset per run, in the dataset's own project. It reads the
dataset's directory in a PanDDA 2 output tree, copies what it finds into its
own job directory (hardlink or copy, never a reference in place: a project
directory stays self-contained, section 8.3), and declares it as typed
outputs -- the apo model of record, the Z-map, and every event as a
``CPanddaEvent`` carrying its own map and pose.

Its entire logic budget is to check that what was delivered matches what
PanDDA declared, and to say so loudly when it does not (section 7.2). Loudly
means ``UNSATISFACTORY``, never ``FAILED``: everything that arrived is still
gleaned and usable, and the shortfall is in the report, the KPIs and the job
status. ``FAILED`` is for a tree with no directory for this dataset at all.

No external program, no database access, no CCP4 binary.
"""
from __future__ import annotations

import logging
import os
from pathlib import Path
from xml.etree import ElementTree as ET

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ErrorHandling import CErrorReport, SEVERITY_WARNING, SEVERITY_ERROR
from ccp4i2.core.CCP4XtalData import CMapDataFile
from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.wrappers.pandda_campaign.script.pandda_staging import link_or_copy

from .pandda_tree import DatasetNotFound, read_dataset

logger = logging.getLogger(f"ccp4i2:{__name__}")


class pandda_events(CPluginScript):

    TASKNAME = 'pandda_events'
    TASKVERSION = 1.0
    TASKMODULE = 'ligands'
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS = False
    PERFORMANCECLASS = 'CPanddaReceiptPerformance'

    ERROR_CODES = {
        201: {'description': 'No directory for this dataset in the PanDDA output tree'},
        202: {'severity': SEVERITY_WARNING,
              'description': 'Receipt is short: PanDDA declared outputs that are not in the tree'},
        203: {'severity': SEVERITY_WARNING, 'description': 'Could not copy a file from the tree'},
        204: {'severity': SEVERITY_WARNING,
              'description': 'This receipt comes from a run that did not finish'},
    }

    # -- validation ---------------------------------------------------------

    def runTimeValidity(self) -> CErrorReport:
        error = super().runTimeValidity()
        if error.maxSeverity() >= SEVERITY_ERROR:
            return error
        tree = str(self.container.inputData.PANDDA_OUT_DIR)
        if not os.path.isdir(tree):
            error.append(klass=self.TASKNAME, code=201,
                         details=f'Not a directory: {tree}',
                         name=f'{self.TASKNAME}.container.inputData.PANDDA_OUT_DIR',
                         severity=SEVERITY_ERROR)
        return error

    # -- the whole task ----------------------------------------------------

    def makeCommandAndScript(self):
        return None

    def startProcess(self):
        inp = self.container.inputData
        out = self.container.outputData
        tree = Path(str(inp.PANDDA_OUT_DIR))
        dtag = str(inp.DTAG)

        try:
            dataset = read_dataset(tree, dtag)
        except DatasetNotFound as e:
            self.appendErrorReport(201, str(e))
            return CPluginScript.FAILED

        self._take(dataset.apo_model, out.XYZIN_APO)
        self._take(dataset.zmap, out.ZMAP)
        self._take(dataset.pandda_model, out.PANDDA_MODEL)

        best_score = None
        for event in dataset.events:
            item = out.EVENTS.makeItem()
            item.EVENT_IDX.set(event.idx)
            self._set_float(item.BDC, event.bdc)
            self._set_float(item.SCORE, event.score)
            self._set_float(item.HIT_PROBABILITY, event.hit_probability)
            if event.site_idx is not None:
                item.SITE_IDX.set(event.site_idx)
            if event.centroid is not None:
                item.CENTROID.x.set(event.centroid[0])
                item.CENTROID.y.set(event.centroid[1])
                item.CENTROID.z.set(event.centroid[2])
            if event.build is not None:
                self._set_float(item.BUILD_SCORE, event.build.build_score)
                self._set_float(item.RSCC, event.build.rscc)
                self._set_float(item.OPTIMAL_CONTOUR, event.build.optimal_contour)
                if event.build.build_score is not None and event.build.path is not None:
                    best_score = (event.build.build_score if best_score is None
                                  else max(best_score, event.build.build_score))
            if event.event_map is not None:
                dst = os.path.join(self.workDirectory, f'event_{event.idx}_map.map')
                if self._copy(event.event_map, dst):
                    item.EVENT_MAP.setFullPath(dst)
                    item.EVENT_MAP.subType.set(CMapDataFile.SUBTYPE_NORMAL)
                    item.EVENT_MAP.annotation.set(
                        f'{dtag} event {event.idx} map (BDC {self._fmt(event.bdc)})')
            if event.pose is not None:
                dst = os.path.join(self.workDirectory, f'event_{event.idx}_pose.pdb')
                if self._copy(event.pose, dst):
                    item.POSE.setFullPath(dst)
                    item.POSE.subType.set(CPdbDataFile.SUBTYPE_FRAGMENT)
                    item.POSE.annotation.set(
                        f'{dtag} event {event.idx} candidate pose '
                        f'(build score {self._fmt(event.build.build_score if event.build else None)})')
            out.EVENTS.append(item)

        perf = out.PERFORMANCE
        perf.nEventsExpected.set(dataset.n_events)
        perf.nEventsDelivered.set(dataset.n_event_maps)
        perf.nPosesExpected.set(dataset.n_poses_expected)
        perf.nPosesDelivered.set(dataset.n_poses)
        if best_score is not None:
            perf.bestBuildScore.set(best_score)

        shortfalls = dataset.shortfalls()
        incomplete = bool(inp.RUN_INCOMPLETE)
        if incomplete:
            self.appendErrorReport(204, 'RUN_INCOMPLETE was set by fan-out', stack=False)
        if shortfalls:
            self.appendErrorReport(
                202, f'{len(shortfalls)} missing: ' + '; '.join(shortfalls), stack=False)

        self._write_program_xml(dataset, shortfalls, incomplete)

        # UNSATISFACTORY returns straight to reportStatus, and track_job
        # gleans that status: what arrived is published, and the job is
        # unmistakably not clean.
        if shortfalls:
            return CPluginScript.UNSATISFACTORY
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        return CPluginScript.SUCCEEDED

    # -- helpers ------------------------------------------------------------

    def _take(self, src, target):
        """Copy ``src`` (if it exists) to the path ``checkOutputData`` gave
        ``target``. Unset outputs are dropped by the gleaner."""
        if src is None:
            return
        self._copy(src, str(target.fullPath))

    def _copy(self, src, dst) -> bool:
        try:
            link_or_copy(src, dst)
            return True
        except OSError as e:
            self.appendErrorReport(203, f'{src} -> {dst}: {e}', stack=False)
            return False

    @staticmethod
    def _set_float(field, value):
        if value is not None:
            field.set(float(value))

    @staticmethod
    def _fmt(value):
        return '?' if value is None else f'{value:.2f}'

    def _write_program_xml(self, dataset, shortfalls, incomplete):
        root = ET.Element('pandda_events')
        ET.SubElement(root, 'dtag').text = dataset.dtag
        ET.SubElement(root, 'tree').text = str(dataset.directory.parent.parent)
        ET.SubElement(root, 'run_incomplete').text = str(bool(incomplete))
        ET.SubElement(root, 'events_table_present').text = str(dataset.events_table_present)
        for tag, present in (('apo_model', dataset.apo_model), ('zmap', dataset.zmap),
                             ('pandda_model', dataset.pandda_model)):
            ET.SubElement(root, tag).text = 'present' if present else 'absent'
        counts = ET.SubElement(root, 'counts')
        ET.SubElement(counts, 'events_expected').text = str(dataset.n_events)
        ET.SubElement(counts, 'event_maps_delivered').text = str(dataset.n_event_maps)
        ET.SubElement(counts, 'poses_expected').text = str(dataset.n_poses_expected)
        ET.SubElement(counts, 'poses_delivered').text = str(dataset.n_poses)
        missing = ET.SubElement(root, 'shortfalls')
        for text in shortfalls:
            ET.SubElement(missing, 'item').text = text
        events = ET.SubElement(root, 'events')
        for event in dataset.events:
            node = ET.SubElement(events, 'event', idx=str(event.idx))
            for tag, value in (('bdc', event.bdc), ('score', event.score),
                               ('hit_probability', event.hit_probability),
                               ('site_idx', event.site_idx)):
                if value is not None:
                    ET.SubElement(node, tag).text = str(value)
            if event.centroid is not None:
                ET.SubElement(node, 'centroid').text = ' '.join(f'{c:.2f}' for c in event.centroid)
            if event.build is not None:
                for tag, value in (('build_score', event.build.build_score),
                                   ('rscc', event.build.rscc),
                                   ('optimal_contour', event.build.optimal_contour)):
                    if value is not None:
                        ET.SubElement(node, tag).text = str(value)
            ET.SubElement(node, 'event_map').text = 'present' if event.event_map else 'absent'
            ET.SubElement(node, 'pose').text = (
                'present' if event.pose else ('absent' if event.build else 'none built'))
        tree = ET.ElementTree(root)
        ET.indent(tree, space='  ')
        tree.write(self.makeFileName('PROGRAMXML'), encoding='utf-8', xml_declaration=True)
