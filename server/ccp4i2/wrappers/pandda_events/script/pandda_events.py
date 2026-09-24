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

from .pandda_tree import PANDDA_RESIDUE_NAME, DatasetNotFound, read_dataset

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
        # PanDDA writes its maps as P1 in the crystal's cell. A P1 map with all
        # angles 90 degrees is what coot takes for an EM map, and it then draws
        # it as a box at its origin instead of a periodic, symmetry-expanded
        # crystal map around the model. The apo model still carries the
        # crystal's space group, so every map copy gets it back: metadata
        # PanDDA dropped, restored. A boxed event map is read into the cell
        # through the same fold, so it too is placed by its grid start.
        spacegroup = self._crystal_spacegroup(dataset.apo_model)
        self._take(dataset.zmap, out.ZMAP, spacegroup=spacegroup)
        # PanDDA names every residue it builds LIG, whatever the dictionary
        # said (autobuild/inbuilt.py). The copies that become CCP4i2 data get
        # the true component code back, so a pose can meet its dictionary in
        # refinement; PanDDA's own tree is left as it is.
        ligand_id = dataset.ligand_id
        rename = ligand_id if ligand_id and ligand_id != PANDDA_RESIDUE_NAME else None
        self._take(dataset.pandda_model, out.PANDDA_MODEL, rename=rename)
        # The dictionary the poses were built with, so a viewer draws them with
        # their bond orders and a refinement can take them: the job decides
        # which dictionary goes with its molecules, and this job carries it.
        self._take(dataset.dictionary, out.DICT)
        if dataset.dictionary is not None and ligand_id:
            out.DICT.annotation.set(f'Dictionary for {ligand_id}, as PanDDA used it')
        if rename:
            out.PANDDA_MODEL.annotation.set(
                f"PanDDA's merged model, ligand {ligand_id}: a machine opinion, not the model of record")

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
                if self._copy(event.event_map, dst, spacegroup=spacegroup):
                    item.EVENT_MAP.setFullPath(dst)
                    item.EVENT_MAP.subType.set(CMapDataFile.SUBTYPE_NORMAL)
                    item.EVENT_MAP.annotation.set(
                        f'{dtag} event {event.idx} map (BDC {self._fmt(event.bdc)})')
            if ligand_id:
                item.LIGAND_ID.set(ligand_id)
            if event.pose is not None:
                dst = os.path.join(self.workDirectory, f'event_{event.idx}_pose.pdb')
                if self._copy(event.pose, dst, rename=rename):
                    item.POSE.setFullPath(dst)
                    item.POSE.subType.set(CPdbDataFile.SUBTYPE_FRAGMENT)
                    item.POSE.annotation.set(
                        f'{dtag} event {event.idx} candidate pose'
                        f'{" of " + ligand_id if ligand_id else ""} '
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

    def _take(self, src, target, rename=None, spacegroup=None):
        """Copy ``src`` (if it exists) to the path ``checkOutputData`` gave
        ``target``. Unset outputs are dropped by the gleaner."""
        if src is None:
            return
        self._copy(src, str(target.fullPath), rename=rename, spacegroup=spacegroup)

    def _copy(self, src, dst, rename=None, spacegroup=None) -> bool:
        """Bytes by hardlink-or-copy; with ``rename`` a model is rewritten,
        with ``spacegroup`` a full-cell P1 map is rewritten with that space
        group -- real copies, never links, since a link would change the file
        inside PanDDA's tree too."""
        try:
            if rename:
                self._write_renamed(src, dst, rename)
            elif spacegroup is not None and self._write_with_spacegroup(src, dst, spacegroup):
                pass
            else:
                link_or_copy(src, dst)
            return True
        except (OSError, RuntimeError, ValueError) as e:
            self.appendErrorReport(203, f'{src} -> {dst}: {e}', stack=False)
            return False

    @staticmethod
    def _crystal_spacegroup(apo_model):
        """The crystal's space group from the apo model's CRYST1, or None
        when there is none or it is P1 (nothing to restore)."""
        if apo_model is None:
            return None
        import gemmi
        try:
            structure = gemmi.read_structure(str(apo_model))
        except Exception:      # noqa: BLE001 - no CRYST1 is no space group
            return None
        sg = structure.find_spacegroup()
        if sg is None or sg.number == 1:
            return None
        return sg

    @staticmethod
    def _write_with_spacegroup(src, dst, spacegroup) -> bool:
        """Write ``src`` to ``dst`` carrying ``spacegroup``, if it is a P1
        map (full cell or a box within it; the cell and sampling are the
        crystal's either way). Returns False when it already has a space
        group, and nothing was written."""
        import gemmi
        m = gemmi.read_ccp4_map(str(src))
        current = m.grid.spacegroup
        if current is not None and current.number != 1:
            return False
        # The header as read is complete and valid; only the ISPG word
        # changes. (update_ccp4_header would want a full-cell setup() first
        # for a boxed map, expanding it to the whole cell for nothing.)
        m.grid.spacegroup = spacegroup
        m.set_header_i32(23, spacegroup.ccp4)
        m.write_ccp4_map(str(dst))
        return True

    @staticmethod
    def _write_renamed(src, dst, ligand_id):
        import gemmi
        structure = gemmi.read_structure(str(src))
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.name == PANDDA_RESIDUE_NAME:
                        residue.name = ligand_id
        structure.setup_entities()
        structure.write_pdb(str(dst))

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
        ET.SubElement(root, 'ligand_id').text = dataset.ligand_id or ''
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
