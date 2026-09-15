"""``molrep_map`` -- fast cryo-EM map -> model placement, prepared for refinement.

Place a model into a cryo-EM map, resolving the hand ambiguity by running the MR
engine on *both* the map and its origin-inverted copy, and emit a refinement-ready
package -- a trimmed, Moorhen/clipper-registered map, a mask, and the placed model
-- for each hand. Half maps, if supplied, are carried through the recommended
hand's transform for servalcat's cross-validated ``--halfmaps`` path.

Design and crystallographic rationale: ``docs/molrep-map-design.md``. The map is
P1, which has no origin for an unphased translation function, so the engine
searches against the map's own *phased* density (molrep ``-f <map>``); the trimmed
output keeps the full cell and encodes position via ``nxstart`` (what coot/clipper
/Moorhen read). All map surgery is gemmi (``preprocess_map``); no ``chapi``.
"""

from __future__ import annotations

import os
import re
import shutil
import logging
from xml.etree import ElementTree as ET

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING

logger = logging.getLogger(__name__)


class molrep_map(CPluginScript):

    TASKNAME = 'molrep_map'
    TASKVERSION = 1.0
    TASKMODULE = 'molecular_replacement'
    WHATNEXT = ['servalcat_pipe']
    ASYNCHRONOUS = False
    RUNEXTERNALPROCESS = False

    ERROR_CODES = {
        201: {'description': 'Failed to read or prepare the input map'},
        202: {'description': 'Failed to prepare the input model'},
        203: {'severity': SEVERITY_WARNING, 'description': 'Placement produced no model for a hand'},
        204: {'severity': SEVERITY_WARNING, 'description': 'Failed to trim/mask a hand map'},
        205: {'severity': SEVERITY_WARNING, 'description': 'Failed to prepare half maps'},
        206: {'severity': SEVERITY_WARNING, 'description': 'Failed to write program XML'},
        207: {'description': 'The selected MR engine is not available'},
    }

    # Hands: (label, output-model attr, output-map attr, output-mask attr)
    HANDS = (
        ('Original', 'ORIGINALMODEL', 'ORIGINALTRIMMEDMAP', 'ORIGINALMASK'),
        ('Flipped', 'FLIPPEDMODEL', 'FLIPPEDTRIMMEDMAP', 'FLIPPEDMASK'),
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._results = {}          # hand label -> engines.PlacementResult
        self._boxes = {}            # hand label -> gemmi.FractionalBox
        self._recommended = None

    # ---- pipeline hooks -------------------------------------------------

    def processInputFiles(self):
        from . import preprocess_map as pp

        self.workDir = str(self.getWorkDirectory())
        inp = self.container.inputData
        par = self.container.controlParameters

        try:
            self._full = {
                'Original': pp.read_map(str(inp.MAPIN.fullPath)),
            }
            self._full['Flipped'] = pp.flip_hand(self._full['Original'])
        except Exception as e:
            self.appendErrorReport(201, str(e))
            return CPluginScript.FAILED

        # Disposable, down-sampled + blurred search maps -- one per hand.
        try:
            downs = float(par.DOWNSAMPLE) if par.DOWNSAMPLE.isSet() else 1.0
            badd = float(par.BADD) if par.BADD.isSet() else 0.0
            self._search = {}
            for hand in ('Original', 'Flipped'):
                search = pp.prepare_for_search(self._full[hand], downs, badd)
                path = os.path.join(self.workDir, f'search_{hand}.map')
                pp.write_map(search, path)
                self._search[hand] = path
        except Exception as e:
            self.appendErrorReport(201, str(e))
            return CPluginScript.FAILED

        # Ensure the model is a PDB molrep can read.
        try:
            self._model = self._ensure_pdb(inp.XYZIN)
        except Exception as e:
            self.appendErrorReport(202, str(e))
            return CPluginScript.FAILED

        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        # This task drives the engine itself in startProcess(); nothing for the
        # single-command framework to build.
        return CPluginScript.SUCCEEDED

    def startProcess(self):
        from . import engines
        from . import preprocess_map as pp

        par = self.container.controlParameters
        out = self.container.outputData

        engine = str(par.ENGINE) if par.ENGINE.isSet() else engines.MOLREP
        nmon = int(par.NMON) if par.NMON.isSet() else 1
        np_peaks = int(par.NP) if par.NP.isSet() else 5
        badd = float(par.BADD) if par.BADD.isSet() else 0.0
        resmax = float(par.SEARCH_RESOLUTION) if par.SEARCH_RESOLUTION.isSet() else None
        tlim = float(par.TIME_LIMIT) if par.TIME_LIMIT.isSet() else None
        border = float(par.BORDER) if par.BORDER.isSet() else 5.0
        mask_radius = float(par.MASK_RADIUS) if par.MASK_RADIUS.isSet() else 3.0

        for hand, model_attr, map_attr, mask_attr in self.HANDS:
            hand_dir = os.path.join(self.workDir, hand)
            try:
                res = engines.place(
                    engine, self._search[hand], self._model, hand_dir,
                    nmon=nmon, np_peaks=np_peaks, b_add=badd,
                    resmax=resmax, time_limit=tlim)
            except NotImplementedError as e:
                self.appendErrorReport(207, str(e))
                return CPluginScript.FAILED
            self._results[hand] = res

            if not res.placed:
                self.appendErrorReport(203, f'{hand} hand', stack=False)
                continue

            try:
                self._write_model(res.model_path, str(getattr(out, model_attr).fullPath))
                cell = self._full[hand].grid.unit_cell
                box = pp.model_frac_box(res.model_path, cell, border)
                self._boxes[hand] = box
                # Mask is built on the full grid (correct sampling) then cropped
                # to the same box, so it aligns with the trimmed map voxel-for-voxel.
                mask = pp.atom_mask(self._full[hand], res.model_path, mask_radius, box)
                pp.write_map(mask, str(getattr(out, mask_attr).fullPath))
                density = pp.crop(self._full[hand], box)
                pp.write_map(density, str(getattr(out, map_attr).fullPath))
            except Exception as e:
                self.appendErrorReport(204, f'{hand}: {e}')

        self._recommended = self._choose_hand()
        self._prepare_half_maps(pp)
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        out = self.container.outputData

        # Annotate whichever outputs actually landed on disk (the gleaner only
        # persists set output files that exist, so unproduced ones drop quietly).
        for hand, model_attr, map_attr, mask_attr in self.HANDS:
            tag = 'original-hand' if hand == 'Original' else 'inverted-hand'
            rec = ' (recommended)' if hand == self._recommended else ''
            if os.path.exists(str(getattr(out, model_attr).fullPath)):
                getattr(out, model_attr).annotation = f'Model placed in the {tag} map{rec}'
            if os.path.exists(str(getattr(out, map_attr).fullPath)):
                getattr(out, map_attr).annotation = f'{tag.capitalize()} map trimmed to the model{rec}'
            if os.path.exists(str(getattr(out, mask_attr).fullPath)):
                getattr(out, mask_attr).annotation = f'Mask around the {tag} model{rec}'

        for attr, n in (('HALFMAPOUT1', 1), ('HALFMAPOUT2', 2)):
            if os.path.exists(str(getattr(out, attr).fullPath)):
                getattr(out, attr).annotation = (
                    f'Half map {n} trimmed to the {self._recommended.lower()} hand')

        try:
            self._write_program_xml()
        except Exception as e:
            self.appendErrorReport(206, str(e))

        return CPluginScript.SUCCEEDED

    # ---- helpers --------------------------------------------------------

    def _ensure_pdb(self, xyzin):
        """Return a filesystem path to a PDB form of the input model."""
        if not xyzin.isSet():
            raise ValueError('XYZIN is not set')
        try:
            ext = xyzin.getExt()
        except Exception:
            ext = os.path.splitext(str(xyzin.fullPath))[1]
        if ext == '.pdb':
            return str(xyzin.fullPath)
        converted = os.path.join(self.workDir, 'model_input.pdb')
        xyzin.convertFormat('pdb', converted)
        return converted

    def _write_model(self, src, dst):
        """Copy molrep's placed model to ``dst``, stripping its #MOLECULE tags."""
        with open(src) as istream:
            content = istream.read()
        content = re.sub(r'\n#MOLECULE\s+[0-9]+\s*', '\n', content)
        with open(dst, 'w') as ostream:
            ostream.write(content)

    def _choose_hand(self):
        """The higher-scoring hand, defaulting to Original when scores are absent."""
        scores = {h: self._results[h].score
                  for h in ('Original', 'Flipped') if h in self._results}
        placed = {h: self._results[h].placed
                  for h in ('Original', 'Flipped') if h in self._results}
        # Prefer a hand that actually placed.
        candidates = [h for h in ('Original', 'Flipped') if placed.get(h)]
        if not candidates:
            return 'Original'
        if all(scores.get(h) is None for h in candidates):
            return candidates[0]
        return max(candidates,
                   key=lambda h: scores[h] if scores.get(h) is not None else float('-inf'))

    def _prepare_half_maps(self, pp):
        """Carry the recommended hand's flip+trim onto the half maps, if given."""
        inp = self.container.inputData
        out = self.container.outputData
        if not (inp.HALFMAP1.isSet() and inp.HALFMAP2.isSet()):
            return
        if self._recommended not in self._boxes:
            return
        box = self._boxes[self._recommended]
        flip = self._recommended == 'Flipped'
        try:
            for hm_in, hm_out in ((inp.HALFMAP1, out.HALFMAPOUT1),
                                  (inp.HALFMAP2, out.HALFMAPOUT2)):
                hmap = pp.read_map(str(hm_in.fullPath))
                if flip:
                    hmap = pp.flip_hand(hmap)
                hmap = pp.crop(hmap, box)
                pp.write_map(hmap, str(hm_out.fullPath))
        except Exception as e:
            self.appendErrorReport(205, str(e))

    # ---- report XML -----------------------------------------------------

    def _write_program_xml(self):
        root = ET.Element('molrep_map')
        rec = ET.SubElement(root, 'recommendation')
        rec.set('hand', self._recommended or 'Original')
        for hand in ('Original', 'Flipped'):
            res = self._results.get(hand)
            el = ET.SubElement(root, hand)
            if res is None:
                el.set('placed', 'false')
                continue
            el.set('placed', 'true' if res.placed else 'false')
            el.set('timed_out', 'true' if res.timed_out else 'false')
            if res.score is not None:
                el.set('score', f'{res.score:.4f}')
            if res.doc_path:
                el.append(self._scrape_doc(res.doc_path))
        with open(str(self.makeFileName('PROGRAMXML')), 'w') as fh:
            ET.indent(root)
            fh.write(ET.tostring(root).decode('utf-8'))

    @staticmethod
    def _scrape_doc(doc_path):
        """Scrape molrep.doc's (phased) TF peak table into a MolrepResult element."""
        results = ET.Element('MolrepResult')
        titles = []
        rf = None
        in_tf = False
        try:
            with open(doc_path) as fh:
                lines = fh.read().split('\n')
        except OSError:
            return results
        for line in lines:
            stripped = line.strip()
            if stripped in ('--- Translation function ---',
                            '--- phased translation function ---'):
                in_tf = True
                continue
            if not in_tf:
                continue
            if stripped.startswith('RF '):
                titles = (line.replace('(', ' ').replace(')', '')
                          .replace('/', '_').split())
                rf = ET.SubElement(results, 'RFpeaks')
                continue
            words = (line.replace('(', ' ').replace(')', '')
                     .replace('-', ' -').split())
            if titles and rf is not None and len(words) == len(titles):
                try:
                    int(words[0]); int(words[1])
                    for i in range(2, len(words)):
                        float(words[i])
                except (ValueError, IndexError):
                    continue
                peak = ET.SubElement(rf, 'RFpeak')
                for key, value in zip(titles, words):
                    child = ET.SubElement(peak, key)
                    child.text = str(float(value))
        return results
