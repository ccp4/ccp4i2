"""Rigid-body refinement of a placed model with Phaser (MR_RNP), then
optionally shift-field refinement and refmac -- a pipeline that hosts
Phaser's PHIL.

The model is a parent structure cut into rigid bodies by atom selections.
Each becomes an ensemble already placed at the origin of its own
coordinates (Phaser's solution_at_origin), so there is no solution file to
hand over. Before Phaser the data are reindexed to match the parent, so
model and data agree on the indexing convention.

The parameters are the MR_RNP task's: the same master PHIL filtered to the
same mode, handed over whole. Resolution limits are Phaser's own to set,
where the classic pipeline hard-coded 25-3 A.
"""
import os

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.PhilPluginScript import PhilPluginScript
from ccp4i2.pipelines.phaser_pipeline_phil.script.phaser_pipeline_phil import phaser_pipeline_phil
from ccp4i2.wrappers.phaser_mr_rnp_phil.script.phaser_mr_rnp_phil import phaser_mr_rnp_phil


class phaser_rnp_pipeline_phil(phaser_pipeline_phil):

    TASKNAME = "phaser_rnp_pipeline_phil"
    MR_TASK = "phaser_mr_rnp_phil"
    PHIL_MODE = phaser_mr_rnp_phil.PHIL_MODE
    WHATNEXT = ["prosmart_refmac", "modelcraft", "coot_rebuild"]
    #: The identity a fragment of the parent model claims to the target: the
    #: classic pipeline's value, for a structure taken to be isomorphous
    FRAGMENT_IDENTITY = 0.9

    ERROR_CODES = {
        **phaser_pipeline_phil.ERROR_CODES,
        208: {"description": "Extracting a rigid body from the parent model failed"},
        210: {"description": "Selection yields no atoms"},
        212: {"description": "Overlapping atom selections"},
    }

    def get_phil_exclude_scopes(self):
        return list(phaser_mr_rnp_phil.PHIL_EXCLUDE_SCOPES) + phaser_mr_rnp_phil.phil_shim_targets()

    # -- validation -----------------------------------------------------------
    def validity(self):
        error = PhilPluginScript.validity(self)
        # The task's composition rule, on the same inputs; its ensemble rules
        # do not apply, the ensembles being built here at run time
        probe = phaser_mr_rnp_phil(workDirectory=str(self.getWorkDirectory()), name=self.MR_TASK)
        self._copy_inputs(probe)
        for report in probe.validity()._reports:
            if report["code"] == 113:
                error.append(klass=self.TASKNAME, code=113, details=report["details"],
                             name=str(report.get("name", "")).replace(f"{self.MR_TASK}.", f"{self.TASKNAME}."),
                             severity=report["severity"])
        self._check_selections(error)
        return error

    def _check_selections(self, error):
        """Every selection must pick some atoms of the parent, and no two may
        share any: a rigid body cannot be in two places."""
        inp = self.container.inputData
        if not inp.XYZIN_PARENT.isSet() or len(inp.SELECTIONS) == 0:
            return
        try:
            content = inp.XYZIN_PARENT.fileContent
        except Exception:
            content = None
        if content is None:
            return
        atom_sets = []
        for index, selection in enumerate(inp.SELECTIONS):
            text = str(selection.text)
            try:
                n_atoms, selected = content.interpretSelection(text)
            except Exception:
                n_atoms, selected = 0, []
            if n_atoms == 0:
                error.append(klass=self.TASKNAME, code=210,
                             details=f'Selection {index + 1} ("{text}") matches no atoms in the parent model',
                             name=f"{self.TASKNAME}.container.inputData.SELECTIONS[{index}]",
                             severity=CCP4ErrorHandling.SEVERITY_ERROR)
                continue
            atom_sets.append((index, text, {
                (chain.name, residue.seqid.num, residue.seqid.icode, residue.name, atom.name, atom.altloc)
                for _model, chain, residue, atom in selected}))
        for i in range(len(atom_sets)):
            for j in range(i + 1, len(atom_sets)):
                index_i, text_i, set_i = atom_sets[i]
                index_j, text_j, set_j = atom_sets[j]
                shared = len(set_i & set_j)
                if shared:
                    error.append(klass=self.TASKNAME, code=212,
                                 details=f'Selections {index_i + 1} ("{text_i}") and {index_j + 1} '
                                         f'("{text_j}") share {shared} atom(s)',
                                 name=f"{self.TASKNAME}.container.inputData.SELECTIONS",
                                 severity=CCP4ErrorHandling.SEVERITY_ERROR)

    # -- the run --------------------------------------------------------------
    def process(self):
        invalid = self.checkInputData()
        if invalid:
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        self.checkOutputData()
        self.xmlroot = etree.Element("PhaserPipeline")
        inp, out = self.container.inputData, self.container.outputData
        # The data indexed as the parent model is
        if self.runPointless() != CPluginScript.SUCCEEDED:
            return CPluginScript.FAILED
        f_sigf = out.F_SIGF_OUT
        freer = out.FREERFLAG_OUT if inp.FREERFLAG.isSet() else None
        self.phaserPlugin = self.makePluginObject(self.MR_TASK)
        self._copy_inputs(self.phaserPlugin)
        self.phaserPlugin.container.inputData.F_SIGF.set(f_sigf)
        if not self.createEnsembles(self.phaserPlugin):
            return CPluginScript.FAILED
        self.hand_phil_to(self.phaserPlugin)
        self.phaserPlugin.xml_responders.append(self.phaserXMLUpdated)
        rv = self.phaserPlugin.process()
        self.phaserXMLUpdated(self.phaserPlugin.xmlroot)
        if rv != CPluginScript.SUCCEEDED:
            self.errorReport.extend(self.phaserPlugin.getErrorReport())
            self.appendErrorReport(200, severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        if len(self.phaserPlugin.container.outputData.XYZOUT) == 0:
            self.appendErrorReport(207, severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        self.harvestPhaser()
        xyzin = out.XYZOUT[0]
        if inp.RUNSHEETBEND:
            if self.runSheetbend(f_sigf, freer, xyzin) != CPluginScript.SUCCEEDED:
                return CPluginScript.FAILED
            xyzin = out.XYZOUT_SHEETBEND
        if inp.RUNREFMAC:
            if self.runRefmac(f_sigf, freer, xyzin) != CPluginScript.SUCCEEDED:
                return CPluginScript.FAILED
        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    def createEnsembles(self, plugin):
        """One ensemble per selection, cut from the parent model and placed at
        the origin of its own coordinates. No selection: the whole parent as
        one rigid body."""
        inp, sub = self.container.inputData, plugin.container.inputData
        sub.ENSEMBLES.clear()
        sub.FIXENSEMBLES.clear()
        texts = [str(selection.text) for selection in inp.SELECTIONS
                 if str(selection.text).strip()] or [None]
        for index, text in enumerate(texts):
            label = f"Fragment_{index + 1}"
            path = os.path.join(str(self.workDirectory), f"{label}.pdb")
            try:
                fragment = CPdbDataFile()
                fragment.setFullPath(str(inp.XYZIN_PARENT.fullPath))
                if text:
                    fragment.selection.text.set(text)
                rc = fragment.getSelectedAtomsPdbFile(path)
                if rc != 0:
                    raise RuntimeError(f"writing the selected atoms returned {rc}")
            except Exception as err:
                self.appendErrorReport(208, f"{label} ({text or 'whole model'}): {err}",
                                       severity=CCP4ErrorHandling.SEVERITY_ERROR)
                self.reportStatus(CPluginScript.FAILED)
                return False
            sub.ENSEMBLES.append(sub.ENSEMBLES.makeItem())
            ensemble = sub.ENSEMBLES[-1]
            ensemble.label.set(label)
            ensemble.number.set(0)
            ensemble.use.set(False)
            item = ensemble.pdbItemList.makeItem()
            ensemble.pdbItemList.append(item)
            item.structure.setFullPath(path)
            item.identity_to_target.set(self.FRAGMENT_IDENTITY)
            sub.FIXENSEMBLES.append(label)
        return True

    def runPointless(self):
        """Reindex the data to match the parent model, and record it."""
        try:
            plugin = self.makePluginObject("pointless_reindexToMatch")
            plugin.container.controlParameters.REFERENCE = "XYZIN_REF"
            plugin.container.inputData.XYZIN_REF.set(self.container.inputData.XYZIN_PARENT)
            plugin.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF)
            if self.container.inputData.FREERFLAG.isSet():
                plugin.container.inputData.FREERFLAG.set(self.container.inputData.FREERFLAG)
            rv = plugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(205, severity=CCP4ErrorHandling.SEVERITY_ERROR)
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
            self.harvestFile(plugin.container.outputData.F_SIGF_OUT, self.container.outputData.F_SIGF_OUT)
            if self.container.inputData.FREERFLAG.isSet():
                self.harvestFile(plugin.container.outputData.FREERFLAG_OUT, self.container.outputData.FREERFLAG_OUT)
            self.appendXML(plugin.makeFileName("PROGRAMXML"), "POINTLESS")
        except Exception as err:
            self.appendErrorReport(205, str(err), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED
