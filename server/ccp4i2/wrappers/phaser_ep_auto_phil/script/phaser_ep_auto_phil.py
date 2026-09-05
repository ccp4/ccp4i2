"""SAD experimental phasing with Phaser, driven from its PHIL.

EP_AUTO over Phaser's own parameters (see phaser_mr_auto_phil for the
pattern). One crystal, one wavelength -- all the driver supports -- with
the anomalous data, the substructure and the wavelength written as a
phaser.crystal block; log-likelihood-gradient completion from the elements
asked for; optionally a partial model. Phaser runs in-process with the
recorder as its callback; afterwards the hands, sites and figures of merit
come from the ResultEP object and the completion cycles from Phaser's own
summary.
"""
import os

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.utils.phil_shims import FixedPhilShim
from ccp4i2.wrappers.phaser_phil.script.phaser_phil import phaser_phil
from ccp4i2.wrappers.phaser_phil.script.phaser_shims import (
    CompositionShim, EpCrystalShim, LlgCompletionShim, PartialModelShim)
from ccp4i2.wrappers.phaser_phil.script import phaser_run


class phaser_ep_auto_phil(phaser_phil):

    TASKNAME = "phaser_ep_auto_phil"
    TASKCOMMAND = None
    PHIL_MODE = "EP_AUTO"
    PHIL_MODE_PATH = "phaser.mode"
    WHATNEXT = ["parrot", "modelcraft", "coot_rebuild"]

    PHIL_EXCLUDE_SCOPES = [
        "phaser.sad_mode", "phaser.run_control", "phaser.output_dir",
        # The driver's own switches sit at the top level of its PHIL
        "dry_run", "show_script", "show_defaults", "test_mode", "verbose", "phaser.hklin", "phaser.labin",
        "phaser.keywords.partial_labin", "phaser.keywords.cluster_pdb",
    ]

    ERROR_CODES = {
        201: {"description": "Expected output file missing"},
        202: {"description": "Phaser reported an error"},
        204: {"description": "Failed to prepare the reflection data"},
    }

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.xmlroot = etree.Element("PhaserEpResults")
        self.resultObject = None

    @property
    def _crystal_shim(self):
        # Made on first use: the base constructor asks for the shims before
        # __init__ here runs
        if getattr(self, "_crystal_shim_instance", None) is None:
            self._crystal_shim_instance = EpCrystalShim(self)
        return self._crystal_shim_instance

    def validity(self):
        error = super().validity()
        inp = self.container.inputData
        name = f"{self.TASKNAME}.container.inputData"
        if not inp.XYZIN_HA.isSet() and not (
                str(inp.PARTIAL_BY) == "MODEL" and inp.XYZIN_PARTIAL.isSet()):
            error.append(klass=self.TASKNAME, code=113,
                         details="Give heavy-atom sites to start from, or a partial model.",
                         name=f"{name}.XYZIN_HA", severity=CCP4ErrorHandling.SEVERITY_ERROR)
        if not inp.WAVELENGTH.isSet():
            error.append(klass=self.TASKNAME, code=114,
                         details="The wavelength is needed for the anomalous scattering factors.",
                         name=f"{name}.WAVELENGTH", severity=CCP4ErrorHandling.SEVERITY_ERROR)
        if str(inp.COMP_BY) == "ASU" and not inp.ASUFILE.isSet():
            error.append(klass=self.TASKNAME, code=113,
                         details="Composition by AsuContent file: choose the file.",
                         name=f"{name}.ASUFILE", severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error

    def get_shim_definitions(self):
        root = os.path.join(str(self.getWorkDirectory()), "PHASER")
        return [
            self._crystal_shim,
            PartialModelShim(),
            LlgCompletionShim(),
            CompositionShim(),
            FixedPhilShim({"phaser.keywords.general.root": root}),
        ]

    def processInputFiles(self):
        error = self._crystal_shim.prepare()
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            self.appendErrorReport(204, str(error), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        self._phil_path = self.build_working_phil()
        return CPluginScript.SUCCEEDED

    def startProcess(self):
        recorder = phaser_run.PhaserRecorder(self.xmlroot, flush=self.flushXML)
        try:
            self.resultObject = self._run(recorder)
        except phaser_run.PhaserInputError as err:
            # Phaser refused the input; its sentence is the whole story
            self.appendErrorReport(202, str(err), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.flushXML(self.xmlroot)
            return CPluginScript.FAILED
        if not self.resultObject.Success():
            self.appendErrorReport(
                202, f"{self.resultObject.ErrorName()}: {self.resultObject.ErrorMessage()}",
                severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.flushXML(self.xmlroot)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def _run(self, recorder):
        return phaser_run.run_mode(
            self.get_master_phil(), self._phil_path, self.PHIL_MODE,
            str(self.getWorkDirectory()), recorder, self.makeFileName("LOG"))

    def processOutputFiles(self):
        result = self.resultObject
        out = self.container.outputData
        work_dir = str(self.getWorkDirectory())
        hands = ["", ".hand"] if getattr(result, "second_hand", False) else [""]
        for i, hand in enumerate(hands):
            xyz = os.path.join(work_dir, f"PHASER.1{hand}.pdb")
            hkl = os.path.join(work_dir, f"PHASER.1{hand}.mtz")
            for path in (xyz, hkl):
                if not os.path.exists(path):
                    self.appendErrorReport(201, path, severity=CCP4ErrorHandling.SEVERITY_ERROR)
                    return CPluginScript.FAILED
            label = "original hand" if i == 0 else "reversed hand"
            out.XYZOUT.append(out.XYZOUT.makeItem())
            out.XYZOUT[-1].setFullPath(xyz)
            out.XYZOUT[-1].annotation.set(f"Substructure sites - {label}")
            out.HKLOUT.append(out.HKLOUT.makeItem())
            out.HKLOUT[-1].setFullPath(hkl)
        if len(out.HKLOUT) > 0:
            from ccp4i2.core import CCP4XtalData
            self.splitHkloutList(
                miniMtzsOut=["ABCDOUT", "MAPOUT"],
                programColumnNames=["HLA,HLB,HLC,HLD", "FWT,PHWT"],
                outputBaseName=["ABCDOUT", "MAPOUT"],
                outputContentFlags=[CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL, 1],
                infileList=out.HKLOUT)
            for i in range(len(out.ABCDOUT)):
                label = "original hand" if i == 0 else "reversed hand"
                out.ABCDOUT[i].annotation.set(f"Phase estimates - {label}")
                out.MAPOUT[i].annotation.set(f"Phased map - {label}")
                out.MAPOUT[i].contentFlag.set(1)
                out.MAPOUT[i].subType.set(1)

        phaser_run.ep_results_xml(result, self.xmlroot)
        blocks = phaser_run.summary_blocks(result.summary())
        summaries = etree.SubElement(self.xmlroot, "Summaries")
        for name, text in blocks:
            node = etree.SubElement(summaries, "Summary")
            node.set("module", name)
            node.text = text
        cycles, final = phaser_run.ep_cycles(blocks)
        phaser_run.ep_cycles_xml(cycles, final, self.xmlroot)
        warnings = self.xmlroot.find("PhaserWarnings")
        if warnings is None:
            warnings = etree.SubElement(self.xmlroot, "PhaserWarnings")
        seen = {w.text for w in warnings}
        for text in result.warnings():
            if text not in seen:
                etree.SubElement(warnings, "Warning").text = text
        self.flushXML(self.xmlroot)
        return CPluginScript.SUCCEEDED

    def flushXML(self, xmlroot):
        target = self.makeFileName("PROGRAMXML")
        tmp = target + "_tmp"
        etree.ElementTree(xmlroot).write(tmp, pretty_print=True, xml_declaration=True,
                                         encoding="utf-8")
        os.replace(tmp, target)
