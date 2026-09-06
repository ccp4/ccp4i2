"""Molecular replacement with Phaser, then optionally shift-field
refinement and refmac -- a pipeline that hosts Phaser's PHIL.

The pipeline's parameters *are* the MR_AUTO task's: the same master PHIL,
filtered to the same mode, so the tree it shows is the tree the sub-job
runs with, handed over whole (hand_phil_to). Its typed inputs are the
task's plus a Free-R set, a reference structure and two switches. It keeps
no keyword snapshot of its own, and the sub-job's live record is embedded
in the pipeline's program.xml as it is written.
"""
import os
import shutil
from copy import deepcopy
from pathlib import Path

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.PhilPluginScript import PhilPluginScript
from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil import phaser_mr_auto_phil
from ccp4i2.wrappers.phaser_phil.script.phaser_phil import phaser_phil

MR_INPUTS = ("F_SIGF", "ENSEMBLES", "FIXENSEMBLES", "COMP_BY", "ASUFILE", "SEQUENCES",
             "ASU_PROTEIN_MW", "ASU_NUCLEICACID_MW", "SOLVENT_FRACTION", "SOLIN")


class phaser_pipeline_phil(PhilPluginScript):

    TASKNAME = "phaser_pipeline_phil"
    TASKCOMMAND = None
    MR_TASK = "phaser_mr_auto_phil"
    PHIL_PARAMS_FILE = phaser_mr_auto_phil.PHIL_PARAMS_FILE
    PHIL_MODE = phaser_mr_auto_phil.PHIL_MODE
    PHIL_MODE_PATH = phaser_mr_auto_phil.PHIL_MODE_PATH
    WHATNEXT = ["prosmart_refmac", "modelcraft", "coot_rebuild"]
    ASYNCHRONOUS = False
    LEGACY_PHIL_VALUES = phaser_phil.LEGACY_PHIL_VALUES

    ERROR_CODES = {
        200: {"description": "Phaser failed"},
        201: {"description": "Phaser found no solution"},
        204: {"description": "Sheetbend failed"},
        205: {"description": "Reindexing the data to the solution failed"},
        206: {"description": "Csymmatch failed"},
        207: {"description": "Phaser produced no coordinates"},
        210: {"description": "Refmac failed"},
        211: {"description": "Copying a sub-job output failed"},
    }

    def get_phil_exclude_scopes(self):
        # The same tree as the task's: what its shims write is left out here too
        return list(phaser_mr_auto_phil.PHIL_EXCLUDE_SCOPES) + phaser_mr_auto_phil.phil_shim_targets()

    def get_command_target(self):
        return None

    def validity(self):
        error = super().validity()
        # The task's own rules, on the same inputs
        probe = phaser_mr_auto_phil(workDirectory=str(self.getWorkDirectory()), name=self.MR_TASK)
        self._copy_inputs(probe)
        for report in probe.validity()._reports:
            if report["code"] in (111, 112, 113, 115):
                error.append(klass=self.TASKNAME, code=report["code"], details=report["details"],
                             name=str(report.get("name", "")).replace(f"{self.MR_TASK}.", f"{self.TASKNAME}."),
                             severity=report["severity"])
        return error

    # -- the run -----------------------------------------------------------
    def process(self):
        invalid = self.checkInputData()
        if invalid:
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        self.checkOutputData()
        self.xmlroot = etree.Element("PhaserPipeline")
        self.phaserPlugin = self.makePluginObject(self.MR_TASK)
        self._copy_inputs(self.phaserPlugin)
        self.hand_phil_to(self.phaserPlugin)
        self.phaserPlugin.xml_responders.append(self.phaserXMLUpdated)
        rv = self.phaserPlugin.process()
        self.phaserXMLUpdated(self.phaserPlugin.xmlroot)
        if rv != CPluginScript.SUCCEEDED:
            self.errorReport.extend(self.phaserPlugin.getErrorReport())
            self.appendErrorReport(200, severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        outputs = self.phaserPlugin.container.outputData
        if len(outputs.XYZOUT) == 0:
            self.appendErrorReport(201, "No solution to refine", severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.UNSATISFACTORY)
            return CPluginScript.UNSATISFACTORY
        self.harvestPhaser()
        f_sigf, freer = self.container.inputData.F_SIGF, self.container.inputData.FREERFLAG
        xyzin = self.container.outputData.XYZOUT[0]
        if bool(outputs.dataReindexed):
            if self.runPointless() != CPluginScript.SUCCEEDED:
                return CPluginScript.FAILED
            f_sigf, freer = self.container.outputData.F_SIGF_OUT, self.container.outputData.FREERFLAG_OUT
        if self.container.inputData.XYZIN_TARGET.isSet():
            if self.runCsymmatch() != CPluginScript.SUCCEEDED:
                return CPluginScript.FAILED
            xyzin = self.container.outputData.XYZOUT_CSYMMATCH
        if self.container.inputData.RUNSHEETBEND:
            if self.runSheetbend(f_sigf, freer, xyzin) != CPluginScript.SUCCEEDED:
                return CPluginScript.FAILED
            xyzin = self.container.outputData.XYZOUT_SHEETBEND
        if self.container.inputData.RUNREFMAC:
            if self.runRefmac(f_sigf, freer, xyzin) != CPluginScript.SUCCEEDED:
                return CPluginScript.FAILED
        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    def _copy_inputs(self, plugin):
        for name in MR_INPUTS:
            src = getattr(self.container.inputData, name, None)
            dst = getattr(plugin.container.inputData, name, None)
            if src is None or dst is None:
                continue
            if hasattr(src, "isSet") and not src.isSet():
                continue
            dst.set(src)

    # -- the record --------------------------------------------------------
    def phaserXMLUpdated(self, xmlroot):
        for old in self.xmlroot.findall("PhaserMrResults"):
            self.xmlroot.remove(old)
        self.xmlroot.append(deepcopy(xmlroot))
        self._writeXML()

    def appendXML(self, path, tag):
        for old in self.xmlroot.findall(tag):
            self.xmlroot.remove(old)
        try:
            self.xmlroot.append(CCP4Utils.openFileToEtree(path))
        except Exception:
            self.xmlroot.append(etree.Element(tag))
        self._writeXML()

    def _writeXML(self):
        target = self.makeFileName("PROGRAMXML")
        tmp = target + "_tmp"
        etree.ElementTree(self.xmlroot).write(tmp, pretty_print=True, xml_declaration=True,
                                              encoding="utf-8")
        os.replace(tmp, target)

    # -- the steps after Phaser ------------------------------------------------
    def harvestPhaser(self):
        outputs, mine = self.phaserPlugin.container.outputData, self.container.outputData
        self.harvestFile(outputs.SOLOUT, mine.SOLOUT)
        for name in ("XYZOUT", "MAPOUT", "DIFMAPOUT", "PHASEOUT"):
            for item in getattr(outputs, name):
                target = getattr(mine, name)
                target.append(target.makeItem())
                target[-1].fullPath = os.path.join(str(self.workDirectory), os.path.basename(str(item.fullPath)))
                self.harvestFile(item, target[-1])
        mine.dataReindexed.set(bool(outputs.dataReindexed))

    def harvestFile(self, src_item, dst_item):
        try:
            src, dst = str(src_item.fullPath), str(dst_item.fullPath)
            if Path(src).suffix.lower() != Path(dst).suffix.lower() and Path(src).suffix.lower() in (".pdb", ".cif", ".mmcif"):
                dst = str(Path(dst).with_suffix(Path(src).suffix.lower()))
                dst_item.setFullPath(dst)
            shutil.copyfile(src, dst)
            dst_item.annotation.set(src_item.annotation)
            dst_item.contentFlag.set(src_item.contentFlag)
            dst_item.subType.set(src_item.subType)
        except Exception as err:
            self.appendErrorReport(211, f"{src_item.fullPath} -> {dst_item.fullPath}: {err}",
                                   severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)

    def runPointless(self):
        try:
            plugin = self.makePluginObject("pointless_reindexToMatch")
            plugin.container.controlParameters.REFERENCE = "HKLIN_FMAP_REF"
            plugin.container.inputData.HKLIN_FMAP_REF.set(self.container.outputData.MAPOUT[0])
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
        except Exception as err:
            self.appendErrorReport(205, str(err), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def runCsymmatch(self):
        try:
            plugin = self.makePluginObject("csymmatch")
            plugin.container.inputData.XYZIN_QUERY.set(self.container.outputData.XYZOUT[0])
            plugin.container.inputData.XYZIN_TARGET.set(self.container.inputData.XYZIN_TARGET)
            rv = plugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(206, severity=CCP4ErrorHandling.SEVERITY_ERROR)
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
            self.harvestFile(plugin.container.outputData.XYZOUT, self.container.outputData.XYZOUT_CSYMMATCH)
            self.container.outputData.XYZOUT_CSYMMATCH.annotation.set("Coordinates moved to match the reference structure")
            self.appendXML(plugin.makeFileName("PROGRAMXML"), "Csymmatch")
        except Exception as err:
            self.appendErrorReport(206, str(err), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def runSheetbend(self, f_sigf, freer, xyzin):
        try:
            plugin = self.makePluginObject("sheetbend")
            plugin.container.inputData.XYZIN.set(xyzin)
            plugin.container.inputData.F_SIGF.set(f_sigf)
            if freer is not None and freer.isSet():
                plugin.container.inputData.FREERFLAG.set(freer)
            rv = plugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(204, severity=CCP4ErrorHandling.SEVERITY_ERROR)
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
            self.harvestFile(plugin.container.outputData.XYZOUT, self.container.outputData.XYZOUT_SHEETBEND)
            self.container.outputData.XYZOUT_SHEETBEND.annotation.set("Model after shift-field refinement")
            self.appendXML(plugin.makeFileName("PROGRAMXML"), "SheetbendResult")
        except Exception as err:
            self.appendErrorReport(204, str(err), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def runRefmac(self, f_sigf, freer, xyzin):
        try:
            plugin = self.makePluginObject("refmac")
            plugin.container.inputData.XYZIN.set(xyzin)
            plugin.container.inputData.F_SIGF.set(f_sigf)
            if freer is not None and freer.isSet():
                plugin.container.inputData.FREERFLAG.set(freer)
            cp = plugin.container.controlParameters
            cp.HYDROGENS.set("NO")
            cp.NCYCLES.set(10)
            cp.PHOUT.set(False)
            cp.USE_JELLY.set(True)
            cp.JELLY_SIGMA.set(0.05)
            cp.MAKE_NEW_LIGAND_EXIT.set(False)
            rv = plugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(210, severity=CCP4ErrorHandling.SEVERITY_ERROR)
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
            outputs, mine = plugin.container.outputData, self.container.outputData
            self.harvestFile(outputs.FPHIOUT, mine.MAPOUT_REFMAC)
            self.harvestFile(outputs.DIFFPHIOUT, mine.DIFMAPOUT_REFMAC)
            self.harvestFile(outputs.XYZOUT, mine.XYZOUT_REFMAC)
            self.appendXML(plugin.makeFileName("PROGRAMXML"), "REFMAC")
            mine.PERFORMANCEINDICATOR.set(outputs.PERFORMANCEINDICATOR)
        except Exception as err:
            self.appendErrorReport(210, str(err), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED
