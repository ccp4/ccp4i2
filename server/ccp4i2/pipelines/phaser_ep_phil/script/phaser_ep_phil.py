"""SAD phasing with Phaser, then density modification and optionally
model building -- a pipeline that hosts Phaser's PHIL.

Its parameters are the EP_AUTO task's (same master, same mode, same
exclusions), handed to the sub-job whole. Its typed inputs are the task's
plus a Free-R set, the choice of a SHELX substructure search, which hands
to phase, and the steps afterwards: parrot on each hand phased, and
ModelCraft on parrot's phases if asked for.
"""
import os
import shutil
from copy import deepcopy

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.PhilPluginScript import PhilPluginScript
from ccp4i2.wrappers.phaser_ep_auto_phil.script.phaser_ep_auto_phil import phaser_ep_auto_phil

EP_INPUTS = ("F_SIGF", "WAVELENGTH", "XYZIN_HA", "PARTIAL_BY", "XYZIN_PARTIAL", "ELEMENTS",
             "LLGC_CYCLES", "PURE_ANOMALOUS", "COMP_BY", "ASUFILE", "SEQUENCES", "SOLVENT_FRACTION")


class phaser_ep_phil(PhilPluginScript):

    TASKNAME = "phaser_ep_phil"
    TASKCOMMAND = None
    EP_TASK = "phaser_ep_auto_phil"
    PHIL_PARAMS_FILE = phaser_ep_auto_phil.PHIL_PARAMS_FILE
    PHIL_MODE = phaser_ep_auto_phil.PHIL_MODE
    LEGACY_PHIL_VALUES = phaser_ep_auto_phil.LEGACY_PHIL_VALUES
    LEGACY_INPUT_RENAMES = phaser_ep_auto_phil.LEGACY_INPUT_RENAMES
    PHIL_MODE_PATH = phaser_ep_auto_phil.PHIL_MODE_PATH
    WHATNEXT = ["parrot", "modelcraft", "coot_rebuild"]
    ASYNCHRONOUS = False

    ERROR_CODES = {
        201: {"description": "SHELXC/D failed"},
        203: {"description": "Phaser failed"},
        206: {"description": "Parrot failed"},
        207: {"description": "ModelCraft failed"},
        211: {"description": "Copying a sub-job output failed"},
    }

    def get_phil_exclude_scopes(self):
        # The task's tree exactly: its exclusions plus what its shims write, and
        # the hand, which the pipeline decides
        targets = ["phaser.crystal", "phaser.keywords.partial", "phaser.keywords.llgcompletion",
                   "phaser.composition.chain", "phaser.composition.solvent",
                   "phaser.keywords.general.root"]
        return list(phaser_ep_auto_phil.PHIL_EXCLUDE_SCOPES) + targets + ["phaser.keywords.hand"]

    def get_command_target(self):
        return None

    def validity(self):
        error = super().validity()
        inp = self.container.inputData
        name = f"{self.TASKNAME}.container.inputData"
        if str(inp.SUBSTRUCTURE_BY) == "SITES" and not inp.XYZIN_HA.isSet() and not (
                str(inp.PARTIAL_BY) == "MODEL" and inp.XYZIN_PARTIAL.isSet()):
            error.append(klass=self.TASKNAME, code=113,
                         details="Give heavy-atom sites, a partial model, or ask for a SHELX search.",
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

    # -- the run -----------------------------------------------------------
    def process(self):
        invalid = self.checkInputData()
        if str(self.container.inputData.SUBSTRUCTURE_BY) == "SEARCH" and "XYZIN_HA" in invalid:
            invalid.remove("XYZIN_HA")
        if invalid:
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        self.checkOutputData()
        self.xmlroot = etree.Element("PhaserEP")
        if str(self.container.inputData.SUBSTRUCTURE_BY) == "SEARCH":
            if self.runShelx() != CPluginScript.SUCCEEDED:
                return CPluginScript.FAILED
        if self.runPhaser() != CPluginScript.SUCCEEDED:
            return CPluginScript.FAILED
        outputs, mine = self.phaserPlugin.container.outputData, self.container.outputData
        for i in range(len(outputs.XYZOUT)):
            self.harvestInto(outputs.XYZOUT[i], mine.XYZOUT)
            self.harvestInto(outputs.ABCDOUT[i], mine.ABCDOUT)
            self.harvestInto(outputs.MAPOUT[i], mine.MAPOUT)
        hands = self.hands()
        self._parrot = {}   # underscore: a bare name is coerced to CData on a plugin
        if self.container.controlParameters.RUNPARROT:
            for hand in hands:
                if self.runParrot(hand) != CPluginScript.SUCCEEDED:
                    return CPluginScript.FAILED
        if self.container.controlParameters.RUNMODELCRAFT and self.container.controlParameters.RUNPARROT:
            for hand in hands:
                if self.runModelCraft(hand) != CPluginScript.SUCCEEDED:
                    return CPluginScript.FAILED
        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    def hands(self):
        """Which of the phased hands to carry on with: the task wrote one
        (index 0) or, for 'both', two (0 original, 1 inverted)."""
        n = len(self.phaserPlugin.container.outputData.ABCDOUT)
        return list(range(n))

    def runShelx(self):
        try:
            plugin = self.makePluginObject("ShelxCD")
            plugin.container.inputData.SAD.set(self.container.inputData.F_SIGF)
            plugin.container.controlParameters.MODE.set("SAD")
            plugin.container.controlParameters.SFAC.set(str(self.container.inputData.SFAC))
            plugin.container.controlParameters.NTRY.set(int(self.container.inputData.NTRY))
            plugin.container.controlParameters.FIND.set(int(self.container.inputData.FIND))
            rv = plugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(201, severity=CCP4ErrorHandling.SEVERITY_ERROR)
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
            self.container.inputData.XYZIN_HA.set(plugin.container.outputData.XYZOUT)
            self.appendXML(plugin.makeFileName("PROGRAMXML"), "ShelxCD")
        except Exception as err:
            self.appendErrorReport(201, str(err), severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def runPhaser(self):
        self.phaserPlugin = self.makePluginObject(self.EP_TASK)
        for name in EP_INPUTS:
            src = getattr(self.container.inputData, name, None)
            dst = getattr(self.phaserPlugin.container.inputData, name, None)
            if src is None or dst is None or (hasattr(src, "isSet") and not src.isSet()):
                continue
            dst.set(src)
        self.hand_phil_to(self.phaserPlugin)
        self.phaserPlugin.set_phil("phaser.keywords.hand", str(self.container.inputData.HAND))
        self.phaserPlugin.xml_responders.append(self.phaserXMLUpdated)
        rv = self.phaserPlugin.process()
        self.phaserXMLUpdated(self.phaserPlugin.xmlroot)
        if rv != CPluginScript.SUCCEEDED:
            self.errorReport.extend(self.phaserPlugin.getErrorReport())
            self.appendErrorReport(203, severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def runParrot(self, hand):
        label = "original" if hand == 0 else "inverted"
        try:
            plugin = self.makePluginObject("parrot")
            plugin.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF)
            if self.container.inputData.ASUFILE.isSet():
                plugin.container.inputData.ASUIN.set(self.container.inputData.ASUFILE)
            plugin.container.inputData.ABCD.set(self.phaserPlugin.container.outputData.ABCDOUT[hand])
            rv = plugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(206, f"{label} hand", severity=CCP4ErrorHandling.SEVERITY_ERROR)
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
            self._parrot[hand] = plugin
            self.appendXML(plugin.makeFileName("PROGRAMXML"), "ParrotResult", hand=label)
            self.harvestInto(plugin.container.outputData.ABCDOUT, self.container.outputData.ABCDOUT,
                             annotation=f"Phases from density modification - {label} hand")
            self.harvestInto(plugin.container.outputData.FPHIOUT, self.container.outputData.FPHIOUT,
                             annotation=f"Map coefficients from density modification - {label} hand")
        except Exception as err:
            self.appendErrorReport(206, f"{label} hand: {err}", severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def runModelCraft(self, hand):
        label = "original" if hand == 0 else "inverted"
        try:
            plugin = self.makePluginObject("modelcraft")
            plugin.container.inputData.F_SIGF.set(self.container.inputData.F_SIGF)
            if self.container.inputData.FREERFLAG.isSet():
                plugin.container.inputData.FREERFLAG.set(self.container.inputData.FREERFLAG)
            if self.container.inputData.ASUFILE.isSet():
                plugin.container.inputData.ASUIN.set(self.container.inputData.ASUFILE)
            plugin.container.inputData.PHASES.set(self._parrot[hand].container.outputData.ABCDOUT)
            cp = plugin.container.controlParameters
            cp.BASIC.set(True)
            cp.USE_MODEL_PHASES.set(False)
            cp.CYCLES.set(int(self.container.controlParameters.MODELCRAFT_ITERATIONS))
            cp.STOP_CYCLES.set(2)
            rv = plugin.process()
            if rv != CPluginScript.SUCCEEDED:
                self.appendErrorReport(207, f"{label} hand", severity=CCP4ErrorHandling.SEVERITY_ERROR)
                self.reportStatus(CPluginScript.FAILED)
                return CPluginScript.FAILED
            node = etree.Element("ModelCraft")
            node.text = str(plugin.getWorkDirectory())
            self.updateXml(node, "ModelCraft", hand=label)
            self.harvestInto(plugin.container.outputData.XYZOUT, self.container.outputData.XYZOUT,
                             annotation=f"Autobuilt model - {label} hand")
        except Exception as err:
            self.appendErrorReport(207, f"{label} hand: {err}", severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    # -- files and the record --------------------------------------------------
    def harvestInto(self, src_item, target_list, annotation=None):
        try:
            target_list.append(target_list.makeItem())
            dst_item = target_list[-1]
            src = str(src_item.fullPath)
            # Named by list and position: two hands' parrot outputs share a basename
            suffix = os.path.splitext(src)[1]
            dst_item.fullPath = os.path.join(
                str(self.workDirectory), f"{target_list.objectName()}_{len(target_list)}{suffix}")
            shutil.copyfile(src, str(dst_item.fullPath))
            dst_item.annotation.set(annotation if annotation else src_item.annotation)
            dst_item.contentFlag.set(src_item.contentFlag)
            dst_item.subType.set(src_item.subType)
        except Exception as err:
            self.appendErrorReport(211, f"{src_item.fullPath}: {err}", severity=CCP4ErrorHandling.SEVERITY_ERROR)
            self.reportStatus(CPluginScript.FAILED)

    def phaserXMLUpdated(self, xmlroot):
        self.updateXml(deepcopy(xmlroot), "PhaserEpResults")

    def appendXML(self, path, tag, hand=None):
        try:
            node = CCP4Utils.openFileToEtree(path)
        except Exception:
            node = etree.Element(tag)
        self.updateXml(node, tag, hand)

    def updateXml(self, node, tag, hand=None):
        if hand is None:
            for old in self.xmlroot.findall(tag):
                self.xmlroot.remove(old)
            self.xmlroot.append(node)
        else:
            hand_node = self.xmlroot.find(hand)
            if hand_node is None:
                hand_node = etree.SubElement(self.xmlroot, hand)
            for old in hand_node.findall(tag):
                hand_node.remove(old)
            hand_node.append(node)
        target = self.makeFileName("PROGRAMXML")
        tmp = target + "_tmp"
        etree.ElementTree(self.xmlroot).write(tmp, pretty_print=True, xml_declaration=True,
                                              encoding="utf-8")
        os.replace(tmp, target)
