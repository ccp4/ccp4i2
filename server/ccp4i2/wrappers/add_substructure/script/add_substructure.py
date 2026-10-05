"""Add an anomalous substructure to a built model.

Model building drops the heavy atoms (ModelCraft rebuilds through Buccaneer,
which outputs protein only), and refining a model without strong scatterers
distorts it to explain them: for a mercury soak of HypF, R-free 0.295 with the
two sites and 0.437 without. This puts them back, once each: a site already a
model atom (S-SAD sulphurs, SeMet selenium) is not doubled, Se on a methionine
makes it MSE, sites from another job are first put in the model's origin and hand, and
residues model building traced into a heavy atom's density are removed (HypF:
ModelCraft's model 0.409 -> 0.267 with the three removed and the two Hg
added). The work is complete_model, shared with Crank2's XYZOUT_COMPLETE.
"""
import os

from ccp4i2.core import CCP4Utils
from ccp4i2.core.CCP4PluginScript import CPluginScript


class add_substructure(CPluginScript):
    TASKNAME = "add_substructure"
    ERROR_CODES = {
        201: {"description": "The substructure could not be added"},
        202: {"description": "The sites could not be put in the model's frame; none were added"},
    }

    def startProcess(self):
        from ccp4i2.pipelines.crank2.script.complete_model import complete_model
        inp, ctrl = self.container.inputData, self.container.controlParameters
        target = str(self.container.outputData.XYZOUT.fullPath)
        if not target.lower().endswith(".pdb"):
            target = os.path.splitext(target)[0] + ".pdb"
        try:
            # underscore: CData.__setattr__ would wrap a dict in a CData
            self._merge_report = complete_model(
                str(inp.XYZIN.fullPath), str(inp.XYZIN_SUB.fullPath), target,
                all_met_to_mse=bool(ctrl.ALL_MET_TO_MSE), find_origin=bool(ctrl.FIND_ORIGIN),
                remove_clashing=bool(ctrl.REMOVE_CLASHING))
        except Exception as err:  # noqa: BLE001 - reported, then the job fails
            self.appendErrorReport(201, f"{type(err).__name__}: {err}")
            return CPluginScript.FAILED
        self.container.outputData.XYZOUT.setFullPath(target)
        return CPluginScript.SUCCEEDED

    def postProcessCheck(self, processId=None):
        ok = os.path.isfile(str(self.container.outputData.XYZOUT.fullPath))
        return (CPluginScript.SUCCEEDED if ok else CPluginScript.FAILED), (0 if ok else 1), 0

    def processOutputFiles(self):
        from lxml import etree
        from ccp4i2.pipelines.crank2.script.complete_model import report_element
        report = self._merge_report
        root = etree.Element("AddSubstructure")
        root.append(report_element(report))
        with open(self.makeFileName("PROGRAMXML"), "w") as out:
            CCP4Utils.writeXML(out, etree.tostring(root, pretty_print=True))
        n_added, n_converted = len(report["added"]), len(report["converted"])
        parts = []
        if n_added:
            parts.append(f"{n_added} site{'s' if n_added != 1 else ''} added")
        if report["removed"]:
            parts.append(f"{len(report['removed'])} residues under them removed")
        if n_converted or report["all_mse"]:
            parts.append(f"{n_converted + report['all_mse']} MET made MSE")
        self.container.outputData.XYZOUT.annotation.set(
            (", ".join(parts) or "no sites added") + " - " + str(self.container.inputData.XYZIN.annotation))
        origin = report.get("origin") or {}
        if any("frame not established" in s for s in report["not_placed"]):
            self.appendErrorReport(202, origin.get("note", ""))
            return CPluginScript.UNSATISFACTORY
        return CPluginScript.SUCCEEDED
