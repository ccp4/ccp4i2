"""The EP pipeline's report: SHELX if it ran, the EP task's record, then
parrot and ModelCraft per hand."""
import json
import os

from ccp4i2.report import Report
from ccp4i2.wrappers.phaser_ep_auto_phil.script.phaser_ep_auto_phil_report import phaser_ep_auto_phil_report


class phaser_ep_phil_report(Report):
    TASKNAME = "phaser_ep_phil"
    RUNNING = True
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is None or jobStatus.lower() == "nooutput":
            return
        self.drawContent(jobStatus, self)

    def drawContent(self, jobStatus=None, parent=None):
        if parent is None:
            parent = self
        shelx = self.xmlnode.find("ShelxCD")
        if shelx is not None:
            try:
                from ccp4i2.wrappers.ShelxCD.script.ShelxCD_report import ShelxCD_report
                fold = parent.addFold(label="Substructure search (SHELXC/D)", initiallyOpen=False)
                ShelxCD_report(xmlnode=shelx, jobStatus="nooutput").drawContent(jobStatus=jobStatus, parent=fold)
            except Exception:
                parent.addText(text="SHELXC/D ran (report unavailable)")
        node = self.xmlnode.find("PhaserEpResults")
        if node is None:
            parent.addText(text="Phaser has not reported yet" if (jobStatus or "").lower() == "running"
                           else "No Phaser record", style="color:orange;")
            return
        phaser_ep_auto_phil_report(xmlnode=node, jobStatus="nooutput").drawContent(jobStatus=jobStatus, parent=parent)
        for hand, label in (("original", "Original hand"), ("inverted", "Inverted hand")):
            parrot = self.xmlnode.find(f"{hand}/ParrotResult")
            if parrot is not None:
                from ccp4i2.wrappers.parrot.script.parrot_report import parrot_report
                fold = parent.addFold(label=f"Density modification: {label}", initiallyOpen=False)
                parrot_report(xmlnode=parrot, jobStatus="nooutput").defaultReport(parent=fold)
                fold.addDiv(style="clear:both;")
            mc = self.xmlnode.find(f"{hand}/ModelCraft")
            if mc is not None and mc.text:
                from ccp4i2.wrappers.modelcraft.script.modelcraft_report import modelcraft_report
                fold = parent.addFold(label=f"Model building: {label}", initiallyOpen=False)
                report = modelcraft_report(jobInfo={"fileroot": ""})
                path = os.path.join(mc.text, "modelcraft", "modelcraft.json")
                if os.path.exists(path):
                    with open(path, encoding="utf-8") as stream:
                        report.json = json.load(stream)
                    report.add_running_job(parent=fold)
                    report.add_table(parent=fold)
                    report.add_message(jobStatus, parent=fold)
                    report.add_graph(parent=fold)
                fold.addDiv(style="clear:both;")
