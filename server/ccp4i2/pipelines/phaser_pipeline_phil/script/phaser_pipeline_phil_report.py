"""The pipeline's report: the MR task's record, then whatever ran after."""
from ccp4i2.report import Report
from ccp4i2.wrappers.csymmatch.script.csymmatch_report import csymmatch_report
from ccp4i2.wrappers.refmac.script.refmac_report import refmac_report
from ccp4i2.wrappers.sheetbend.script.sheetbend_report import sheetbend_report
from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil_report import phaser_mr_auto_phil_report


class phaser_pipeline_phil_report(Report):
    TASKNAME = "phaser_pipeline_phil"
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
        node = self.xmlnode.find("PhaserMrResults")
        if node is None:
            parent.addText(text="Phaser has not reported yet" if (jobStatus or "").lower() == "running"
                           else "No Phaser record", style="color:orange;")
            return
        phaser_mr_auto_phil_report(xmlnode=node, jobStatus="nooutput").drawContent(jobStatus=jobStatus, parent=parent)
        if self.xmlnode.find("Csymmatch") is not None:
            fold = parent.addFold(label="Output from csymmatch", initiallyOpen=False, brief="Symmatch")
            csymmatch_report(xmlnode=self.xmlnode.find("Csymmatch"), jobStatus="nooutput").drawContent(
                jobStatus=jobStatus, parent=fold)
        if self.xmlnode.find("SheetbendResult") is not None:
            fold = parent.addFold(label="Shift-field refinement", initiallyOpen=True, brief="Sheetbend")
            sheetbend_report(xmlnode=self.xmlnode.find("SheetbendResult"), jobStatus="nooutput").defaultReport(parent=fold)
        if self.xmlnode.find("REFMAC") is not None:
            refmac_report(xmlnode=self.xmlnode.find("REFMAC"), jobStatus="nooutput").addSummary(parent=parent)
