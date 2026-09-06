"""The pipeline's report: the reindexing, the RNP task's record, then
whatever ran after."""
from ccp4i2.pipelines.phaser_pipeline_phil.script.phaser_pipeline_phil_report import phaser_pipeline_phil_report
from ccp4i2.wrappers.pointless.script.pointless_report import pointless_report


class phaser_rnp_pipeline_phil_report(phaser_pipeline_phil_report):
    TASKNAME = "phaser_rnp_pipeline_phil"

    def drawContent(self, jobStatus=None, parent=None):
        if parent is None:
            parent = self
        node = self.xmlnode.find("POINTLESS")
        if node is not None:
            fold = parent.addFold(label="Reindexing to match the parent model", initiallyOpen=False,
                                  brief="Reindex")
            pointless_report(xmlnode=node, jobStatus="nooutput").BestReindex(parent=fold)
        super().drawContent(jobStatus=jobStatus, parent=parent)
