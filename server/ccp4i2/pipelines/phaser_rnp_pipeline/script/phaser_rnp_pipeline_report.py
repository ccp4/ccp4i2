from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script.phaser_MR_AUTO_report import (
    phaser_MR_AUTO_report,
)
from ccp4i2.report import Report
from ccp4i2.wrappers.pointless.script.pointless_report import pointless_report
from ccp4i2.wrappers.refmac.script.refmac_report import refmac_report


class phaser_rnp_pipeline_report(Report):
    TASKNAME = 'phaser_rnp_pipeline'
    RUNNING = True
    
    def __init__(self,jobStatus=None,**kw):
        Report.__init__(self, jobStatus=jobStatus, **kw)

        if jobStatus == None or jobStatus.lower() =='nooutput': return
        self.addDiv(style='clear:both;')

        try:
            pointlessNode = self.xmlnode.findall('POINTLESS')[0]
        except:
            return
        pointlessFold = self.addFold(label='Reindexing results from pointless',brief='Reindex')
        pointlessReport = pointless_report(xmlnode=pointlessNode, jobStatus='nooutput')
        pointlessReport.BestReindex(parent=pointlessFold)
        
        try:
            pmaNode = self.xmlnode.findall('PhaserMrResults')[0]
        except:
            return
        
        phaser_MRAReport = phaser_MR_AUTO_report(xmlnode=pmaNode, jobStatus='nooutput')
        
        phaser_MRAReport.drawContent(jobStatus=jobStatus, parent=self)

        if len(self.xmlnode.findall('REFMAC')) == 0: return
        refmacNode = self.xmlnode.findall('REFMAC')[0]

        rM_report = refmac_report(xmlnode=refmacNode, jobStatus='nooutput')

        rM_report.addSummary(parent=self)
