from ccp4i2.report import Report
from ccp4i2.wrappers.csymmatch.script.csymmatch_report import csymmatch_report
from ccp4i2.wrappers.refmac.script.refmac_report import refmac_report
from ccp4i2.wrappers.sheetbend.script.sheetbend_report import sheetbend_report


class phaser_pipeline_report(Report):
    TASKNAME = 'phaser_pipeline'
    RUNNING = True
    SEPARATEDATA=True
    
    def __init__(self,jobStatus=None,**kw):
        Report.__init__(self, jobStatus=jobStatus, **kw)
        
        try:
            # filenames is 'the data keyed by task parameter name- not necessarilly all files but should all have __str__ methods'
            self.MODE_TY = self.jobInfo['filenames']['MODE_TY']
        except KeyError:
            # old versions of task
            self.MODE_TY = None 
        if self.MODE_TY == "MR_FRF":
            from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_FRF.script.phaser_MR_FRF_report import (
                phaser_MR_FRF_report,
            )
        elif self.MODE_TY == "MR_FTF":
            from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_FTF.script.phaser_MR_FTF_report import (
                phaser_MR_FTF_report,
            )
        elif self.MODE_TY == "MR_PAK":
            from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_PAK.script.phaser_MR_PAK_report import (
                phaser_MR_PAK_report,
            )
        else:
            from ccp4i2.pipelines.phaser_pipeline.wrappers.phaser_MR_AUTO.script.phaser_MR_AUTO_report import (
                phaser_MR_AUTO_report,
            )

        if jobStatus == None or jobStatus.lower() =='nooutput': return
        self.addDiv(style='clear:both;')

        if len(self.xmlnode.findall('PhaserMrResults')):
            pmaNode = self.xmlnode.findall('PhaserMrResults')[0]
        else:
            return
        
        if self.MODE_TY == "MR_FRF":
            phaser_MRFRFReport = phaser_MR_FRF_report(xmlnode=pmaNode, jobStatus='nooutput')
            phaser_MRFRFReport.drawContent(jobStatus=jobStatus, parent=self)
        elif self.MODE_TY == "MR_FTF":
            phaser_MRFTFReport = phaser_MR_FTF_report(xmlnode=pmaNode, jobStatus='nooutput')
            phaser_MRFTFReport.drawContent(jobStatus=jobStatus, parent=self)
        elif self.MODE_TY == "MR_PAK":
            phaser_MRPAKReport = phaser_MR_PAK_report(xmlnode=pmaNode, jobStatus='nooutput')
            phaser_MRPAKReport.drawContent(jobStatus=jobStatus, parent=self)
        else:
            phaser_MRAReport = phaser_MR_AUTO_report(xmlnode=pmaNode, jobStatus='nooutput')
            phaser_MRAReport.drawContent(jobStatus=jobStatus, parent=self)

        if len(self.xmlnode.findall('Csymmatch'))>0:
            csymmatchNode = self.xmlnode.findall('Csymmatch')[0]
            csymmatchReport = csymmatch_report(xmlnode=csymmatchNode, jobStatus='nooutput')
            csymmatchFold = self.addFold(label='Output from CSYMMATCH',initiallyOpen=False,brief='Symmatch')
            csymmatchReport.drawContent(jobStatus=jobStatus, parent=csymmatchFold)
        
        if len(self.xmlnode.findall('SheetbendResult'))>0:
            sheetbendNode = self.xmlnode.findall('SheetbendResult')[0]
            sheetbendReport = sheetbend_report(xmlnode=sheetbendNode, jobStatus='nooutput')
            sheetbendFold = self.addFold(label='Shift field refinement',initiallyOpen=True,brief='Sheetbend')
            sheetbendReport.defaultReport(parent=sheetbendFold)
        
        if len(self.xmlnode.findall('REFMAC'))>0:
            refmacNode = self.xmlnode.findall('REFMAC')[0]
        else:
            return

        rM_report = refmac_report(xmlnode=refmacNode, jobStatus='nooutput')

        rM_report.addSummary(parent=self)


