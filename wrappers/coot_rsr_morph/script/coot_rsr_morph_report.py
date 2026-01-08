from report.CCP4ReportParser import Report

class coot_rsr_morph_report(Report):
    TASKNAME = 'coot_rsr_morph'
    RUNNING = False
    USEPROGRAMXML = False

    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        super().__init__(xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        self.addText(text='Coot real space morphing finished. Full reporting is not yet available in this task.')
