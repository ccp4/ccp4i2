from report.CCP4ReportParser import Report

class coot1_report(Report):
    TASKNAME = 'coot1'
    RUNNING = False
    USEPROGRAMXML = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super(). __init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.addDiv(style="clear:both;")
        self.addText(text='Happily finished')
