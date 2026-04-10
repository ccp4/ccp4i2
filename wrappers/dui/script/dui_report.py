from report.CCP4ReportParser import Report

class dui_report(Report):
    TASKNAME = 'dui'
    USEPROGRAMXML = False
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addResults()
