import os

from ccp4i2.report import Report


class editbfac_report(Report):
    TASKNAME = 'editbfac'
    RUNNING = False

    def __init__(self, xmlnode=None,jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.defaultReport()
    
    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.addResults()
        parent.append("Edit B-factors finished.")
        logf = os.path.join(self.getJobFolder(), "log.txt")
        bfacHtmlFold = parent.addFold(label='Report from cctbx/i2', initiallyOpen=True)
        for lline in open(logf, 'r'):
            bfacHtmlFold.append('<span style="font-size:95%">' + lline + '</span>')
            bfacHtmlFold.append('<br></br>')
