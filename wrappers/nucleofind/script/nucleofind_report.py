from report.CCP4ReportParser import *
import sys

class nucleofind_report(Report):
    TASKNAME = 'nucleofind'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        self.addDiv(style="clear:both;")
        # self.addReference()
        
        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")


