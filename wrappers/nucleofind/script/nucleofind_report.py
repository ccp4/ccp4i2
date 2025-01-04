from report.CCP4ReportParser import *
import sys

class nucleofind_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'nucleofind'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        print("In nucleofind report")
        self.addFold()
        self.addDiv("Nucleofind finished.")

