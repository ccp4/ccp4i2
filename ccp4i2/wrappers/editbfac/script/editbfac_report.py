from report.CCP4ReportParser import *
import sys
from lxml import etree

class editbfac_report(Report):
    # Specify which gui task and/or pluginscript this applies to
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
