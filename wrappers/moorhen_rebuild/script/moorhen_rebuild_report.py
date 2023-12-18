from report.CCP4ReportParser import Report
import sys

class moorhen_rebuild_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'moorhen_rebuild'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.addText(text='Happily finished')
