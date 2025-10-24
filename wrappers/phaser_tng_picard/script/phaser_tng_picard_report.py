from report.CCP4ReportParser import Report


class phaser_tng_picard_report(Report):
    TASKNAME = "phaser_tng_picard"
    RUNNING = False
    USEPROGRAMXML = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.addDiv(style="clear:both;")
