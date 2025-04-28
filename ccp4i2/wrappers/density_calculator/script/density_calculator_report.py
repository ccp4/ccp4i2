from ....report.CCP4ReportParser import Report


class density_calculator_report(Report):
    TASKNAME = "density_calculator"
    RUNNING = False
    USEPROGRAMXML = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
