from ....report.CCP4ReportParser import Report


class add_fractional_coords_report(Report):
    TASKNAME = "add_fractional_coords"
    RUNNING = True
    USEPROGRAMXML = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")
