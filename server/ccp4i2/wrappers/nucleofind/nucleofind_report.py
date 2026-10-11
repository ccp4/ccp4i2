from ccp4i2.report.CCP4ReportParser import Report


class nucleofind_report(Report):
    TASKNAME = 'nucleofind'
    RUNNING = True
    USEPROGRAMXML = False
    WATCHED_FILE= 'log_err.txt'

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        self.addDiv()

        if jobStatus in ["Running", "Running remotely"]:
            self.addText(text="The job is currently running.", tone="emphasis")
            self.addDiv()
