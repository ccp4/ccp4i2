from report.CCP4ReportParser import *
import sys
from pathlib import Path
class nucleofind_report(Report):
    TASKNAME = 'nucleofind'
    RUNNING = True
    USEPROGRAMXML = False
    WATCHED_FILE= 'log_err.txt'
    UPDATE_INTERVAL = 1
    
    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        
        self.addDiv(style="clear:both;")

        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")
            self.addDiv(style="clear:both;")

        # log_path = Path(jobInfo["fileroot"], "log_err.txt")
        # print(f"Job info fileroot: {jobInfo['fileroot']}")
        # summaryFold = self.addFold(
        #     label="Prediction", initiallyOpen=True, brief="Prediction"
        # )
        
        # print(f"Log path: {log_path}")

        # if log_path.exists():
        #     print(f"Log path exists: {log_path}")
        #     with open(log_path, encoding="utf-8") as text:
        #         for line in text:
        #             text = line.strip()
        #             summaryFold.append(f"<p>{text}</p>")


