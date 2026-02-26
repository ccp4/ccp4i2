import os
import sys
from pathlib import Path

from ccp4i2.report import Report


class pairef_report(Report):
    TASKNAME = 'pairef'
    USEPROGRAMXML = False
    SEPARATEDATA = True
    RUNNING = True
    WATCHED_FILE = os.path.join("pairef_project", "PAIREF_project.html")

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super().__init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        if jobStatus not in ["Running", "Running remotely"]:
            self.addResults()
            self.append("Paired Refinement complete.")
        else:
            self.append("Paired Refinement is running.")
        pairef_html = Path(self.jobInfo["fileroot"], "pairef_project", "PAIREF_project.html").resolve()
        if pairef_html.exists():
            projectid = self.jobInfo["projectid"]
            jobNumber = self.jobInfo["jobnumber"]
            url = f"/database/projectId/{projectid}/jobNumber/{jobNumber}/fileName/pairef_project/PAIREF_project.html"
            fold = self.addFold(label="pairef report", initiallyOpen=True)
            fold.append('<span style="font-size:110%">Click on the following link to display the pairef report </span>')
            if jobStatus not in ["Running", "Running remotely"]:
                fold.append(f'<a href="{url}">Open Results</a>')
            else:
                if sys.platform == "win32":
                    pairef_html = pairef_html.as_uri()
                fold.append(f'<a href="{pairef_html}">Open Results</a>')
        else:
            self.append("The html report is not ready yet")
