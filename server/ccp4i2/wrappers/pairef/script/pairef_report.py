import os
import sys

from ccp4i2.report import Report


class pairef_report(Report):

    TASKNAME = 'pairef'
    USEPROGRAMXML = True
    SEPARATEDATA = True
    RUNNING = True
    WATCHED_FILE = os.path.join("pairef_project", "PAIREF_project.html")

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        
        if not jobStatus in ["Running", "Running remotely"]:
            self.addResults()
            self.append("Paired Refinement complete.")
        else:
            self.append("Paired Refinement is running.")
        pairef_html = os.path.normpath(os.path.join(self.jobInfo["fileroot"], "pairef_project", "PAIREF_project.html"))
        if os.path.exists(pairef_html):
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)
            pairef_url = (
                "/database/projectId/"
                + projectid
                + "/jobNumber/"
                + jobNumber
                + "/fileName/pairef_project/PAIREF_project.html"
            )
            pairefrFolder = self.addFold(label="pairef report", initiallyOpen=True)
            pairefrFolder.append(
                '<span style="font-size:110%">Click on the '
                "following link to display the pairef report </span>"
            )
            if not jobStatus in ["Running", "Running remotely"]:
                pairefrFolder.append('<a href="{0}">Open Results</a>'.format(pairef_url))
            else:
                if sys.platform == "win32":
                    import pathlib
                    pairef_html = pathlib.Path(pairef_html).as_uri()
                pairefrFolder.append('<a href="{0}">Open Results</a>'.format(pairef_html))
        else:
            self.append("The html report is not ready yet")
        return
