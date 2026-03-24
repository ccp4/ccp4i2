import os

from ccp4i2.core import CCP4Modules
from ccp4i2.report import Report


class prosmart_report(Report):
    TASKNAME='prosmart'

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus in ["Finished"]:
            jobId = self.jobInfo.get("jobid", None)

            jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
            prosmart_html = os.path.join(jobDirectory, "ProSMART_Results.html")
            if os.path.isfile(prosmart_html):
                self.addFileLink(
                    label='Open ProSMART Results',
                    relativePath='ProSMART_Results.html',
                    fileType='html',
                )

        self.addDiv(style="clear:both;")
