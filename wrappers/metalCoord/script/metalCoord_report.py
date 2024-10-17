from report.CCP4ReportParser import *
from core import CCP4Modules
import sys,os

class metalCoord_report(Report):
    TASKNAME='metalCoord'

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus in ["Finished"]:
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)
            jobId = self.jobInfo.get("jobid", None)

            jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
        self.addDiv(style="clear:both;")
