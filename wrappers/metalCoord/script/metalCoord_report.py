from report.CCP4ReportParser import *
from core import CCP4Modules
import sys, os


class metalCoord_report(Report):
    TASKNAME='metalCoord'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return

        # projectid = self.jobInfo.get("projectid", None)
        # jobNumber = self.jobInfo.get("jobnumber", None)
        jobId = self.jobInfo.get("jobid", None)
        jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)
        self.jobLog = os.path.join(jobDirectory, "log.txt")
        if jobStatus is not None and jobStatus.lower() == "running":
            self.runningReport(parent=self)
        else:  # elif jobStatus in ["Finished"]:
            self.defaultReport(parent=self)


    def runningReport(self, parent=None):
        if parent is None:
            parent = self
        if os.path.isfile(self.jobLog):
            jobLogFold = parent.addFold(label="MetalCoord log", initiallyOpen=True)
            jobLogFold.addPre("MetalCoord is running...")


    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        self.addDiv(style="clear:both;")  # gives space for the title

        for i, site in enumerate(self.xmlnode.findall(".//site")):
            for j, symmClass in enumerate(site.findall(".//ligands")):
                table = parent.addTable(xmlnode=symmClass)  # transpose=True)
                table.addData(title="Symmetry class", select='class')
                table.addData(title="Procrustes", select='procrustes')
                table.addData(title="Coordination", select='coordination')
                table.addData(title="Count", select='count')
                table.addData(title="Description", select='description')
                table = parent.addTable(xmlnode=symmClass)  # transpose=True)
                for k, entry in enumerate(symmClass.findall(".//base")):
                    table = parent.addTable(xmlnode=entry)  # transpose=True)
                    table.addData(title="Distance", select='distance')
                    table.addData(title="St. dev.", select='std')
                for k, entry in enumerate(symmClass.findall(".//pdb")):
                    table = parent.addTable(xmlnode=entry)  # transpose=True)
                    table.addData(title="Distance", select='distance')
                    table.addData(title="St. dev.", select='std')
                for k, entry in enumerate(symmClass.findall(".//angles")):
                    table = parent.addTable(xmlnode=entry)  # transpose=True)
                    table.addData(title="Angle", select='angle')
                    table.addData(title="St. dev.", select='std')

        self.addDiv(style="clear:both;")
