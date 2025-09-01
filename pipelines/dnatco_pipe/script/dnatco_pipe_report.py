"""
    dnatco_pipe_report.py: CCP4 GUI Project
    Copyright (C) 2025 MRC-LMB
    Author: Martin Maly
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
"""

from report.CCP4ReportParser import *
from core import CCP4Modules
from wrappers.dnatco.script.dnatco_report import dnatco_report
import os


class dnatco_pipe_report(Report):
    TASKNAME= 'dnatco_pipe'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return

        projectid = self.jobInfo.get("projectid", None)
        jobNumber = self.jobInfo.get("jobnumber", None)
        jobId = self.jobInfo.get("jobid", None)
        jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId)
        self.jobLog = os.path.join(jobDirectory, "log.txt")
        if jobStatus is not None and jobStatus.lower() == "running":
            self.runningReport(parent=self)
        else:
            self.defaultReport(parent=self)


    def runningReport(self, parent=None):
        if parent is None:
            parent = self
        if os.path.isfile(self.jobLog):
            jobLogFold = parent.addFold(label="DNATCO log", initiallyOpen=True)
            jobLogFold.addPre("DNATCO is running...")


    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        self.addDiv(style="clear:both;")  # gives space for the title

        jobId = self.jobInfo.get("jobid", None)
        jobDirectory = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId = jobId)

        jobDirectory1 = os.path.join(jobDirectory, "job_1")  # TODO: derive from jobId?
        jobDirectory2 = os.path.join(jobDirectory, "job_2")  # TODO: derive from jobId?
        jobDirectories = [jobDirectory1]
        if os.path.isdir(jobDirectory2):
            jobDirectories.append(jobDirectory2)

        dnatco_report1 = dnatco_report()
        dnatco_report1.defaultReport(parent, jobDirectories)