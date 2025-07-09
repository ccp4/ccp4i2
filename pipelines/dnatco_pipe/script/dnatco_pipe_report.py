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
from pathlib import Path


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
        self.jobLog = str(Path(jobDirectory) / "log.txt")
        if jobStatus is not None and jobStatus.lower() == "running":
            self.runningReport(parent=self)
        else:
            self.defaultReport(parent=self)


    def runningReport(self, parent=None):
        if parent is None:
            parent = self
        if Path(self.jobLog).is_file():
            jobLogFold = parent.addFold(label="DNATCO log", initiallyOpen=True)
            jobLogFold.addPre("DNATCO is running...")


    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        self.addDiv(style="clear:both;")  # gives space for the title

        ciffile1Path = str(self.jobInfo['filenames']['CIFOUT1'])
        ciffile2Path = str(self.jobInfo['filenames']['CIFOUT2'])
        ciffilePaths = [ciffile1Path]
        if Path(ciffile2Path).is_file():
            ciffilePaths.append(ciffile2Path)
        dnatco_report1 = dnatco_report()
        dnatco_report1.defaultReport(parent, ciffilePaths)