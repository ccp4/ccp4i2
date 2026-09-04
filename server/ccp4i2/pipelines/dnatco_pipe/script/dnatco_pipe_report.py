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

import os

from ccp4i2.report import Report
from ccp4i2.wrappers.dnatco.script.dnatco_report import draw_dnatco_report


class dnatco_pipe_report(Report):
    TASKNAME = 'dnatco_pipe'
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        if jobStatus.lower() == 'running':
            self.runningReport(parent=self)
        else:
            self.defaultReport(parent=self)

    def runningReport(self, parent=None):
        if parent is None:
            parent = self
        fold = parent.addFold(label="DNATCO log", initiallyOpen=True)
        fold.addPre("DNATCO is running...")

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.addDiv(style="clear:both;")  # gives space for the title
        filenames = self.jobInfo.get('filenames', {}) if self.jobInfo else {}
        cif_paths = [filenames.get('CIFOUT1')]
        json_paths = [filenames.get('JSONOUT1')]
        # The second model's files exist only when the job compared two models.
        if filenames.get('CIFOUT2') and os.path.isfile(str(filenames['CIFOUT2'])):
            cif_paths.append(filenames['CIFOUT2'])
            json_paths.append(filenames.get('JSONOUT2'))
        draw_dnatco_report(parent, cif_paths, json_paths)
