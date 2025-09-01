"""
    dnatco_pipe.py: CCP4 GUI Project
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
import shutil
import xml.etree.ElementTree as ET
from PySide2 import QtCore
from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import *
from core import CCP4Utils


class dnatco_pipe(CPluginScript):

    TASKMODULE = 'pipelines'  # Where this plugin will appear on gui
    TASKNAME = 'dnatco_pipe'  # Task name - should be same as class name
    TASKVERSION = 0.1         # Version of this plugin
    MAINTAINER = 'martin.maly@mrc-lmb.cam.ac.uk'

    ERROR_CODES = { 201 : { 'description' : 'No output restraint file from DNATCO' },
                    202 : { 'description' : 'No output extended mmCIF file from DNATCO' },
                    203 : { 'description' : 'No output PDF report file from DNATCO' },
                    204 : { 'description' : 'Failed to dnatcoify structure' },
                    }


    def process(self):
        self.dnatco1 = self.dnatcoCreate(self.container.inputData.XYZIN1)
        self.dnatco1.process()
        if (
            bool(self.container.controlParameters.TOGGLE_XYZIN2)
            and os.path.isfile(str(self.container.inputData.XYZIN2.fullPath))
        ):
            self.dnatco2 = self.dnatcoCreate(self.container.inputData.XYZIN2)
            self.dnatco2.process()
        self.process_finish()


    def dnatcoCreate(self, model):
        dnatco = self.makePluginObject('dnatco')
        dnatco.container.inputData.XYZIN.set(model)
        dnatco.container.controlParameters.copyData(self.container.controlParameters)
        # dnatco.container.controlParameters.GENERATE_RESTRAINTS.set(self.container.controlParameters.GENERATE_RESTRAINTS)
        # dnatco.container.controlParameters.MAX_RMSD.set(self.container.controlParameters.MAX_RMSD)
        # dnatco.container.controlParameters.RESTRAINTS_SIGMA.set(self.container.controlParameters.RESTRAINTS_SIGMA)
        self.connectSignal(dnatco, 'finished', self.dnatcoFinished)
        dnatco.waitForFinished = -1
        return dnatco


    @QtCore.Slot(dict)
    def dnatcoFinished(self, statusDict):
        status = statusDict['finishStatus']
        if status == CPluginScript.FAILED:
            self.reportStatus(status)


    def process_finish(self):
        xmlText = ""
        xmlroot = ET.fromstringlist(["<DNATCO_PIPE>", xmlText, "</DNATCO_PIPE>"])
        ET.indent(xmlroot, space="\t", level=0)
        with open(self.makeFileName('PROGRAMXML'), 'w', encoding="utf-8") as programXML:
            CCP4Utils.writeXML(programXML, ET.tostring(xmlroot))

        if bool(self.container.controlParameters.GENERATE_RESTRAINTS):
            dnatco_obj = self.dnatco2 if bool(self.container.controlParameters.TOGGLE_XYZIN2) else self.dnatco1
            self.container.outputData.RESTRAINTS.annotation.set(dnatco_obj.container.outputData.RESTRAINTS.annotation)
            restraintsPathJob = str(dnatco_obj.container.outputData.RESTRAINTS.fullPath)
            restraintsPathPipeline = str(self.container.outputData.RESTRAINTS.fullPath)
            if os.path.isfile(restraintsPathJob):
                shutil.copyfile(restraintsPathJob, restraintsPathPipeline)

        self.reportStatus(CPluginScript.SUCCEEDED)
