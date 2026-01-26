"""
    dnatco.py: CCP4 GUI Project
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

from pathlib import Path
import xml.etree.ElementTree as ET
from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import *
from core import CCP4Utils


class dnatco(CPluginScript):

    TASKMODULE = 'wrappers'  # Where this plugin will appear on gui
    TASKNAME = 'dnatco'      # Task name - should be same as class name
    TASKVERSION= 0.1         # Version of this plugin
    TASKCOMMAND = '/opt/ccp4-10/bin/dnatco'   # The command to run the executable
    MAINTAINER = 'martin.maly@mrc-lmb.cam.ac.uk'

    ERROR_CODES = { 201 : { 'description' : 'No output restraint file from DNATCO' },
                    202 : { 'description' : 'No output extended mmCIF file from DNATCO' },
                    203 : { 'description' : 'No output PDF report file from DNATCO' },
                    204 : { 'description' : 'Failed to dnatcoify structure' },
                    }


    def makeCommandAndScript(self):
        self.appendCommandLine(['--coords', str(self.container.inputData.XYZIN.fullPath)])
        self.appendCommandLine(['--outputDir', self.workDirectory])
        self.appendCommandLine(['--prefix', "dnatco"])
        self.appendCommandLine(['--extendedCIF'])
        # if self.container.inputData.HKLIN.isSet():
        #     self.appendCommandLine(['--reflns', str(self.container.inputData.HKLIN.fullPath)])
        if bool(self.container.controlParameters.GENERATE_RESTRAINTS):
            self.appendCommandLine(['--refmacRestraints'])
            if float(self.container.controlParameters.MAX_RMSD) != 0.5:
                self.appendCommandLine(['--restraintsRmsd', str(float(self.container.controlParameters.MAX_RMSD))])
            if float(self.container.controlParameters.RESTRAINTS_SIGMA) != 1.0:
                self.appendCommandLine(['--restraintsSigmaFactor', str(float(self.container.controlParameters.RESTRAINTS_SIGMA))])


    def processOutputFiles(self):
        # outputFilenamePrefix = str(os.path.splitext(
        #     os.path.basename(str(
        #         self.container.inputData.XYZIN.fullPath)))[0]).upper().replace("_1", "")  # TODO
        # outputFilenamePrefix = "1Q93"

        # check that DNATCO has finished successfully
        logText = self.logFileText()
        pyListLogLines = logText.split("\n")
        for j, pyStrLine in enumerate(pyListLogLines):
            if "Failed to dnatcoify structure" in pyStrLine:
                self.appendErrorReport(204)
                return CPluginScript.FAILED

        # outputCifFilename = outputFilenamePrefix + '_extended.cif'
        # outputCifPath = os.path.normpath(os.path.join(self.getWorkDirectory(), outputCifFilename))
        work_dir = Path(self.workDirectory)
        cif_file = "dnatco_extended.cif"
        if Path(str(cif_file)).is_file():
            self.appendErrorReport(202, str(cif_file))
            return CPluginScript.FAILED
        self.container.outputData.CIFOUT.setFullPath(str(cif_file))
        self.container.outputData.CIFOUT.annotation.set('Extended model (mmCIF format)')

        if bool(self.container.controlParameters.GENERATE_RESTRAINTS):
            # outputRestraintsFilename = outputFilenamePrefix + '_restraints_refmac.txt'
            # outputRestraintsPath = os.path.normpath(os.path.join(self.getWorkDirectory(), outputRestraintsFilename))
            restraint_file = "dnatco_restraints_refmac.txt"
            if Path(str(restraint_file)).is_file():
                outputRestraintsPath = str(restraint_file)
            else:
                outputRestraintsPath = ""
            if Path(outputRestraintsPath).is_file():
                self.container.outputData.RESTRAINTS.setFullPath(outputRestraintsPath)
                self.container.outputData.RESTRAINTS.annotation.set('DNATCO restraints')
            else:
                self.appendErrorReport(201, outputRestraintsPath)
                return CPluginScript.FAILED

        xmlText = ""
        xmlroot = ET.fromstringlist(["<DNATCO>", xmlText, "</DNATCO>"])
        ET.indent(xmlroot, space="\t", level=0)
        with open(self.makeFileName('PROGRAMXML'), 'w', encoding="utf-8") as programXML:
            CCP4Utils.writeXML(programXML, ET.tostring(xmlroot))

        return CPluginScript.SUCCEEDED
