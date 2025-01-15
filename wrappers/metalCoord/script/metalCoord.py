"""
    metalCoord.py: CCP4 GUI Project
    Martin Maly, MRC-LMB
"""

import os
import json
from xml.etree import ElementTree as ET
from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import *
from core import CCP4Utils
from wrappers.servalcat.script.json2xml import json2xml


class metalCoord(CPluginScript):
    
    TASKMODULE = 'wrappers'        # Where this plugin will appear on gui
    TASKNAME = 'metalCoord'        # Task name - should be same as class name
    TASKVERSION = 0.1              # Version of this plugin
    TASKCOMMAND = 'metalCoord'     # The command to run the executable
    MAINTAINER = 'martin.maly@soton.ac.uk'

    ERROR_CODES = { 201 : { 'description' : 'No output JSON file from metalCoord' },
                    202 : { 'description' : 'Log file does not report successful job completion' , 'severity' : SEVERITY_WARNING },
                    }


    def processInputFiles(self):
        return CPluginScript.SUCCEEDED


    def makeCommandAndScript(self):
        self.appendCommandLine(['--no-progress'])
        self.appendCommandLine(['stats'])
        self.appendCommandLine(['-p', str(self.container.inputData.XYZIN.fullPath)])
        if self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER.isSet():
            self.appendCommandLine(['-c', str(self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER)])
        if self.container.controlParameters.MINIMUM_SAMPLE_SIZE.isSet():
            self.appendCommandLine(['-m', str(self.container.controlParameters.MINIMUM_SAMPLE_SIZE)])
        if self.container.controlParameters.DISTANCE_THRESHOLD.isSet():
            self.appendCommandLine(['-d', str(self.container.controlParameters.DISTANCE_THRESHOLD)])
        if self.container.controlParameters.PROCRUSTES_DISTANCE_THRESHOLD.isSet():
            self.appendCommandLine(['-t', str(self.container.controlParameters.PROCRUSTES_DISTANCE_THRESHOLD)])
        if self.container.controlParameters.IDEAL_ANGLES:
            self.appendCommandLine(['--ideal-angles'])
        if self.container.controlParameters.SIMPLE:
            self.appendCommandLine(['--simple'])
        if self.container.controlParameters.USE_PDB:
            self.appendCommandLine(['--use-pdb'])
        self.appendCommandLine(['-l', self.container.inputData.LIGAND_CODE])
        self.outputJsonFilename = str(self.container.inputData.LIGAND_CODE) + ".json"
        self.outputJsonPath = os.path.join(self.getWorkDirectory(), self.outputJsonFilename)
        self.appendCommandLine(['-o', self.outputJsonFilename])


    def processOutputFiles(self):
        # sanity check that metalCoord finished successfully - from .json.status.json
        self.outputJsonStatusFilename = str(self.container.inputData.LIGAND_CODE) + ".json.status.json"
        self.outputJsonStatusPath = os.path.join(self.getWorkDirectory(), self.outputJsonStatusFilename)
        with open(self.outputJsonStatusPath, 'r') as file:
            data = json.load(file)
        if data['status'] != "Success":
            self.appendErrorReport(202)

        # sanity check that metalCoord finished successfully - from log
        ok = False
        logText = self.logFileText()
        pyListLogLines = logText.split("\n")
        for j, pyStrLine in enumerate(pyListLogLines):
            if "Report written" in pyStrLine:
                ok = True
        if not ok: self.appendErrorReport(202)

        # Load JSON file content
        if os.path.isfile(self.outputJsonPath):
            self.container.outputData.JSON.setFullPath(self.outputJsonPath)
            self.container.outputData.JSON.annotation = 'Full analysis for monomer ' + str(self.container.inputData.LIGAND_CODE)
            with open(self.outputJsonPath, 'r') as outputJsonFile:
                outputJsonText = outputJsonFile.read()
        else:
            self.appendErrorReport(201, str(self.container.outputData.JSON))
            return CPluginScript.FAILED

        # Convert JSON to external restraint keywords
        self.outputRestraintsPrefix = str(self.container.inputData.LIGAND_CODE) + "_restraints"
        self.outputRestraintsFilename = self.outputRestraintsPrefix + ".txt"
        self.outputRestraintsPathPrefix = os.path.join(self.getWorkDirectory(), self.outputRestraintsPrefix)
        self.outputRestraintsPath = os.path.join(self.getWorkDirectory(), self.outputRestraintsFilename)
        if self.container.controlParameters.SAVE_PDBMMCIF:
            stPath = str(self.container.inputData.XYZIN.fullPath)
        else:
            stPath = None
        from . import json2restraints
        json2restraints.main(
            jsonPaths=[self.outputJsonPath],
            stPath=stPath,
            outputPrefix=self.outputRestraintsPathPrefix,
            jsonEquivalentsPath=None,
            keep_links=bool(self.container.controlParameters.KEEP_LINKS))
        if os.path.isfile(self.outputRestraintsPath):
            self.container.outputData.RESTRAINTS.setFullPath(self.outputRestraintsPath)
            self.container.outputData.RESTRAINTS.annotation = 'Restraints for ' + str(self.container.inputData.LIGAND_CODE)

        # Convert JSON to program.xml for i2 report
        outputJsonStats = json.loads(outputJsonText)
        xmlText = json2xml(list(outputJsonStats), tag_name_subroot="site")
        self.xmlroot = ET.fromstringlist(["<METALCOORD>", xmlText, "</METALCOORD>"])
        ET.indent(self.xmlroot, space="\t", level=0)
        with open(self.makeFileName('PROGRAMXML'), 'w') as programXML:
            CCP4Utils.writeXML(programXML, ET.tostring(self.xmlroot))

        return CPluginScript.SUCCEEDED
