"""
Martin Maly, MRC-LMB
"""

import os
import json
import xml.etree.ElementTree as ET
from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import SEVERITY_WARNING
from core import CCP4Utils
from wrappers.servalcat.script.json2xml import json2xml
from . import json2restraints


class metalCoord(CPluginScript):

    TASKMODULE = 'wrappers'        # Where this plugin will appear on gui
    TASKNAME = 'metalCoord'        # Task name - should be same as class name
    TASKVERSION = 0.2              # Version of this plugin
    TASKCOMMAND = 'metalCoord'     # The command to run the executable
    MAINTAINER = 'martin.maly@mrc-lmb.cam.ac.uk'

    ERROR_CODES = {
        201: {'description': 'No output JSON file from metalCoord'},
        202: {'description': 'Log file does not report successful job completion', 'severity': SEVERITY_WARNING},
    }

    def makeCommandAndScript(self):
        self.appendCommandLine('--no-progress')
        self.appendCommandLine('stats')
        self.appendCommandLine(['-p', self.container.inputData.XYZIN.fullPath])
        if self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER.isSet():
            if str(self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER) != "auto":
                self.appendCommandLine(['-c', self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER])
                option = f"COORD{int(str(self.container.controlParameters.MAXIMUM_COORDINATION_NUMBER)):02d}"
                if hasattr(self.container.coordination, option):
                    if getattr(self.container.coordination, option) != "auto":
                        self.appendCommandLine(['--cl', getattr(self.container.coordination, option)])
        if self.container.controlParameters.MINIMUM_SAMPLE_SIZE.isSet():
            self.appendCommandLine(['-m', self.container.controlParameters.MINIMUM_SAMPLE_SIZE])
        if self.container.controlParameters.DISTANCE_THRESHOLD.isSet():
            self.appendCommandLine(['-d', self.container.controlParameters.DISTANCE_THRESHOLD])
        if self.container.controlParameters.PROCRUSTES_DISTANCE_THRESHOLD.isSet():
            self.appendCommandLine(['-t', self.container.controlParameters.PROCRUSTES_DISTANCE_THRESHOLD])
        if self.container.controlParameters.IDEAL_ANGLES:
            self.appendCommandLine('--ideal-angles')
        if self.container.controlParameters.SIMPLE:
            self.appendCommandLine('--simple')
        if self.container.controlParameters.USE_PDB:
            self.appendCommandLine('--use-pdb')
        self.appendCommandLine(['-l', self.container.inputData.LIGAND_CODE])
        outputJsonFilename = str(self.container.inputData.LIGAND_CODE) + ".json"
        self.outputJsonPath = os.path.join(self.getWorkDirectory(), outputJsonFilename)
        self.appendCommandLine(['-o', outputJsonFilename])

    def processOutputFiles(self):
        # sanity check that metalCoord finished successfully - from .json.status.json
        outputJsonStatusFilename = str(self.container.inputData.LIGAND_CODE) + ".json.status.json"
        outputJsonStatusPath = os.path.join(self.getWorkDirectory(), outputJsonStatusFilename)
        with open(outputJsonStatusPath, encoding="utf-8") as file:
            data = json.load(file)
        if data['status'] != "Success":
            self.appendErrorReport(202)

        # sanity check that metalCoord finished successfully - from log
        for line in self.logFileText().splitlines():
            if "Report written" in line:
                break
        else:
            self.appendErrorReport(202)

        # Load JSON file content
        if os.path.isfile(self.outputJsonPath):
            self.container.outputData.JSON.setFullPath(self.outputJsonPath)
            self.container.outputData.JSON.annotation = 'Full analysis for monomer ' + str(self.container.inputData.LIGAND_CODE)
            with open(self.outputJsonPath, encoding="utf-8") as outputJsonFile:
                outputJsonText = outputJsonFile.read()
        else:
            self.appendErrorReport(201, str(self.container.outputData.JSON))
            return CPluginScript.FAILED

        # Convert JSON to external restraint keywords
        outputRestraintsPrefix = str(self.container.inputData.LIGAND_CODE) + "_restraints"
        outputRestraintsFilename = outputRestraintsPrefix + ".txt"
        outputRestraintsMmcifFilename = outputRestraintsPrefix + ".mmcif"
        outputRestraintsPathPrefix = os.path.join(self.getWorkDirectory(), outputRestraintsPrefix)
        outputRestraintsPath = os.path.join(self.getWorkDirectory(), outputRestraintsFilename)
        outputRestraintsMmcifPath = os.path.join(self.getWorkDirectory(), outputRestraintsMmcifFilename)
        if self.container.controlParameters.SAVE_PDBMMCIF:
            stPath = str(self.container.inputData.XYZIN.fullPath)
        else:
            stPath = None
        json2restraints.main(
            jsonPaths=[self.outputJsonPath],
            stPath=stPath,
            outputPrefix=outputRestraintsPathPrefix,
            jsonEquivalentsPath=None,
            keep_links=bool(self.container.controlParameters.KEEP_LINKS))
        if os.path.isfile(outputRestraintsPath):
            self.container.outputData.RESTRAINTS.setFullPath(outputRestraintsPath)
            self.container.outputData.RESTRAINTS.annotation = 'Restraints for ' + str(self.container.inputData.LIGAND_CODE)
        if self.container.controlParameters.SAVE_PDBMMCIF:
            if os.path.isfile(outputRestraintsMmcifPath) and stPath:
                self.container.outputData.XYZOUT.setFullPath(outputRestraintsMmcifPath)
                self.container.outputData.XYZOUT.annotation = 'Structure model with links from MetalCoord (mmCIF format)'

        # Convert JSON to program.xml for i2 report
        outputJsonStats = json.loads(outputJsonText)
        xmlText = json2xml(list(outputJsonStats), tag_name_subroot="site")
        xmlroot = ET.fromstringlist(["<METALCOORD>", xmlText, "</METALCOORD>"])
        ET.indent(xmlroot, space="\t", level=0)
        with open(self.makeFileName('PROGRAMXML'), 'w', encoding="utf-8") as programXML:
            CCP4Utils.writeXML(programXML, ET.tostring(xmlroot))

        return CPluginScript.SUCCEEDED
