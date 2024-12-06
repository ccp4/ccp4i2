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
        ##import os
        ##import shutil
        # Use temp input filename from which prosmart takes output restraints filename
        ##self.tempFile = os.path.splitext(str(self.container.outputData.RESTRAINTS))[0]+'_TARGET.pdb'
        ##print('prosmart tempFile',self.tempFile)
        ##shutil.copyfile(self.container.inputData.TARGET_MODEL.__str__(), self.tempFile)
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

        '''xmlPath = self.makeFileName('PROGRAMXML')
        from lxml import etree
        xmlRoot = etree.Element('PROSMART')
        xmlString = etree.tostring(xmlRoot,pretty_print=True)
        xmlFile=open( xmlPath,'w')
        xmlFile.write( xmlString )
        xmlFile.close()'''

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
            CCP4Utils.writeXML(programXML, ET.tostring(self.xmlroot)) # CCP4Utils.writeXML(programXML, etree.tostring(self.xmlroot, pretty_print=True))

        return CPluginScript.SUCCEEDED

#======================================================
# PLUGIN TESTS
# See Python documentation on unittest module
"""
import unittest

class testprosmart(unittest.TestCase):
    
    def setUp(self):
        from core import CCP4Modules
        self.app = CCP4Modules.QTAPPLICATION()
        # make all background jobs wait for completion
        # this is essential for unittest to work
        CCP4Modules.PROCESSMANAGER().setWaitForFinished(10000)
    
    def tearDown(self):
        from core import CCP4Modules
        CCP4Modules.PROCESSMANAGER().setWaitForFinished(-1)
    
    def test_1(self):
        from core import CCP4Modules, CCP4Utils
        import os
        
        workDirectory = CCP4Utils.getTestTmpDir()
        # this needs to agree with name attribute below
        logFile = os.path.join(workDirectory,'prosmart_test1.log')
        # Delete any existing log file
        if os.path.exists(logFile): os.remove(logFile)
        
        self.wrapper = prosmart(parent=CCP4Modules.QTAPPLICATION(),name='prosmart_test1',workDirectory=workDirectory)
        self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','prosmart','test_data','prosmart_test1.data.xml'))
        
        self.wrapper.setWaitForFinished(1000000)
        pid = self.wrapper.process()
        self.wrapper.setWaitForFinished(-1)
        if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())
    #self.assertTrue(os.path.exists(logFile),'No log file found')
    
    def test_2(self):
        from core import CCP4Modules, CCP4Utils
        import os
        
        workDirectory = CCP4Utils.getTestTmpDir()
        # this needs to agree with name attribute below
        logFile = os.path.join(workDirectory,'prosmart_test2.log')
        # Delete any existing log file
        if os.path.exists(logFile): os.remove(logFile)
        
        self.wrapper = prosmart(parent=CCP4Modules.QTAPPLICATION(),name='prosmart_test2',workDirectory=workDirectory)
        self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','prosmart','test_data','prosmart_test2.data.xml'))
        
        self.wrapper.setWaitForFinished(1000000)
        pid = self.wrapper.process()
        self.wrapper.setWaitForFinished(-1)
        if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())
#self.assertTrue(os.path.exists(logFile),'No log file found')


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testprosmart)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
"""
