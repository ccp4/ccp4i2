from __future__ import print_function
import unittest

# from lxml import etree
from xml.etree import ElementTree as ET

from core.CCP4PluginScript import CPluginScript
import pathlib
import json

# import coot_headless_api as chapi


class moorhen_node_tools(CPluginScript):

    # Where this plugin will appear on the gui
    TASKMODULE = 'model_building'
    # A short title for gui menu
    TASKTITLE = 'Access the coot toolbox through moorhen javascript'
    # Task name - should be same as class name
    TASKNAME = 'moorhen_node_tools'
    TASKVERSION = 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    RUNEXTERNALPROCESS = True
    TASKCOMMAND = 'node'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def makeCommandAndScript(self):
        import os
        self.dropDir = pathlib.Path(self.workDirectory) / 'COOT_FILE_DROP'
        self.dropDir.mkdir(parents=True, exist_ok=True)
        inputJsonPath = pathlib.Path(self.workDirectory) / "input.json"
        self.appendCommandLine([str(pathlib.Path(__file__).parents[0] / 'script.js'),
                                inputJsonPath,
                                str(pathlib.Path(__file__).parents[0]),
                               self.getWorkDirectory()])
        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        containerContent = self.container.jsonString()
        inputJsonPath = pathlib.Path(self.workDirectory) / "input.json"
        with open(inputJsonPath, 'w') as inputJson:
            inputJson.write(containerContent)

    def processOutputFiles(self):
        outputJsonPath = pathlib.Path(self.workDirectory) / "output.json"
        with open(outputJsonPath, 'r') as outputJsonFile:
            outputJson = outputJsonFile.read()
            outputData = json.loads(outputJson)
            for outputType in self.container.outputData.contents():
                if outputType in outputData:
                    outputDataList = getattr(
                        self.container.outputData, outputType)
                    for outputFile in outputData[outputType]:
                        outputDataList.append(outputDataList.makeItem())
                        newOutput = outputDataList[-1]
                        newOutput.setFullPath(outputFile['filePath'])
                        newOutput.annotation = outputFile['annotation']
                        newOutput.setContentFlag(reset=True)
        return CPluginScript.SUCCEEDED


# ====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module


class testmoorhen_node_tools(unittest.TestCase):

    def setUp(self):
        # make all background jobs wait for completion
        # this is essential for unittest to work
        from core.CCP4Modules import QTAPPLICATION, PROCESSMANAGER
        self.app = QTAPPLICATION()
        PROCESSMANAGER().setWaitForFinished(10000)

    def tearDown(self):
        from core.CCP4Modules import PROCESSMANAGER
        PROCESSMANAGER().setWaitForFinished(-1)

    def test_1(self):
        from core.CCP4Modules import QTAPPLICATION
        wrapper = moorhen_node_tools(
            parent=QTAPPLICATION(), name='moorhen_node_tools_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testmoorhen_node_tools)
    return suite


def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
