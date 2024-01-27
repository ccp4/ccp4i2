from __future__ import print_function

#from lxml import etree
from xml.etree import ElementTree as ET

from core.CCP4PluginScript import CPluginScript
import pathlib
import uuid

#import coot_headless_api as chapi

class moorhen_node_tools(CPluginScript):
    
    TASKMODULE = 'model_building'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Access the coot toolbox through moorhen javascript'     # A short title for gui menu
    TASKNAME = 'moorhen_node_tools'                                  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
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
        self.appendCommandLine([str(pathlib.Path(__file__).parents[0] / 'script.js'), self.container.controlParameters.STARTPOINT.__str__()])
        return CPluginScript.SUCCEEDED
    
    def processOutputFiles(self):

        return CPluginScript.SUCCEEDED
    

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testmoorhen_node_tools(unittest.TestCase):
    
    def setUp(self):
        # make all background jobs wait for completion
        # this is essential for unittest to work
        from core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
        self.app = QTAPPLICATION()
        PROCESSMANAGER().setWaitForFinished(10000)
    
    def tearDown(self):
        from core.CCP4Modules import PROCESSMANAGER
        PROCESSMANAGER().setWaitForFinished(-1)
    
    def test_1(self):
        from core.CCP4Modules import QTAPPLICATION
        wrapper = moorhen_node_tools(parent=QTAPPLICATION(),name='moorhen_node_tools_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testmoorhen_node_tools)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
