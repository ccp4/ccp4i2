from __future__ import print_function

from lxml import etree

from core.CCP4PluginScript import CPluginScript
from PySide6 import QtCore

class coot_find_waters(CPluginScript):
    
    TASKMODULE = 'model_building'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Find waters with coot'     # A short title for gui menu
    TASKNAME = 'coot_find_waters'  # Task name - should be same as class name
    TASKCOMMAND = 'coot'
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    
    
    def makeCommandAndScript(self):
        import os
        cootScriptPath = os.path.join(self.workDirectory,'script.py')
        self.appendCommandLine(['--no-state-script','--no-graphics','--python','--pdb',self.container.inputData.XYZIN.fullPath,'--script',cootScriptPath])

        cootScript = open(cootScriptPath,"w")
        cootScript.write("make_and_draw_map(r'" + str(self.container.inputData.FPHIIN.fullPath)+"', 'F', 'PHI', 'PHI', 0, 0)\n")
        cootScript.write("set_ligand_water_to_protein_distance_limits("+str(self.container.controlParameters.MINDIST)+","+str(self.container.controlParameters.MAXDIST)+")\n")
        cootScript.write("execute_find_waters_real(1,0,0,"+str(self.container.controlParameters.THRESHOLD)+")\n")
        cootScript.write("write_pdb_file(0,r'"+str(self.container.outputData.XYZOUT.fullPath)+"')\n")
        cootScript.write("coot_real_exit(0)\n")
        cootScript.close()
        
        return CPluginScript.SUCCEEDED


    
    def processOutputFiles(self):
        import os
        status = CPluginScript.FAILED
        if os.path.exists(self.container.outputData.XYZOUT.__str__()): status = CPluginScript.SUCCEEDED
        print('coot_find_waters.handleFinish',self.container.outputData.XYZOUT, os.path.exists(self.container.outputData.XYZOUT.__str__()))
        # Create a trivial xml output file
        from core import CCP4File
        root = etree.Element('coot_find_waters')

        cootlines = open(self.makeFileName('LOG')).readlines()
        for line in cootlines:
            if line.startswith('INFO:: found'):
                nWatersElement = etree.SubElement(root,'WatersFound')
                lineElements = line.split()
                nWatersElement.text = lineElements[2]
    
        self.container.outputData.XYZOUT.subType = 1
        
        f = CCP4File.CXmlDataFile(fullPath=self.makeFileName('PROGRAMXML'))
        f.saveFile(root)
        return status



#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testcoot_find_waters(unittest.TestCase):
    
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
        wrapper = coot_find_waters(parent=QTAPPLICATION(),name='coot_find_waters_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testcoot_find_waters)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
