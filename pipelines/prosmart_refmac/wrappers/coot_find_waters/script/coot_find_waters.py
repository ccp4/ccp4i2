import os
import pathlib
from lxml import etree
from core.CCP4PluginScript import CPluginScript
from core.CCP4ModelData import CPdbDataFile

class coot_find_waters(CPluginScript):
    
    TASKMODULE = 'model_building'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Find waters with Coot API'     # A short title for gui menu
    TASKNAME = 'coot_find_waters'  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    
    def addWaters(self):
        outFormat = "cif" if self.container.inputData.XYZIN.isMMCIF() else "pdb"
        oldFullPath = pathlib.Path(str(self.container.outputData.XYZOUT.fullPath))
        if outFormat == "cif":
            self.container.outputData.XYZOUT.setFullPath(str(oldFullPath.with_suffix('.cif')))
            self.container.outputData.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_MMCIF)

        xyzoutPath = str(self.container.outputData.XYZOUT.fullPath)
        import chapi
        mc = chapi.molecules_container_py(True)
        mc.set_use_gemmi(True)
        imol = mc.read_pdb(str(self.container.inputData.XYZIN.fullPath))
        imap = mc.read_mtz(str(self.container.inputData.FPHIIN.fullPath),"F","PHI","",False,False)
        #mc.set_add_waters_water_to_protein_distance_lim_min(2.4) #Seemingly not in CCP4 release (nor Moorhen)
        #mc.set_add_waters_water_to_protein_distance_lim_min(3.4)
        nwaters = mc.add_waters(imol,imap)
        mc.write_coordinates(imol,xyzoutPath)
        return nwaters

    def process(self, command=None, handler=None, **kw):
        nwaters = self.add_waters()
        if os.path.exists(self.container.outputData.XYZOUT.__str__()): status = CPluginScript.SUCCEEDED
        from core import CCP4File
        root = etree.Element('coot_find_waters')

        nWatersElement = etree.SubElement(root, "WatersFound")
        nWatersElement.text = str(nwaters)
    
        self.container.outputData.XYZOUT.subType = 1
        
        f = CCP4File.CXmlDataFile(fullPath=self.makeFileName('PROGRAMXML'))
        f.saveFile(root)

        self.reportStatus(status)
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
