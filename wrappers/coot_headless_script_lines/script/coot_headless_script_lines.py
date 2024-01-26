from __future__ import print_function

#from lxml import etree
from xml.etree import ElementTree as ET

from core.CCP4PluginScript import CPluginScript
import pathlib
import uuid
import sys
import os
import importlib
import gemmi

import sys
import platform
if platform.platform().startswith('macOS'):
    sys.path.insert(0,'/Applications/ccp4-9.0/coot_py3/Frameworks/Python.framework/Versions/3.9/lib/python3.9/site-packages')
elif platform.platform().startswith('Linux'):
    sys.path.insert(0,str(pathlib.Path(os.environ['CCP4']) / 'coot_py3' / 'lib' / 'python3.9' / 'site-packages'))
os.environ.setdefault('COOT_PREFIX', str(pathlib.Path(os.environ['CCP4']) / 'coot_py3'))
os.environ.setdefault('COOT_DATA_DIR', str(pathlib.Path(os.environ['CCP4']) / 'coot_py3' / 'share' / 'coot'))

import coot_headless_api as chapi

class coot_headless_script_lines(CPluginScript):
    
    TASKMODULE = 'model_building'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Execute arbitrary script code within Coot'     # A short title for gui menu
    TASKNAME = 'coot_headless_script_lines'                                  # Task name - should be same as class name
    TASKVERSION= 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    RUNEXTERNALPROCESS = False
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.scriptFileRoot = f'script_{uuid.uuid4()}'
        self.scriptFilePath = pathlib.Path(self.workDirectory) / f'{self.scriptFileRoot}.py'
        self.mc = chapi.molecules_container_t(True)

    def makeCommandAndScript(self):
        import os
        self.dropDir = pathlib.Path(self.workDirectory) / 'COOT_FILE_DROP'
        self.dropDir.mkdir(parents=True, exist_ok=True)

        with open(self.scriptFilePath, 'w') as scriptFile:
            scriptFile.write(self.container.controlParameters.SCRIPT.__str__())

        self.appendCommandLine([self.scriptFilePath.__str__()])

        return CPluginScript.SUCCEEDED

    def startProcess(self, command, **kw):
        sys.path.insert(0, pathlib.Path(self.workDirectory).absolute().__str__())
        importedModule = importlib.import_module(self.scriptFileRoot)
        self.result=importedModule.exec(self, self.mc)
        return CPluginScript.SUCCEEDED
    
    def processOutputFiles(self):

        return CPluginScript.SUCCEEDED
    
def fit_ligands(plugin, mc, **kwargs):
    inputData = plugin.container.inputData
    mc.import_cif_dictionary(str(inputData.DICTIN[0]), -999999)
    structure = gemmi.read_structure(str(inputData.XYZIN[0]))
    gemmifiedCoordPath = str(pathlib.Path(str(plugin.workDirectory)) / 'new.cif')
    structure.make_mmcif_document().write_file(gemmifiedCoordPath)
    MolHandle_1 = mc.read_pdb(gemmifiedCoordPath)
    print(f"Read model {MolHandle_1}")
    sys.stdout.flush()
    MapHandle_1 = mc.read_mtz(str(inputData.FPHIIN[0]), 
                              'F', 'PHI', '', False, False)
    print(f"Read map {MapHandle_1}")
    sys.stdout.flush()
    imol_lig = mc.get_monomer('DRG')
    if mc.is_valid_model_molecule(MolHandle_1):
        if mc.is_valid_map_molecule(MapHandle_1):
            solutions = mc.fit_ligand(MolHandle_1, MapHandle_1, imol_lig, 1.0, True, 30)
            print("Found", len(solutions), "solutions")
            if solutions:
                sl = ":".join([f'{iSol}' for ii, iSol in enumerate(solutions)])
                mc.merge_molecules(MolHandle_1, sl)
    mc.write_coordinates(MolHandle_1, str(pathlib.Path() / 'COOT_FILE_DROP' / 'output.pdb'))
    return(True)

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testcoot_headless_script_lines(unittest.TestCase):
    
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
        wrapper = coot_headless_script_lines(parent=QTAPPLICATION(),name='coot_headless_script_lines_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testcoot_headless_script_lines)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
