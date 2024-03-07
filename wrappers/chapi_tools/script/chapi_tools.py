from __future__ import print_function
import unittest

# from lxml import etree
from xml.etree import ElementTree as ET

from core.CCP4PluginScript import CPluginScript
from core import CCP4ModelData
import pathlib
import json
import shutil
import chapi

# import coot_headless_api as chapi


class chapi_tools(CPluginScript):

    # Where this plugin will appear on the gui
    TASKMODULE = 'model_building'
    # A short title for gui menu
    TASKTITLE = 'Access the coot toolbox through coot headless api'
    # Task name - should be same as class name
    TASKNAME = 'chapi_tools'
    TASKVERSION = 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    RUNEXTERNALPROCESS = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.outputJsonPath = pathlib.Path(self.getWorkDirectory()) / "output.json"

    def makeCommandAndScript(self):
        import os
        self.dropDir = pathlib.Path(self.workDirectory) / 'COOT_FILE_DROP'
        self.dropDir.mkdir(parents=True, exist_ok=True)
        return CPluginScript.SUCCEEDED
    
    def startProcess(self, command=None, handler=None, **kw):
        if self.container.controlParameters.STARTPOINT.__str__() == "FIT_LIGAND":
            self.fit_ligand()
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        with open(self.outputJsonPath, 'r') as outputJsonFile:
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

                        fullPath = pathlib.Path(outputFile['filePath'])
                        # We are an mmcif organisation, rename file if it is one, and convert to crate 
                        # an mmcif version (using gemmi) if not
                        if isinstance(newOutput, CCP4ModelData.CPdbDataFile) and fullPath.suffix == '.pdb':
                            cifFullPath = fullPath.with_suffix('.cif')
                            if newOutput.isMMCIF():
                                shutil.copy(fullPath.__str__(), cifFullPath.__str__())
                                newOutput.setFullPath(cifFullPath.__str__())
                            else:
                                gemmiStructure = gemmi.read_structure(fullPath.__str__())
                                gemmiStructure.make_mmcif_document().write_file(cifFullPath.__str__())
                                outputDataList.append(outputDataList.makeItem())
                                cifOutput = outputDataList[-1]
                                cifOutput.setFullPath(cifFullPath.__str__())
                                cifOutput.annotation = f"CIF copy of {outputFile['annotation']}"
                                cifOutput.setContentFlag(reset=True)
                                
        return CPluginScript.SUCCEEDED

    def fit_ligand(self):
        try:
            XYZIN_0, FPHIIN_0, DICTIN_0, TLC, MAXCOPIES  = [self.container.inputData.XYZIN[0].fullPath.__str__(),
                                                            self.container.inputData.FPHIIN[0].fullPath.__str__(),
                                                            self.container.inputData.DICTIN[0].fullPath.__str__(),
                                                            self.container.controlParameters.FIT_LIGAND.TLC,
                                                            self.container.controlParameters.FIT_LIGAND.MAXCOPIES,
                                                            ]
            mc = chapi.molecules_container_t(False)
            print({"mc":mc})
            iMol = mc.read_pdb(XYZIN_0)
            print({"iMol":iMol})
            iMap = mc.read_mtz(FPHIIN_0, 'F', 'PHI', '', False, False)
            print({"iMap":iMap})
            mc.set_imol_refinement_map(iMap)
            readDictResult = mc.import_cif_dictionary(DICTIN_0, iMol)
            print({"readDictResult":readDictResult})
            iDictMol = mc.get_monomer_from_dictionary(TLC, iMol, False)
            print({"iDictMol":iDictMol})
            solutions = mc.fit_ligand(iMol, iMap, iDictMol, 1.0, True, 30)
            print({"solutions":solutions})
            nSolutions = solutions.size()
            print({ "nSolutions":nSolutions })
            solutionMols = []
            for iSolution in range(nSolutions):
                solutionMols.append({ "imol": solutions.get(iSolution).imol, 
                                     "cluster_idx": solutions.get(iSolution).cluster_idx })
            
            print({ "solutionsMol":solutionMols })
            representedClusters = []
            filteredSolutions = []
            for solution in solutionMols:
                if not solution.cluster_idx in representedClusters:
                    representedClusters.append(solution.cluster_idx)
                    filteredSolutions.append(solution.imol)                
            
            print({ "representedClusters":representedClusters })
            print({ "filteredSolutions":filteredSolutions })
            useSolutions = []
            if len(filteredSolutions) > MAXCOPIES:
                useSolutions = filteredSolutions[:MAXCOPIES-1]
            else:
                useSolutions = filteredSolutions
            
            print({"useSolutions":useSolutions}, ":".join(useSolutions))
            mergeResult = mc.merge_molecules(iMol,  ":".join(useSolutions))
            print({ "mergeResults":mergeResult })
            outputFilePath = pathlib.Path(self.getWorkDirectory(), "result.pdb")
            print({ "outputFilePath":outputFilePath })
            writeResult = mc.writePDBASCII(iMol, outputFilePath.__str__())
            print({ "writeResult":writeResult })
            result = { "XYZOUT": [{ 
                "filePath": outputFilePath, 
                "annotation": "Model with ligands fit" 
                }] }
            with open (self.outputJsonPath, "w") as outputJson:
                outputJson.write(json.dumps(result))
            return CPluginScript.SUCCEEDED
        
        except Exception as err:
            return CPluginScript.FAILED

# ====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module


class testchapi_tools(unittest.TestCase):

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
        wrapper = chapi_tools(
            parent=QTAPPLICATION(), name='chapi_tools_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testchapi_tools)
    return suite


def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
