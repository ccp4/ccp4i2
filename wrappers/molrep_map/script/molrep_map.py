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
import subprocess

# import coot_headless_api as chapi


class molrep_map(CPluginScript):

    # Where this plugin will appear on the gui
    TASKMODULE = 'molecular_replacement'
    # A short title for gui menu
    # Task name - should be same as class name
    TASKNAME = 'molrep_map'
    TASKVERSION = 0.0                                     # Version of this plugin
    WHATNEXT = ['prosmart_refmac']
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    RUNEXTERNALPROCESS = False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.jobDir = pathlib.Path(self.getWorkDirectory())
        self.originalModel = self.jobDir / 'Original.pdb'
        self.flippedModel = self.jobDir / 'Flipped.pdb'
        self.flippedMapPath = self.jobDir / 'Flipped.map'
        self.originalBlurredMapPath = self.jobDir / 'OriginalBlurred.map'
        self.flippedBlurredMapPath = self.jobDir / 'FlippedBlurred.map'
        self.originalTrimmedMapPath = self.jobDir / 'OriginalTrimmed.map'
        self.flippedTrimmedMapPath = self.jobDir / 'FlippedTrimmed.map'

    def makeCommandAndScript(self, container=None):
        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        mc = chapi.molecules_container_t(False)
        iMap = mc.read_ccp4_map(
            self.container.inputData.MAPIN.__str__(), False)
        iMapFlipped = mc.flip_hand(iMap)
        mc.write_map(iMapFlipped, self.flippedMapPath.__str__())        
        iMapDownSampled = mc.sharpen_blur_map(iMap, 50,  False)
        mc.write_map(iMapDownSampled, self.originalBlurredMapPath.__str__())
        iMapFlippedDown = mc.flip_hand(iMapDownSampled)
        mc.write_map(iMapFlippedDown, self.flippedBlurredMapPath.__str__())
        return CPluginScript.SUCCEEDED

    def startProcess(self, command=None, handler=None, **kw):
        self.makeCommandAndScript()
        jobDir = pathlib.Path(self.getWorkDirectory())

        with open(jobDir / "molrep_com.txt", "w") as comFile:
            if self.container.controlParameters.MOLREP.NMON.isSet():                
                comFile.write(f"NMON {self.container.controlParameters.MOLREP.NMON.__str__()}\n")
            if self.container.controlParameters.MOLREP.NP.isSet():                
                comFile.write(f"NP {self.container.controlParameters.MOLREP.NP.__str__()}\n")

        with open(jobDir / "mapextend_com.txt", "w") as comFile:
            comFile.write("BORDER 15\n")
            
        firstHandDir = self.jobDir / 'FirstHand'
        firstHandDir.mkdir()
        with open(jobDir / "molrep_com.txt") as comIn:
            subprocess.run(['molrep', 
                            '-f', self.originalBlurredMapPath, 
                            '-m', self.container.inputData.XYZIN.fullPath.__str__(), 
                            '-i'],
                           cwd=firstHandDir.__str__(), 
                           stdin=comIn)

        outputCoordPath = firstHandDir / "molrep.pdb"

        if outputCoordPath.exists():
            shutil.copyfile(outputCoordPath, self.originalModel)
            with open(jobDir/"mapextend_com.txt", "r") as comIn:
                subprocess.run(['mapmask', 
                                'MAPIN', self.originalBlurredMapPath.__str__(), 
                                'XYZIN', outputCoordPath.__str__(),
                                'MAPOUT', self.originalTrimmedMapPath.__str__()],
                            cwd=firstHandDir.__str__(),
                            stdin=comIn)

        secondHandDir = self.jobDir / 'SecondHand'
        secondHandDir.mkdir()
        with open(jobDir / "molrep_com.txt") as comIn:
            subprocess.run(['molrep', 
                            '-f', self.flippedBlurredMapPath.__str__(), 
                            '-m', self.container.inputData.XYZIN.fullPath.__str__(), 
                            '-i'],
                           cwd=secondHandDir.__str__(), 
                           stdin=comIn)

        outputCoordPath = secondHandDir / "molrep.pdb"
        if outputCoordPath.exists():
            shutil.copyfile(outputCoordPath, self.flippedModel)
            with open(jobDir/"mapextend_com.txt", "r") as comIn:
                subprocess.run(['mapmask', 
                                'MAPIN', self.flippedMapPath.__str__(), 
                                'XYZIN', outputCoordPath.__str__(),
                                'MAPOUT', self.flippedTrimmedMapPath.__str__()],
                                cwd=secondHandDir.__str__(),
                                stdin=comIn)

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        if self.originalModel.exists():
            self.container.outputData.ORIGINALMODEL.setFullPath(self.originalModel.__str__())
            self.container.outputData.ORIGINALMODEL.annotation = "Model placed in original map"
            
        if self.flippedModel.exists():
            self.container.outputData.FLIPPEDMODEL.setFullPath(self.flippedModel.__str__())
            self.container.outputData.FLIPPEDMODEL.annotation = "Model placed in flipped map"
            
        if self.originalTrimmedMapPath.exists():
            self.container.outputData.ORIGINALTRIMMEDMAP.setFullPath(self.originalTrimmedMapPath.__str__())
            self.container.outputData.ORIGINALTRIMMEDMAP.annotation = "Original map trimmed to model"
            
        if self.flippedTrimmedMapPath.exists():
            self.container.outputData.FLIPPEDTRIMMEDMAP.setFullPath(self.flippedTrimmedMapPath.__str__())
            self.container.outputData.FLIPPEDTRIMMEDMAP.annotation = "Flipped map trimmed to model"
            
        if self.originalBlurredMapPath.exists():
            self.originalBlurredMapPath.unlink()
            
        if self.flippedBlurredMapPath.exists():
            self.flippedBlurredMapPath.unlink()
            
        return CPluginScript.SUCCEEDED


# ====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module


class testmolrep_map(unittest.TestCase):

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
        wrapper = molrep_map(
            parent=QTAPPLICATION(), name='molrep_map_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testmolrep_map)
    return suite


def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
