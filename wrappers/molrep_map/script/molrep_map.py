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
    TASKCOMMAND = 'molrep'
    RUNEXTERNALPROCESS = True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.blurredMapPath = (pathlib.Path(
            self.getWorkDirectory()) / 'Blurred.map').absolute().__str__()
        self.blurredFlippedMapPath = (pathlib.Path(
            self.getWorkDirectory()) / 'BlurredFlipped.map').absolute().__str__()

    def makeCommandAndScript(self, container=None):
        return CPluginScript.SUCCEEDED

    def processInputFiles(self):
        mc = chapi.molecules_container_t(False)
        iMap = mc.read_ccp4_map(
            self.container.inputData.MAPIN.__str__(), False)
        iMapDownSampled = mc.sharpen_blur_map(iMap, 50,  False)
        print({"blurredMapPath": self.blurredMapPath})
        mc.write_map(iMapDownSampled, self.blurredMapPath)
        iMapFlipped = mc.flip_hand(iMapDownSampled)
        mc.write_map(iMapFlipped, self.blurredFlippedMapPath)
        return CPluginScript.SUCCEEDED

    def startProcess(self, command=None, handler=None, **kw):
        self.makeCommandAndScript()
        jobDir = pathlib.Path(self.getWorkDirectory())
        with open(jobDir / "com.txt", "w") as comFile:
            comFile.write("NMON 1\n")

        firstHandDir = jobDir / 'FirstHand'
        firstHandDir.mkdir()
        with open(jobDir / "com.txt") as comIn:
            subprocess.run(['molrep', '-f', self.blurredMapPath, '-m',
                            self.container.inputData.XYZIN.fullPath.__str__()],
                           cwd=firstHandDir.__str__(), stdin=comIn)
        outputCoordPath = firstHandDir/"molrep.pdb"
        if outputCoordPath.exists():
            subprocess.run(['mapmask', 'MAPIN', self.blurredMapPath, 'XYZIN',
                            outputCoordPath.fullPath.__str__(), 
                            'MAPOUT', (firstHandDir / 'trimmed.map').__str__()],
                           cwd=firstHandDir.__str__())

        secondHandDir = jobDir / 'SecondHand'
        secondHandDir.mkdir()
        with open(jobDir / "com.txt") as comIn:
            subprocess.run(['molrep', '-f', self.blurredFlippedMapPath, '-m',
                            self.container.inputData.XYZIN.fullPath.__str__()],
                           cwd=firstHandDir.__str__(), stdin=comIn)
        outputCoordPath = secondHandDir/"molrep.pdb"
        if outputCoordPath.exists():
            subprocess.run(['mapmask', 'MAPIN', self.blurredFlippedMapPath, 'XYZIN',
                            outputCoordPath.fullPath.__str__(), 
                            'MAPOUT', (secondHandDir / 'trimmed.map').__str__()],
                           cwd=secondHandDir.__str__())

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
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
