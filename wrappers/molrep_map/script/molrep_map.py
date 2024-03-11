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
        self.xmlRoot = ET.Element('molrep_map')

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
                comFile.write(
                    f"NMON {self.container.controlParameters.MOLREP.NMON.__str__()}\n")
            if self.container.controlParameters.MOLREP.NP.isSet():
                comFile.write(
                    f"NP {self.container.controlParameters.MOLREP.NP.__str__()}\n")

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
                           stdin=comIn,
                           stdout=open(firstHandDir / 'log.txt', 'w')
                           )

        outputDocPath = firstHandDir / "molrep.doc"
        if outputDocPath.exists():
            originalHandResults = ET.SubElement(self.xmlRoot,'Original')
            scrapedResults = self.scrapeXML(outputDocPath.__str__())
            originalHandResults.append(scrapedResults)

        outputCoordPath = firstHandDir / "molrep.pdb"
        if outputCoordPath.exists():
            shutil.copyfile(outputCoordPath, self.originalModel)
            with open(jobDir/"mapextend_com.txt", "r") as comIn:
                subprocess.run(['mapmask',
                                'MAPIN', self.originalBlurredMapPath.__str__(),
                                'XYZIN', outputCoordPath.__str__(),
                                'MAPOUT', self.originalTrimmedMapPath.__str__()],
                               cwd=firstHandDir.__str__(),
                               stdin=comIn
                               )

        secondHandDir = self.jobDir / 'SecondHand'
        secondHandDir.mkdir()
        with open(jobDir / "molrep_com.txt") as comIn:
            subprocess.run(['molrep',
                            '-f', self.flippedBlurredMapPath.__str__(),
                            '-m', self.container.inputData.XYZIN.fullPath.__str__(),
                            '-i'],
                           cwd=secondHandDir.__str__(),
                           stdin=comIn,
                           stdout=open(secondHandDir / 'log.txt', 'w')
                           )

        outputDocPath = secondHandDir / "molrep.doc"
        if outputDocPath.exists():
            flippedHandResults = ET.SubElement(self.xmlRoot,'Flipped')
            scrapedResults = self.scrapeXML(outputDocPath.__str__())
            flippedHandResults.append(scrapedResults)

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
            self.container.outputData.ORIGINALMODEL.setFullPath(
                self.originalModel.__str__())
            self.container.outputData.ORIGINALMODEL.annotation = "Model placed in original map"

        if self.flippedModel.exists():
            self.container.outputData.FLIPPEDMODEL.setFullPath(
                self.flippedModel.__str__())
            self.container.outputData.FLIPPEDMODEL.annotation = "Model placed in flipped map"

        if self.originalTrimmedMapPath.exists():
            self.container.outputData.ORIGINALTRIMMEDMAP.setFullPath(
                self.originalTrimmedMapPath.__str__())
            self.container.outputData.ORIGINALTRIMMEDMAP.annotation = "Original map trimmed to model"

        if self.flippedTrimmedMapPath.exists():
            self.container.outputData.FLIPPEDTRIMMEDMAP.setFullPath(
                self.flippedTrimmedMapPath.__str__())
            self.container.outputData.FLIPPEDTRIMMEDMAP.annotation = "Flipped map trimmed to model"

        if self.originalBlurredMapPath.exists():
            self.originalBlurredMapPath.unlink()

        if self.flippedBlurredMapPath.exists():
            self.flippedBlurredMapPath.unlink()

        with open(self.makeFileName('PROGRAMXML'), 'w') as outputXml:
            ET.indent(self.xmlRoot)
            outputXml.write(ET.tostring(self.xmlRoot).decode('utf-8'))
            
        return CPluginScript.SUCCEEDED

    def scrapeXML(self, docFileName):
        from core import CCP4Utils
        titles = []
        status = 0
        results = ET.Element('MolrepResult')
        tf = ET.Element('MR_TF')
        results.append(tf)
        for key, value in [['err_level', '0'],
                           ['err_message', 'normal termination'],
                           ['n_solution', '1'],
                           ['mr_score', '0.0000']]:

            e = ET.Element(key)
            e.text = value
            tf.append(e)

        '''
      table = [ '<?xml version="1.0" encoding="ASCII" standalone="yes"?>' ]
      table.append( "<MolrepResult>" )
      table.append( "<MR_TF>" )
      table.append( "Error <err_level>0</err_level>" )
      table.append( "Message <err_message>normal termination</err_message>" )
      table.append( "nmon_solution <n_solution>1</n_solution>" )
      table.append( "mr_score <mr_score>0.0000</mr_score>" )
      table.append( "</MR_TF>" )
      table.append( " <RFpeaks>" )
      '''

        docfileText = CCP4Utils.readFile(docFileName)
        docfileList = docfileText.split('\n')

        for line in docfileList:
            # print 'line in docfile',status,line
            if status == 0:
                lstr = line.strip()
                lst1 = "--- Translation function ---"
                lst2 = "--- phased translation function ---"
                if lstr == lst1 or lstr == lst2:
                    status = 1

            elif status == 1:
                if line.strip().startswith('RF '):
                    titles = line.replace("(", " ").replace(
                        ")", "").replace("/", "_").split()
                    # print 'titles',titles
                    rf = ET.Element('RFpeaks')
                    results.append(rf)

                else:
                    words = line.replace("(", " ").replace(
                        ")", "").replace("-", " -").split()
                    if len(words) == len(titles):
                        try:
                            for i in (0, 1):
                                ii = int(words[i])
                            for i in range(2, len(words)):
                                ii = float(words[i])
                            peak = ET.Element('RFpeak')
                            for key, value in zip(titles, words):
                                # print ' key,value', key,value
                                e = ET.Element(key)
                                e.text = str(float(value))
                                peak.append(e)
                        except:
                            pass
                        else:
                            rf.append(peak)
        return results

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
