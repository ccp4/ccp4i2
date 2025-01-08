import os
import shutil
import unittest

from lxml import etree

from ...core import CCP4Utils
from ...core.CCP4File import CDataFile
from ...core.CCP4Modules import PROCESSMANAGER
from ...core.CCP4Modules import QTAPPLICATION
from ...core.CCP4PluginScript import CPluginScript


class prosmart_martin(CPluginScript):
    
    TASKMODULE = 'deprecated'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Simple prosmart for use in prosmart_refmac pipeline'     # A short title for gui menu
    TASKNAME = 'prosmart_martin'                                  # Task name - should be same as class name
    TASKCOMMAND = 'prosmart'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    
    '''
        def __init__(self,parent=None,name=None,workDirectory=''):
        CPluginScript. __init__(self,parent=parent,name=name)
        '''
    def processInputFiles(self):
        self.tempFile = os.path.join(self.workDirectory,'XYZIN_temp.pdb')
        shutil.copyfile(self.container.inputData.XYZIN.__str__(), self.tempFile)
        return CPluginScript.SUCCEEDED
    
    def processOutputFiles(self):
        directory,fileName = os.path.split(self.tempFile)
        fileRoot,ext = os.path.splitext(fileName)
        restraintFilePath = os.path.join(self.workDirectory.__str__(), fileRoot + '.txt')
        print('My name is'+restraintFilePath)
        self.container.outputData.RESTRAINTS=CDataFile(restraintFilePath)
        self.container.outputData.RESTRAINTS.annotation = 'Restraints  to ' + fileName
        
        xmlPath = self.makeFileName('PROGRAMXML')
        try:
            xmlRoot= CCP4Utils.openFileToEtree(xmlPath)
        except:
            xmlRoot = etree.Element('PROSMART')
        htmlNode = etree.SubElement(xmlRoot,'htmlPath')
        htmlNode.text = 'ProSMART_Results.html'
        with open(xmlPath,'w') as xmlFile:
            xmlFile.write(etree.tostring(xmlRoot,pretty_print=True))
        
        #try:
        #    restraintSrcPath = #os.path.join(self.workDirectory.__str__(),'Output_Files','Restraints',str(self.container.inputData.MODEL))
        #restraintDestPath = self.container.outputData.RESTRAINTS.__str__()
        #   shutil.copyfile(restraintSrcPath, restraintDestPath)
        #except:
        #    pass
        #self.container.outputData.RESTRAINTS = CDataFile
        #self.container.outputData.RESTRAINTS.annotation = 'Generated restraints'
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        self.appendCommandLine(['-p1',self.tempFile])
        self.appendCommandLine(['-p2',self.container.inputData.TARGET])
        self.appendCommandLine(['-o',self.workDirectory])
        return CPluginScript.SUCCEEDED


#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

class testprosmart_martin(unittest.TestCase):
    
    def setUp(self):
        # make all background jobs wait for completion
        # this is essential for unittest to work
        self.app = QTAPPLICATION()
        PROCESSMANAGER().setWaitForFinished(10000)
    
    def tearDown(self):
        PROCESSMANAGER().setWaitForFinished(-1)
    
    def test_1(self):
        wrapper = prosmart_martin(parent=QTAPPLICATION(),name='prosmart_martin_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testprosmart_martin)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
