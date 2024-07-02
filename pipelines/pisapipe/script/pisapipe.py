from __future__ import print_function

from PySide6 import QtCore
from core.CCP4PluginScript import CPluginScript

  
class pisapipe(CPluginScript):

    TASKMODULE = 'test'                              # Where this plugin will appear on the gui
    TASKTITLE = ' Structure analysis with Pisa'     # A short title for gui menu
    TASKNAME = 'pisapipe'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template


    def process(self):
      #identifier = str(self.container.paramsHeader.project) +'_'+str(self.container.paramsHeader.jobNumber)
      identifier = str(self.container.inputData.PDBIN.dbFileId)
      btask = self.makePluginObject('pisa_analyse')
      btask.container.controlParameters.IDENTIFIER = identifier
      # set synchronous job control
      btask.setWaitForFinished(1000000)
      btask.container.inputData.PDBIN = self.container.inputData.PDBIN
      jobStatus = btask.process()
      if jobStatus == CPluginScript.FAILED:
          self.reportStatus(jobStatus)
          return
      self.xmlTask = self.makePluginObject('pisa_xml')
      self.xmlTask.container.controlParameters.IDENTIFIER = identifier
      self.connectSignal(self.xmlTask,'finished',self.pisa_xmlFinished)
      self.xmlTask.waitForFinished = -1
      self.xmlTask.doAsync=False
      self.xmlTask.process()

    @QtCore.Slot(dict)
    def pisa_xmlFinished(self, statusDict):
        status = statusDict['finishStatus']
        print('pisa_xmlFinished', status)
        if status == CPluginScript.FAILED: self.reportStatus(status)
        from core import CCP4Utils
        from lxml import etree
        self.xmlroot = etree.Element('pisapipe')
        xmlOfTask = CCP4Utils.openFileToEtree(self.xmlTask.makeFileName('PROGRAMXML'))
        self.xmlroot.append(xmlOfTask)
        with open(self.makeFileName('PROGRAMXML'),'w') as outputXMLFile:
            CCP4Utils.writeXML(outputXML,etree.tostring(self.xmlroot,pretty_print=True))
        self.reportStatus(CPluginScript.SUCCEEDED)
     
#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testpisa(unittest.TestCase):

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
     wrapper = pisapipe(parent=QTAPPLICATION(),name='pisa_test1')
     wrapper.container.loadDataFromXml()
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testpisa)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
