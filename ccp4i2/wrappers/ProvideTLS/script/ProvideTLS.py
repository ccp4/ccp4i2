import unittest

from lxml import etree

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript
from ....core.CCP4ProcessManager import PROCESSMANAGER
from ....utils.QApp import QTAPPLICATION


class ProvideTLS(CPluginScript):

    TASKNAME = 'ProvideTLS'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template

    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
      CPluginScript. __init__(self,parent=parent,name=name)
    '''
    
    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        
        self.checkOutputData()

        with open(self.container.outputData.TLSFILE.fullPath.__str__(),"w") as myFile:
            myFile.write(self.container.controlParameters.TLSTEXT.__str__() )

        root = etree.Element('ProvideTLSOutput')
        tlsElement = etree.SubElement(root,'TLSProvided')
        tlsElement.text = self.container.controlParameters.TLSTEXT.__str__()
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
            CCP4Utils.writeXML(xmlFile,etree.tostring(root,pretty_print=True))

        self.reportStatus(CPluginScript.SUCCEEDED)

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

class testProvideTLS(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     wrapper = ProvideTLS(parent=QTAPPLICATION(),name='ProvideTLS_test1')
     wrapper.container.loadDataFromXml()

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testProvideTLS)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
