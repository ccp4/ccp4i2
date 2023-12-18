
from core.CCP4PluginScript import CPluginScript
from core import CCP4Utils
from xml.etree import ElementTree as ET
  
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
        
        import shutil
        with open(self.container.outputData.TLSFILE.fullPath.__str__(),"w") as myFile:
            myFile.write(self.container.controlParameters.TLSTEXT.__str__() )
        
        root = ET.Element('ProvideTLSOutput')
        tlsElement = ET.SubElement(root,'TLSProvided')
        tlsElement.text = self.container.controlParameters.TLSTEXT.__str__()
        with open(self.makeFileName('PROGRAMXML'),'w') as xmlFile:
            ET.indent(root)
            CCP4Utils.writeXML(xmlFile,ET.tostring(root))

        self.reportStatus(CPluginScript.SUCCEEDED)

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testProvideTLS(unittest.TestCase):

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
     wrapper = ProvideTLS(parent=QTAPPLICATION(),name='ProvideTLS_test1')
     wrapper.container.loadDataFromXml()
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testProvideTLS)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
