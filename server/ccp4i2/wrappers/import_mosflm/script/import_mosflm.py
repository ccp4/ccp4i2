
from ccp4i2.core.CCP4PluginScript import CPluginScript

  
class import_mosflm(CPluginScript):

    TASKNAME = 'import_mosflm'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

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
        programXMLPath = self.makeFileName('PROGRAMXML')
        shutil.copyfile(self.container.inputData.MOSFLMXML.fullPath.__str__(), programXMLPath)
        self.reportStatus(CPluginScript.SUCCEEDED)

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testimport_mosflm(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    from ccp4i2.core.CCP4Modules import QTAPPLICATION,PROCESSMANAGER
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from ccp4i2.core.CCP4Modules import PROCESSMANAGER
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     from ccp4i2.core.CCP4Modules import QTAPPLICATION
     wrapper = import_mosflm(parent=QTAPPLICATION(),name='import_mosflm_test1')
     wrapper.container.loadDataFromXml()
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testimport_mosflm)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
