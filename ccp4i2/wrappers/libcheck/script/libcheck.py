
from core.CCP4PluginScript import CPluginScript

  
class libcheck(CPluginScript):

    TASKMODULE = 'test'                               # Where this plugin will appear on the gui
    TASKTITLE = 'Manage dictionary libraries'             # A short title for gui menu
    TASKNAME = 'libcheck'                                  # Task name - should be same as class name
    TASKCOMMAND = 'libcheck'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
      
    def process(self):
      if self.container.controlParameters.RUN_MODE == 'MERGE':
        return self.mergeProcess()

    def mergeProcess(self):
      import os,copy
      errorFiles = []; copies = 0
      lastLib = self.container.inputData.DICTLIB.__str__()
      for indx in range(len(self.container.inputData.MERGELIST)):
        mergeFile = self.container.inputData.MERGELIST[indx]
        newLib = os.path.join(self.workDirectory,'tmp_merged_'+str(indx))
        self.appendCommandScript('_Y',clear=True)
        self.appendCommandScript('_FILE_L '+lastLib)
        self.appendCommandScript('_FILE_L2 '+mergeFile.__str__())
        self.appendCommandScript('_FILE_O '+newLib)
        self.appendCommandScript('_END')        
        status = self.startProcess(command=self.command,reportStatus=False,fileQualifier=str(copies))
        if status == CPluginScript.SUCCEEDED:
          # Beware libcheck sticks a 'lib' extension even if there is already one there
          lastLib = copy.deepcopy(newLib) + '.lib'
          copies += 1          
        else:
          errorFiles.append(mergeFile.__str__())
     
      if copies>0 and os.path.exists(lastLib):
        from core import CCP4Utils
        CCP4Utils.backupFile(self.container.inputData.DICTLIB.__str__(),delete=True)
        #print 'libcheck.process lastLib',lastLib,os.path.exists(lastLib),self.container.inputData.DICTLIB.__str__()
        os.rename(lastLib,self.container.inputData.DICTLIB.__str__())
        return CPluginScript.SUCCEEDED
      else:
        return CPluginScript.FAILED
     
#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testlibcheck(unittest.TestCase):

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
     wrapper = libcheck(parent=QTAPPLICATION(),name='libcheck_test1')
     wrapper.container.loadDataFromXml()
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testlibcheck)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
