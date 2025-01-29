import os
import unittest

from ....core import CCP4Utils
from ....core.CCP4PluginScript import CPluginScript
from ....core.CCP4ProcessManager import PROCESSMANAGER
from ....utils.QApp import QTAPPLICATION


class baverage(CPluginScript):

    TASKMODULE = 'wrappers'  # Where this plugin will appear on the gui
    TASKTITLE = 'Average B over main and side chain atoms'
    TASKNAME = 'baverage'
    TASKCOMMAND = 'baverage'
    TASKVERSION= 0.0
    COMTEMPLATEFILE = None

    def makeCommandAndScript(self):
      inp = self.container.inputData
      out = self.container.outputData
      par = self.container.controlParameters
      self.appendCommandLine(['XYZIN',inp.XYZIN.fullPath])
      self.appendCommandLine(['RMSTAB',out.RMSTAB.fullPath])
      self.appendCommandLine(['XYZOUT',out.XYZOUT.fullPath])

      try:
        if self.container.guiAdmin.jobTitle.isSet(): self.appendCommandScript("TITLE %s" % self.container.guiAdmin.jobTitle.__str__())
      except:
         pass
      if par.SET_BLIMIT:
         self.appendCommandScript("BLIMIT %f %f %f %f" % (par.BLIMIT_MC.start,par.BLIMIT_MC.end,par.BLIMIT_SC.start,par.BLIMIT_SC.end))
      if par.SET_BRANGE:
         self.appendCommandScript("BRANGE %f" % par.BRANGE)
      self.appendCommandScript("END")
      return 0

#====================================================================================================
# PLUGIN TESTS
# See Python documentation on unittest module

class testbaverage(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     workDirectory = os.path.join(CCP4Utils.getTestTmpDir(),'test1')
     if not os.path.exists(workDirectory): os.mkdir(workDirectory)

     self.wrapper = baverage(parent=QTAPPLICATION(),name='baverage_test1',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','baverage','test_data','test1.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testbaverage)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
