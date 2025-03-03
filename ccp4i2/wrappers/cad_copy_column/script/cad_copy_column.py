from __future__ import print_function

"""
     cad_copy_column.py: CCP4 GUI Project
     Copyright (C) 2011 STFC
"""

from core.CCP4PluginScript import CPluginScript

class cad_copy_column(CPluginScript):

    TASKMODULE = None      # Where this plugin will appear on the gui
    TASKTITLE = 'Copy an MTZ column between files' # A short title for gui menu
    TASKNAME = 'cad_copy_column'   # Task name - should be same as class name
    TASKCOMMAND = 'cad'            # The command to run the executable
    TASKVERSION= 0.0               # Version of this plugin
    COMTEMPLATE = None             # The program com file template
    COMTEMPLATEFILE = None         # Name of file containing com file template

    def makeCommandAndScript(self):

      self.appendCommandLine(['HKLIN1',self.container.inputData.HKLIN1.fullPath])
      self.appendCommandLine(['HKLIN2',self.container.inputData.HKLIN2.fullPath])
      self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT.fullPath])

      if self.container.controlParameters.TITLE.isSet():
          self.appendCommandScript("TITLE %s"%(str(self.container.controlParameters.TITLE)))

      # basic functionality assumes we are adding columns to file 1
      self.appendCommandScript('LABI FILE 1 all')

      # this caters for 2 routes of uniqueify - needs generalising
      if self.container.inputData.FREERFLAG.isSet():
         column_choices = 'E1 = %s'%(str(self.container.inputData.FREERFLAG.FREE))
      elif self.container.inputData.FSIGF.isSet():
         column_choices = 'E1 = %s E2 = %s'%(str(self.container.inputData.FSIGF.F),str(self.container.inputData.FSIGF.SIGF))
      self.appendCommandScript('LABI FILE 2 %s'%(column_choices))
      self.appendCommandScript('END')

      return 0

     
#======================================================
# PLUGIN TESTS
# See Python documentation on unittest module

import unittest

class testcad_copy_column(unittest.TestCase):

   def setUp(self):
    from core import CCP4Modules
    self.app = CCP4Modules.QTAPPLICATION()
    # make all background jobs wait for completion
    # this is essential for unittest to work
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    from core import CCP4Modules
    CCP4Modules.PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     import os
     from core import CCP4Modules, CCP4Utils

     workDirectory = CCP4Utils.getTestTmpDir()
     logFile = os.path.join(workDirectory,'cad_copy_column_test1.log')
     # Delete any existing log file
     if os.path.exists(logFile): os.remove(logFile)

     self.wrapper = cad_copy_column(parent=CCP4Modules.QTAPPLICATION(),name='cad_copy_column_test1',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','cad_copy_column','test_data','test1.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())
     #self.assertTrue(os.path.exists(logFile),'No log file found')
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testcad_copy_column)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
