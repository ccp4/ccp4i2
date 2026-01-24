"""
     mtzutils.py: CCP4 GUI Project
     Copyright (C) 2011 STFC
     Author: Martyn Winn
"""

from ccp4i2.core.CCP4PluginScript import CPluginScript

class mtzutils(CPluginScript):

    TASKTITLE = 'Add or delete MTZ columns' # A short title for gui menu
    TASKNAME = 'mtzutils'   # Task name - should be same as class name
    TASKCOMMAND = 'mtzutils'            # The command to run the executable
    TASKVERSION= 0.0               # Version of this plugin

    def makeCommandAndScript(self):

      self.appendCommandLine(['HKLIN1',self.container.inputData.HKLIN1.fullPath])
      if self.container.inputData.HKLIN2.isSet():
         self.appendCommandLine(['HKLIN2',self.container.inputData.HKLIN2.fullPath])
      self.appendCommandLine(['HKLOUT',self.container.outputData.HKLOUT.fullPath])

      if self.container.controlParameters.TITLE.isSet():
          self.appendCommandScript("TITLE %s"%(str(self.container.controlParameters.TITLE)))

      if self.container.controlParameters.MODE=="INCLUDE":
         inc_exc_line = "INCLUDE "
      elif self.container.controlParameters.MODE=="EXCLUDE":
         inc_exc_line = "EXCLUDE "

      # check for column assignments
      # for the moment, we assume these are just for INCL/EXCL keywords
      if self.container.inputData.FSIGF.isSet():
         inc_exc_line = inc_exc_line + str(self.container.inputData.FSIGF.F) + ' ' + str(self.container.inputData.FSIGF.SIGF)

      self.appendCommandScript(inc_exc_line)

      self.appendCommandScript('END')

      return 0

     
