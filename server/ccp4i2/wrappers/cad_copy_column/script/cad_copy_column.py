from __future__ import print_function

"""
     cad_copy_column.py: CCP4 GUI Project
     Copyright (C) 2011 STFC
"""

from ccp4i2.core.CCP4PluginScript import CPluginScript

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
