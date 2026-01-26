"""
     convert2mtz.py: CCP4 GUI Project
     Copyright (C) 2011 STFC
     Author: Martyn Winn

     Wrapper to convert2mtz
"""

from ccp4i2.wrappers.x2mtz.script import x2mtz

class convert2mtz(x2mtz.x2mtz):

    TASKMODULE = 'test' # Where this plugin will appear on gui
    TASKTITLE = 'Convert merged reflection file to MTZ' # A short title for gui menu
    TASKNAME = 'convert2mtz'   # Task name - should be same as class name
    TASKVERSION= 0.1               # Version of this plugin
    TASKCOMMAND = 'convert2mtz'   # The command to run the executable

    def makeCommandAndScript(self):

      self.appendCommandLine(['-stdin'])

      # input parameters
      self.appendCommandScript("hklin %s" % self.container.inputData.HKLIN)
      cell = self.container.inputData.SPACEGROUPCELL.cell
      if cell.a.isSet() and cell.b.isSet() and cell.c.isSet():
         com = "cell %f %f %f" % (cell.a,cell.b,cell.c)
         if cell.alpha.isSet() and cell.beta.isSet() and cell.gamma.isSet():
           com = com + "  %f %f %f" % (cell.alpha,cell.beta,cell.gamma)
         else:
            com = com + " 90.0 90.0 90.0"
         self.appendCommandScript(com)
      if self.container.inputData.SPACEGROUPCELL.spaceGroup.isSet():
         self.appendCommandScript("spacegroup %s" % self.container.inputData.SPACEGROUPCELL.spaceGroup)
      if self.container.inputData.XYZIN.isSet():
         self.appendCommandScript("pdbin %s" % self.container.inputData.XYZIN)

      # output parameters
      if self.container.outputData.HKLOUT.isSet():
         self.appendCommandScript("mtzout %s" % self.container.outputData.HKLOUT)

      # control parameters
      if self.container.controlParameters.NO_COMPLETE:
          self.appendCommandScript("no-complete")
      if self.container.controlParameters.ANOMALOUS:
          self.appendCommandScript("anomalous")
      if self.container.controlParameters.SEED.isSet():
          self.appendCommandScript("seed %f" % self.container.controlParameters.SEED)

      return 0


