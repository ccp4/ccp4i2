from ccp4i2.core.CCP4PluginScript import CPluginScript
  
class baverage(CPluginScript):

    TASKMODULE = 'wrappers'  # Where this plugin will appear on the gui
    TASKTITLE = 'Average B over main and side chain atoms'
    TASKNAME = 'baverage'
    TASKCOMMAND = 'baverage'

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
