
from ccp4i2.core.CCP4PluginScript import CPluginScript


class pisa_analyse(CPluginScript):

    TASKTITLE = 'Structure analysis with Pisa'     # A short title for gui menu
    TASKNAME = 'pisa_analyse'                                  # Task name - should be same as class name
    TASKCOMMAND = 'pisa'                                     # The command to run the executable

    def makeCommandAndScript(self):
      self.appendCommandLine([self.container.controlParameters.IDENTIFIER.__str__(), '-analyse',self.container.inputData.PDBIN.__str__()] )
      return 0
