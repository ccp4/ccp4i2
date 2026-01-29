
from ccp4i2.core.CCP4PluginScript import CPluginScript


class pisa_analyse(CPluginScript):

    TASKTITLE = 'Structure analysis with Pisa'
    TASKNAME = 'pisa_analyse'
    TASKCOMMAND = 'pisa'

    def makeCommandAndScript(self):
      self.appendCommandLine([self.container.controlParameters.IDENTIFIER.__str__(), '-analyse',self.container.inputData.PDBIN.__str__()] )
      return 0
