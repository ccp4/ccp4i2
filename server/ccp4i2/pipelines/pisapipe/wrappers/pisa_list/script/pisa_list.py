
from ccp4i2.core.CCP4PluginScript import CPluginScript

  
class pisa_list(CPluginScript):

    TASKTITLE = 'Structure analysis with Pisa'
    TASKNAME = 'pisa_list'
    TASKCOMMAND = 'pisa'

    def makeCommandAndScript(self):
      self.appendCommandLine([self.container.controlParameters.IDENTIFIER.__str__(), '-list','monomers'] )
