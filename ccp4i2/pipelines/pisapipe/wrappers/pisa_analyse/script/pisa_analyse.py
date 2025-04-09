
from ......core.CCP4PluginScript import CPluginScript


class pisa_analyse(CPluginScript):

    TASKMODULE = None                           # Where this plugin will appear on the gui
    TASKTITLE = 'Structure analysis with Pisa'     # A short title for gui menu
    TASKNAME = 'pisa_analyse'                                  # Task name - should be same as class name
    TASKCOMMAND = 'pisa'                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template

    def makeCommandAndScript(self):
      self.appendCommandLine([self.container.controlParameters.IDENTIFIER.__str__(), '-analyse',self.container.inputData.PDBIN.__str__()] )
      return 0
