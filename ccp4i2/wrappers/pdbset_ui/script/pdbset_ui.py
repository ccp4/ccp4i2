from ....core.CCP4PluginScript import CPluginScript


class pdbset_ui(CPluginScript):

    TASKTITLE='Scripted structure edits - Pdbset'
    TASKNAME = 'pdbset_ui'
    TASKMODULE= 'model_data_utility'
    TASKCOMMAND = 'pdbset'
    TASKVERSION= 0.0
    COMLINETEMPLATE = None
    COMTEMPLATE = None

    def makeCommandAndScript(self):
      inp = self.container.inputData
      out = self.container.outputData

      if inp.XYZIN.fullPath.isSet():
          xyzin_target_file = str( inp.XYZIN.fullPath )
          self.appendCommandLine( [ "XYZIN",xyzin_target_file ] )

      if out.XYZOUT.fullPath.isSet():
          xyzout_target_file = str( out.XYZOUT.fullPath )
          self.appendCommandLine( [ "XYZOUT",xyzout_target_file ] )

      s = ""
      if self.container.controlParameters.EXTRA_PDBSET_KEYWORDS.isSet():
           for kwLine in str(self.container.controlParameters.EXTRA_PDBSET_KEYWORDS).split('\n'):
              kw = kwLine.lstrip().rstrip()
              if len(kw)>0:
                 if str(kw)[0] != '#':
                    if kw == "END":
                        break
                    s += kw + "\n"

      print("Keyword script:",kw)

      self.appendCommandScript( s+"\nEND\n" )

      return CPluginScript.SUCCEEDED
