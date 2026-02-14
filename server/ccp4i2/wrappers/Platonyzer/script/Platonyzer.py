import os

from ccp4i2.core.CCP4PluginScript import CPluginScript


class Platonyzer(CPluginScript):
    TASKNAME = 'Platonyzer'
    TASKCOMMAND="platonyzer"

    def makeCommandAndScript(self):
        if str(self.container.controlParameters.MODE) == 'NA_MG':
           self.appendCommandLine(['--create-na-mg-links'])
           if self.container.controlParameters.RM_VDW:
              self.appendCommandLine(['--delete-vdw-rest'])
        self.appendCommandLine([self.container.inputData.XYZIN])
        self.appendCommandLine([self.container.outputData.XYZOUT])

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        self.container.outputData.RESTRAINTS.setFullPath(os.path.join(self.getWorkDirectory(),"XYZOUT.restraints"))
        if not os.path.exists(str(self.container.outputData.RESTRAINTS)):
           return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED
