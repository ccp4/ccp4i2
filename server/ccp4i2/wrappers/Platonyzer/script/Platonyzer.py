import os

from ccp4i2.core.CCP4PluginScript import CPluginScript


class Platonyzer(CPluginScript):
    TASKNAME = 'Platonyzer'   # Task name - should be same as class name and match pluginTitle in the .def.xml file
    TASKVERSION= 0.1               # Version of this plugin
    MAINTAINER = 'nicholls@mrc-lmb.cam.ac.uk'
    ERROR_CODES = { 201 : {'description' : 'Failed to analyse output files' },
                    202 : {'description' : 'Failed applying selection ot PDB file' }
                    }
    PURGESEARCHLIST = [ [ 'hklin.mtz' , 0 ],
                       ['log_mtzjoin.txt', 0]
                       ]
    TASKCOMMAND="platonyzer"
    
    def __init__(self, *args, **kws):
        super(Platonyzer, self).__init__(*args, **kws)

    def processInputFiles(self):
        return CPluginScript.SUCCEEDED

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
