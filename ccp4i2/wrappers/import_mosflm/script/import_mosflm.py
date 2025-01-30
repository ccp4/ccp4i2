import shutil

from ....core.CCP4PluginScript import CPluginScript


class import_mosflm(CPluginScript):

    TASKNAME = 'import_mosflm'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    COMTEMPLATE = None                                   # The program com file template
    COMTEMPLATEFILE = None                               # Name of file containing com file template
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        
        self.checkOutputData()
        
        programXMLPath = self.makeFileName('PROGRAMXML')
        shutil.copyfile(self.container.inputData.MOSFLMXML.fullPath.__str__(), programXMLPath)
        self.reportStatus(CPluginScript.SUCCEEDED)
