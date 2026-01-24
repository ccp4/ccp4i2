
from ccp4i2.core.CCP4PluginScript import CPluginScript

  
class import_mosflm(CPluginScript):

    TASKNAME = 'import_mosflm'                                  # Task name - should be same as class name
    TASKCOMMAND = ''                                     # The command to run the executable
    TASKVERSION= 0.0                                     # Version of this plugin
    MAINTAINER = 'martin.noble@newcastle.ac.uk'

    '''
    def __init__(self,parent=None,name=None,workDirectory=''):
      CPluginScript. __init__(self,parent=parent,name=name)
    '''
    
    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        
        self.checkOutputData()
        
        import shutil
        programXMLPath = self.makeFileName('PROGRAMXML')
        shutil.copyfile(self.container.inputData.MOSFLMXML.fullPath.__str__(), programXMLPath)
        self.reportStatus(CPluginScript.SUCCEEDED)
