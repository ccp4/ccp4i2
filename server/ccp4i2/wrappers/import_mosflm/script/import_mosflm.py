
from ccp4i2.core.CCP4PluginScript import CPluginScript

  
class import_mosflm(CPluginScript):

    TASKNAME = 'import_mosflm'

    def process(self):
        invalidFiles = self.checkInputData()
        if len(invalidFiles)>0:
            self.reportStatus(CPluginScript.FAILED)
        
        self.checkOutputData()
        
        import shutil
        programXMLPath = self.makeFileName('PROGRAMXML')
        shutil.copyfile(self.container.inputData.MOSFLMXML.fullPath.__str__(), programXMLPath)
        self.reportStatus(CPluginScript.SUCCEEDED)
