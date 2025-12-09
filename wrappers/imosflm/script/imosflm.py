
from core.CCP4PluginScript import CPluginScript
from ccp4i2.baselayer import QtCore
import os,re,time,sys
import platform

class imosflm(CPluginScript):
    TASKMODULE = 'data_processing'                        # Where this plugin will appear on the gui
    TASKTITLE = 'Integrate images - iMosflm'     # A short title for gui menu
    DESCRIPTION = 'Launch iMosflm and capture output'
    TASKNAME = 'imosflm'                                  # Task name - should be same as class name
    TASKCOMMAND = 'imosflm'                                     # The command to run the executable
    if platform.system() == 'Windows':
        TASKCOMMAND = 'imosflm.bat'
    TASKVERSION= 0.0                                     # Version of this plugin
    ASYNCHRONOUS = False
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = 'martin.noble@newcastle.ac.uk'
    
    ERROR_CODES = {  200 : { 'description' : 'imosflm exited with error status' }, 201 : { 'description' : 'Failed in harvest operation' },202 : { 'description' : 'Failed in processOutputFiles' }}

    def makeCommandAndScript(self):
        self.appendCommandLine('--ccp4i2')

        
    def processOutputFiles(self):
        import glob, shutil
        fileList1 = glob.glob(os.path.join(str(self.getWorkDirectory()), "*.mos"))
        fileList2 = glob.glob(os.path.join(str(self.getWorkDirectory()), "*.sav"))
        outputXMLs = self.container.outputData.MOSFLMXMLOUT
        for file in fileList1+fileList2:
            shutil.copyfile(file, self.makeFileName('PROGRAMXML'))
            outputXMLs.append(outputXMLs.makeItem())
            outputXMLs[-1].setFullPath(file)
            outputXMLs[-1].annotation = os.path.basename(file)
        fileList = glob.glob(os.path.join(self.getWorkDirectory(), "*.mtz"))
        if len(fileList)==0: return CPluginScript.UNSATISFACTORY
        outputMERGED = self.container.outputData.MERGEDMTZ
        outputUNMERGED = self.container.outputData.UNMERGEDMTZ
        for file in fileList:
            basename = os.path.basename(file)
            #print "file basename", basename
            # We ought to be processing the session_files.xml manifest, but
            # for now just save ctruncate_*-unique.mtz (merged) and initial
            # mosflm output (unmerged)
            # THe ctruncate output should be processed by cmtzsplit 
            # Not sure what happens if there are multiple sweeps integrated
            # Multiple QuickScale runs seem to overwrite
            if basename.startswith('aimless_') or basename.startswith('pointless_'):
                # don't record these ones
                pass
            elif basename.startswith('ctruncate_'):
                if basename.endswith('unique.mtz'):
                    outputMERGED.append(outputMERGED.makeItem())
                    outputMERGED[-1].setFullPath(file)
                    outputMERGED[-1].annotation = os.path.basename(file)
            else:
                outputUNMERGED.append(outputUNMERGED.makeItem())
                outputUNMERGED[-1].setFullPath(file)
                outputUNMERGED[-1].annotation = os.path.basename(file)
        return CPluginScript.SUCCEEDED
            
