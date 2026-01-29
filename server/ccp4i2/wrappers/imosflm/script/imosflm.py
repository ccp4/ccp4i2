import glob
import os
import shutil

from ccp4i2.core.CCP4PluginScript import CPluginScript


class imosflm(CPluginScript):
    TASKMODULE = 'data_processing'
    TASKTITLE = 'Integrate images - iMosflm'
    DESCRIPTION = 'Launch iMosflm and capture output'
    TASKNAME = 'imosflm'
    TASKCOMMAND = 'imosflm'
    MAINTAINER = 'martin.noble@newcastle.ac.uk'
    
    ERROR_CODES = {  200 : { 'description' : 'imosflm exited with error status' }, 201 : { 'description' : 'Failed in harvest operation' },202 : { 'description' : 'Failed in processOutputFiles' }}

    def makeCommandAndScript(self):
        self.appendCommandLine('--ccp4i2')

        
    def processOutputFiles(self):
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
