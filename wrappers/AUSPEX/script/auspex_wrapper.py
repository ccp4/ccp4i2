from core.CCP4PluginScript import CPluginScript
from core import CCP4XtalData
from core import CCP4ErrorHandling
import os
import glob

class AUSPEX(CPluginScript):

    TASKMODULE = 'data_reduction'      # Gui menu location
    TASKTITLE = 'AUSPEX'             # Short title for Gui
    TASKNAME = 'AUSPEX'              # Task name - same as class name
    TASKCOMMAND = 'auspex'             # The command to run the executable
    TASKVERSION = 1.0                  # plugin version
    COMTEMPLATE = None                 # The program com file template
    COMTEMPLATEFILE = None             # Name of file containing com file template
    PERFORMANCECLASS = 'CExpPhasPerformance'
    ASYNCHRONOUS = True
    MAINTAINER = 'Andrea.Thorn@web.de'

    def __init__(self, *args, **kwargs):
        self.hklin = None
        CPluginScript.__init__(self, *args, **kwargs)

    def process(self):
        CPluginScript.process(self)

    def processInputFiles(self):
        cols1 = []
        self.container.inputData.F_SIGF.loadFile()
        self.container.inputData.F_SIGF.setContentFlag()
        bFData = self.container.inputData.F_SIGF.contentFlag == 2 or self.container.inputData.F_SIGF.contentFlag == 4
        bIData = self.container.inputData.F_SIGF.contentFlag == 1 or self.container.inputData.F_SIGF.contentFlag == 3
        if bIData:
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN])
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN])
        if bFData:
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN])
        self.hklin1, __, error1 = self.makeHklInput(cols1, extendOutputColnames=True, useInputColnames=True)
        if error1.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED

    def processOutputFiles(self):
        outContA = self.container.outputData.IM_OUT
        outFiles = glob.glob(os.path.normpath(os.path.join(self.getWorkDirectory(), '*.png')))
        for outFile in outFiles:
            outContA.append(outContA.makeItem())
            outContA[-1].fullPath = outFile
            outContA[-1].annotation = 'AUSPEX IMG'
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, container=None):
        self.appendCommandLine("--no-filename-in-title")
        self.appendCommandLine("%s"%str(self.hklin1))
        self.appendCommandLine("--ylim %s"%(self.container.inputData.YLIM))
        if self.container.inputData.DLIM.isSet():
            self.appendCommandLine("--dmin %f"%(float(self.container.inputData.DLIM)))
        if self.container.inputData.SINGFIG:
            self.appendCommandLine("--single-figure")
        if not self.container.inputData.FLAGICE:
            self.appendCommandLine("--no-automatic")
