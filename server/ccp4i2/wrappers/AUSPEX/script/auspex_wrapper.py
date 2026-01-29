from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4XtalData
from ccp4i2.core import CCP4ErrorHandling
import os
import glob

class AUSPEX(CPluginScript):

    TASKMODULE = 'data_reduction'
    TASKTITLE = 'AUSPEX'
    TASKNAME = 'AUSPEX'
    TASKCOMMAND = 'auspex'
    PERFORMANCECLASS = 'CExpPhasPerformance'
    ASYNCHRONOUS = True
    MAINTAINER = 'Andrea.Thorn@web.de'

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hklin = None

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
        self.hklin, __, error1 = self.makeHklInput(cols1, extendOutputColnames=True, useInputColnames=True)
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

    def makeCommandAndScript(self):
        self.appendCommandLine("--no-filename-in-title")
        self.appendCommandLine(self.hklin)
        self.appendCommandLine(["--ylim", self.container.inputData.YLIM])
        if self.container.inputData.DLIM.isSet():
            self.appendCommandLine(["--dmin", self.container.inputData.DLIM])
        if self.container.inputData.SINGFIG:
            self.appendCommandLine("--single-figure")
        if not self.container.inputData.FLAGICE:
            self.appendCommandLine("--no-automatic")
