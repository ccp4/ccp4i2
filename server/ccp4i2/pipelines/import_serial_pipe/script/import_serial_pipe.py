import shutil

from ccp4i2.core.CCP4PluginScript import CPluginScript


class import_serial_pipe(CPluginScript):
    TASKTITLE = 'Import Serial Pipeline'
    TASKNAME = 'import_serial_pipe'
    MAINTAINER = 'martin.maly@soton.ac.uk'

    def __init__(self, *args, **kwargs):
        self.hklin = None
        self.hklin1 = None
        self.hklin2 = None
        self.stream = None
        CPluginScript.__init__(self, *args, **kwargs)

    def process(self):
        self.importSerialProcess = self.makePluginObject('import_serial')
        self.importSerialProcess.container.inputParameters = self.container.inputParameters
        self.importSerialProcess.container.inputData.copyData(otherContainer=self.container.inputData)
        self.importSerialProcess.container.outputData.copyData(otherContainer=self.container.outputData)
        self.xmlout = self.makeFileName('PROGRAMXML')
        print("import_serial_pipe: import_serial starting")
        status = self.importSerialProcess.process()
        print("import_serial_pipe: import_serial ended")

        # Run aimless for a report on data quality
        self.aimlessPipe = self.makePluginObject('aimless_pipe', pluginTitle='DR run for data analysis')
        mergedList = self.aimlessPipe.container.inputData.UNMERGEDFILES
        if len(mergedList)==0: mergedList.addItem()
        # Always do analysis on the file which is saved as pipeline output
        mergedList[0].file.set(self.importSerialProcess.container.outputData.HKLOUT.__str__())
        # parameters for Pointless
        self.aimlessPipe.container.controlParameters.MODE = 'CHOOSE'
        self.aimlessPipe.container.controlParameters.CHOOSE_MODE = 'SPACEGROUP'
        self.aimlessPipe.container.controlParameters.CHOOSE_SPACEGROUP = \
            str(self.container.inputParameters.SPACEGROUP).strip()
        # parameters for Aimless
        self.aimlessPipe.container.controlParameters.SCALING_PROTOCOL = 'CONSTANT'
        self.aimlessPipe.container.controlParameters.ONLYMERGE = True
        self.aimlessPipe.container.controlParameters.ANALYSIS_MODE = True
        self.aimlessPipe.container.controlParameters.OUTPUT_UNMERGED = False
        self.aimlessPipe.container.controlParameters.SDCORRECTION_OVERRIDE = True
        self.aimlessPipe.container.controlParameters.SDCORRECTION_REFINE = False
        self.aimlessPipe.container.controlParameters.SDCORRECTION_SET = True
        self.aimlessPipe.container.controlParameters.SDCORRECTION_SDFAC = 1.0
        self.aimlessPipe.container.controlParameters.SDCORRECTION_SDB = 0.0
        self.aimlessPipe.container.controlParameters.SDCORRECTION_SDADD = 0.0

        #  Start data reduction
        print("import_serial_pipe: import_merged aimless_pipe starting")
        self.aimlessPipe.process()
        print("import_serial_pipe: import_merged aimless_pipe end")
        shutil.copyfile(str(self.aimlessPipe.container.outputData.IMEANOUT[0]), str(self.container.outputData.HKLOUT))
        self.container.outputData.HKLOUT.set(self.aimlessPipe.container.outputData.HKLOUT)
        self.container.outputData.HKLOUT.setAnnotation("Merged intensities")
        self.container.outputData.HKLOUT.setContentFlag(reset=True) # this will set contentFlag to 3 (CONTENT_FLAG_IMEAN)

        # Save xml
        self.xmlout = self.makeFileName('PROGRAMXML')
        importSerialProcessXMLpath = self.importSerialProcess.makeFileName('PROGRAMXML')
        importSerialPipelineXMLpath = self.makeFileName('PROGRAMXML')
        shutil.copyfile(str(importSerialProcessXMLpath), str(importSerialPipelineXMLpath))
        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED
