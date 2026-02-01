import os

from ccp4i2.core import CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript


class import_serial(CPluginScript):
    TASKMODULE = 'data_entry'
    TASKTITLE = 'Import Serial Core'
    TASKNAME = 'import_serial'
    TASKCOMMAND = 'import_serial'
    MAINTAINER = 'martin.maly@soton.ac.uk'

    def __init__(self, *args, **kwargs):
        print("import_serial: init")
        self.hklin = None
        self.hklin1 = None
        self.hklin2 = None
        self.stream = None
        CPluginScript.__init__(self, *args, **kwargs)

    def processInputFiles(self):
        print("import_serial: processInputFiles")
        self.hklin = self.container.inputData.HKLIN.fullPath.__str__()
        if self.container.inputData.HKLIN1.isSet() and self.container.inputData.HKLIN2.isSet():
            self.hklin1 = self.container.inputData.HKLIN1.fullPath.__str__()
            self.hklin2 = self.container.inputData.HKLIN2.fullPath.__str__()
        if self.container.inputData.REFERENCEFILE.isSet() and self.container.inputData.SYMMETRY_SOURCE == "reference":
            self.reference = self.container.inputData.REFERENCEFILE.fullPath.__str__()
        if self.container.inputData.CELLFILE.isSet() and self.container.inputData.SYMMETRY_SOURCE == "cellfile":
            self.cellfile = self.container.inputData.CELLFILE.fullPath.__str__()
        if self.container.inputData.STREAMFILE.isSet() and self.container.inputData.SYMMETRY_SOURCE == "streamfile":
            self.streamfile = self.container.inputData.STREAMFILE.fullPath.__str__()
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        print("import_serial: makeCommandAndScript start")
        self.appendCommandLine("--hklin")
        self.appendCommandLine(str(self.hklin))
        if self.container.inputParameters.SPACEGROUP.isSet():
            self.appendCommandLine("--spacegroup")
            self.appendCommandLine(str(self.container.inputParameters.SPACEGROUP).strip())
        if self.container.inputParameters.CELL.isSet():
            self.appendCommandLine("--cell")
            cell_list = str(self.container.inputParameters.CELL).split()
            for param in cell_list:
                self.appendCommandLine(param)
        # self.appendCommandLine(str(self.container.inputParameters.CELL))
        if self.container.inputData.HKLIN1.isSet() and self.container.inputData.HKLIN2.isSet():
            self.appendCommandLine("--half-dataset")
            self.appendCommandLine(str(self.hklin1))
            self.appendCommandLine(str(self.hklin2))
        if self.container.inputParameters.N_BINS.isSet():
            self.appendCommandLine("--nbins")
            self.appendCommandLine(str(self.container.inputParameters.N_BINS))
        if self.container.inputParameters.D_MIN:
            self.appendCommandLine("--dmin")
            self.appendCommandLine(str(self.container.inputParameters.D_MIN))
        if self.container.inputParameters.D_MAX.isSet():
            self.appendCommandLine("--dmax")
            self.appendCommandLine(str(self.container.inputParameters.D_MAX))
        if self.container.inputParameters.WAVELENGTH.isSet():
            self.appendCommandLine("--wavelength")
            self.appendCommandLine(str(self.container.inputParameters.WAVELENGTH))
        # --reference --cellfile --streamfile are not needed because
        # symmetry should be specified using --spacegroup and --cell
        if self.container.inputData.CELLFILE.isSet() and self.container.inputData.SYMMETRY_SOURCE == "cellfile":
            self.appendCommandLine("--cellfile")
            self.appendCommandLine(str(self.cellfile))
        if self.container.inputData.REFERENCEFILE.isSet() and self.container.inputData.SYMMETRY_SOURCE == "reference":
            self.appendCommandLine("--reference")
            self.appendCommandLine(str(self.reference))
        if self.container.inputData.STREAMFILE.isSet() and self.container.inputData.SYMMETRY_SOURCE == "streamfile":
            self.appendCommandLine("--streamfile")
            self.appendCommandLine(str(self.streamfile))

        # XML output 'program.xml' is produced by the command line application
        self.xmlout = self.makeFileName('PROGRAMXML')
        print("import_serial: makeCommandAndScript end")
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        print("import_serial: processOutputFiles start")
        self.container.outputData.HKLOUT.setFullPath(os.path.join(self.getWorkDirectory(), "project_dataset.mtz"))
        self.container.outputData.HKLOUT.setAnnotation("Merged intensities")
        self.container.outputData.HKLOUT.contentFlag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN # not sure if does anything

        print("import_serial: processOutputFiles end")
        return CPluginScript.SUCCEEDED
