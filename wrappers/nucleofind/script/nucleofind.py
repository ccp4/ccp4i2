import os

from core.CCP4PluginScript import CPluginScript


class nucleofind(CPluginScript):
    TASKMODULE = 'model_building'
    TASKCOMMAND = 'nucleofind'
    WHATNEXT = ['coot_rebuild']
    MAINTAINER = 'jordan.dialpuri@york.ac.uk'

    def makeCommandAndScript(self):
        inp = self.container.inputData
        par = self.container.controlParameters
        self.appendCommandLine(["--input", inp.FPHIIN])
        self.appendCommandLine(["--amplitude", "F"])
        self.appendCommandLine(["--phase", "PHI"])
        self.appendCommandLine(["--output", "."])
        self.appendCommandLine(["--model", "core"])
        self.appendCommandLine(["--nthreads", par.THREADS])
        if par.GPU.isSet():
            self.appendCommandLine("--gpu")
        self.appendCommandLine(["--overlap", par.OVERLAP])
        if par.FULL_CELL.isSet():
            self.appendCommandLine("--use-symmetry")
        if par.RESOLUTION.isSet():
            self.appendCommandLine(["--resolution", par.RESOLUTION])
        if par.OUTPUT_TYPE == "RAW":
            self.appendCommandLine("--raw")
        elif par.OUTPUT_TYPE == "VARIANCE":
            self.appendCommandLine("--variance")
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        directory = self.getWorkDirectory()
        phosphate = os.path.join(directory, "nucleofind-phosphate.map")
        sugar = os.path.join(directory, "nucleofind-sugar.map")
        base = os.path.join(directory, "nucleofind-base.map")
        os.rename(phosphate, str(self.container.outputData.PHOSPHATE))
        os.rename(sugar, str(self.container.outputData.SUGAR))
        os.rename(base, str(self.container.outputData.BASE))
