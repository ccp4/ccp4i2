import os
from pathlib import Path

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
        if par.GPU:
            self.appendCommandLine("--gpu")
        self.appendCommandLine(["--overlap", par.OVERLAP])
        if par.FULL_CELL:
            self.appendCommandLine("--use-unit-cell")
        if par.RESOLUTION.isSet():
            self.appendCommandLine(["--resolution", par.RESOLUTION])
        if par.OUTPUT_TYPE == "RAW":
            self.appendCommandLine("--raw")
        elif par.OUTPUT_TYPE == "VARIANCE":
            self.appendCommandLine("--variance")
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        out = self.container.outputData
        directory = Path(self.getWorkDirectory())
        os.rename(directory / "nucleofind-phosphate.map", str(out.PHOSPHATE))
        os.rename(directory / "nucleofind-sugar.map", str(out.SUGAR))
        os.rename(directory / "nucleofind-base.map", str(out.BASE))
        out.PHOSPHATE.annotation.set("NucleoFind Predicted Phosphate Map")
        out.SUGAR.annotation.set("NucleoFind Predicted Sugar Map")
        out.BASE.annotation.set("NucleoFind Predicted Base Map")
        return CPluginScript.SUCCEEDED
