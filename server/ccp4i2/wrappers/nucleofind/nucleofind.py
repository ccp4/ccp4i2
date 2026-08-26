import os
from pathlib import Path

from ccp4i2.core.CCP4PluginScript import CPluginScript


class nucleofind(CPluginScript):
    # Without TASKNAME, CPluginScript.__init__ never looks for a def.xml, so
    # the container comes up empty and i2run builds a parser with no task
    # arguments -- which is why every nucleofind run since the task was added
    # ended in "unrecognized arguments: --FPHIIN" and SystemExit(2).
    TASKNAME = "nucleofind"
    TASKCOMMAND = "nucleofind"
    WHATNEXT = ["coot_rebuild", "coot1"]

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
        par = self.container.controlParameters
        directory = Path(self.getWorkDirectory())
        # str() first: a CString compares equal to a str but does not hash
        # equal to one, so it silently misses as a dict key and every RAW or
        # VARIANCE run looked for the plain ".map" that nucleofind had not
        # written. (The identity hash is deliberate --- see
        # tests/unit/containers/test_hash_collision_fix.py --- so the caller
        # converts.)
        ext_map = {"RAW": ".raw.map", "VARIANCE": ".variance.map"}
        ext = ext_map.get(str(par.OUTPUT_TYPE), ".map")
        os.rename(directory / f"nucleofind-phosphate{ext}", str(out.PHOSPHATE))
        os.rename(directory / f"nucleofind-sugar{ext}", str(out.SUGAR))
        os.rename(directory / f"nucleofind-base{ext}", str(out.BASE))
        out.PHOSPHATE.annotation.set("NucleoFind Predicted Phosphate Map")
        out.SUGAR.annotation.set("NucleoFind Predicted Sugar Map")
        out.BASE.annotation.set("NucleoFind Predicted Base Map")
        return CPluginScript.SUCCEEDED
