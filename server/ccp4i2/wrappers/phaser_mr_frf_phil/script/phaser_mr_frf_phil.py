"""The rotation function with Phaser (MR_FRF), over its PHIL: the MR_AUTO
task with the mode changed, and a rotation list out instead of placed
models. The translation-function task takes the list up.
"""
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil import phaser_mr_auto_phil
from ccp4i2.wrappers.phaser_phil.script import phaser_run


class phaser_mr_frf_phil(phaser_mr_auto_phil):

    TASKNAME = "phaser_mr_frf_phil"
    PHIL_MODE = "MR_FRF"
    WHATNEXT = ["phaser_mr_ftf_phil"]

    def processOutputFiles(self):
        import pickle
        result = self.resultObject
        rotations = result.getDotRlist()
        if len(rotations) > 0:
            out = self.container.outputData
            with open(str(out.RFILEOUT.fullPath), "wb") as handle:
                pickle.dump(rotations, handle)
            out.RFILEOUT.annotation.set("Rotation list from Phaser")
        self.recordRun(with_strategy=False)
        return CPluginScript.SUCCEEDED

    def recordSolutions(self, result):
        phaser_run.rotations_xml(result, self.xmlroot)
