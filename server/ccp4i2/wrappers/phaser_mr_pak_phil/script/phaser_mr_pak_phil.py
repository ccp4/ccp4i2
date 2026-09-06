"""The packing test with Phaser (MR_PAK), over its PHIL: the MR_AUTO task
with the mode changed and a set of solutions to test required. The
ensembles define the components the solutions name; none need ask for
copies.
"""
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil import phaser_mr_auto_phil


class phaser_mr_pak_phil(phaser_mr_auto_phil):

    TASKNAME = "phaser_mr_pak_phil"
    PHIL_MODE = "MR_PAK"
    SEARCHES_ENSEMBLES = False   # tests what SOLIN places
    WHATNEXT = ["phaser_mr_rnp_phil"]
    ERROR_CODES = {**phaser_mr_auto_phil.ERROR_CODES,
                   116: {"description": "Nothing to refine or test"}}

    def validity(self):
        error = super().validity()
        if not self.container.inputData.SOLIN.isSet():
            error.append(klass=self.TASKNAME, code=116,
                         details="The packing test needs the solutions of a translation-function job.",
                         name=f"{self.TASKNAME}.container.inputData.SOLIN",
                         severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error

    def processOutputFiles(self):
        self.writeSolutions()
        self.recordRun(with_strategy=False)
        return CPluginScript.SUCCEEDED
