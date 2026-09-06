"""The translation function with Phaser (MR_FTF), over its PHIL: the MR_AUTO
task with the mode changed and a rotation list to try in place of a set
of solutions. The ensembles define the components the list names; none
need ask for copies.
"""
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil import phaser_mr_auto_phil
from ccp4i2.wrappers.phaser_phil.script.phaser_shims import is_rotation_list


class phaser_mr_ftf_phil(phaser_mr_auto_phil):

    TASKNAME = "phaser_mr_ftf_phil"
    PHIL_MODE = "MR_FTF"
    SEARCHES_ENSEMBLES = False   # the rotation list says what to place
    SOLUTION_INPUT = "RFILEIN"
    WHATNEXT = ["phaser_mr_pak_phil", "phaser_mr_rnp_phil"]

    def solutions_kind_check(self, solutions):
        if not is_rotation_list(solutions):
            return ("This is a set of solutions (rotation and translation both); the translation "
                    "function needs the rotation list of a rotation-function job.")
        return None

    def processOutputFiles(self):
        self.writeSolutions()
        self.recordRun(with_strategy=False)
        return CPluginScript.SUCCEEDED
