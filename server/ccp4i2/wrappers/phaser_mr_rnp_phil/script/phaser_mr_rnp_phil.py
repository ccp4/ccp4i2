"""Rigid-body refinement of placed solutions with Phaser (MR_RNP), over its
PHIL: the MR_AUTO task with the mode changed and a set of solutions to
start from required. The ensembles define the components the solutions
name; none need ask for copies to be searched for.
"""
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil import phaser_mr_auto_phil


class phaser_mr_rnp_phil(phaser_mr_auto_phil):

    TASKNAME = "phaser_mr_rnp_phil"
    PHIL_MODE = "MR_RNP"
    SEARCHES_ENSEMBLES = False   # refines what SOLIN places

    def validity(self):
        error = super().validity()
        if not self.container.inputData.SOLIN.isSet():
            error.append(klass=self.TASKNAME, code=116,
                         details="Rigid-body refinement needs the solutions from a previous "
                                 "Phaser job to refine.",
                         name=f"{self.TASKNAME}.container.inputData.SOLIN",
                         severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error
