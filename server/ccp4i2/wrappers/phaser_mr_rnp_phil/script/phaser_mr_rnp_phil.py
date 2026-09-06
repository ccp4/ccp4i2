"""Rigid-body refinement of placed solutions with Phaser (MR_RNP), over its
PHIL: the MR_AUTO task with the mode changed and something to refine
required. That is either the solutions of a previous Phaser job (SOLIN) or
ensembles already placed at the origin of their own coordinates
(FIXENSEMBLES, Phaser's solution_at_origin) -- a model cut into rigid
bodies, say. The ensembles define the components; none need ask for copies
to be searched for.
"""
from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil import phaser_mr_auto_phil


class phaser_mr_rnp_phil(phaser_mr_auto_phil):

    TASKNAME = "phaser_mr_rnp_phil"
    PHIL_MODE = "MR_RNP"
    SEARCHES_ENSEMBLES = False   # refines what SOLIN places

    def validity(self):
        error = super().validity()
        inp = self.container.inputData
        if not inp.SOLIN.isSet() and len(inp.FIXENSEMBLES) == 0:
            error.append(klass=self.TASKNAME, code=116,
                         details="Rigid-body refinement needs something placed to refine: the "
                                 "solutions of a previous Phaser job, or search models already "
                                 "placed at the origin of their coordinates.",
                         name=f"{self.TASKNAME}.container.inputData.SOLIN",
                         severity=CCP4ErrorHandling.SEVERITY_ERROR)
        return error
