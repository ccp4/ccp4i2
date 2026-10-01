from ccp4i2.wrappers.phaser_phil.script.phaser_report_base import PhaserReportBase
from ccp4i2.wrappers.phaser_mr_auto_phil.script.phaser_mr_auto_phil_report import phaser_mr_auto_phil_report


class phaser_mr_rnp_phil_report(phaser_mr_auto_phil_report):
    TASKNAME = "phaser_mr_rnp_phil"
    # These modes (and FTF and PAK, which inherit this report) work on what
    # they are given and record no search of their own: "No search attempt
    # recorded" would read as something missing.
    SEARCHES = False

    def drawVerdict(self, parent, running):
        # Refinement of placed solutions: Phaser gives the refined LLGs, no
        # verdict on a search, so none is missing when there is none
        PhaserReportBase.drawVerdict(self, parent, running)
