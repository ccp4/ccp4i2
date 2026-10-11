
from ccp4i2.report import Report


class dials_rlattice_report(Report):
    TASKNAME = 'dials_rlattice'

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.addText(text="Running the Dials Reciprocal Lattice Viewer")
