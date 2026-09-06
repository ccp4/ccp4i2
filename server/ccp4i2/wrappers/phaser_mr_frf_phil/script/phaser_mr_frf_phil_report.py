"""Report for phaser_mr_frf_phil: the rotation peaks, no verdict."""
from ccp4i2.wrappers.phaser_phil.script.phaser_report_base import PhaserReportBase


class phaser_mr_frf_phil_report(PhaserReportBase):
    TASKNAME = "phaser_mr_frf_phil"

    def drawResults(self, parent, running):
        rotations = self.xmlnode.find("Rotations")
        if rotations is None:
            parent.addText(text="No rotation list recorded yet" if running else "No rotation list recorded")
            return
        fold = parent.addFold(label=f"Rotation peaks ({len(rotations)})", brief="Rotations", initiallyOpen=True)
        table = fold.addTable(xmlnode=rotations, select="Rotation", style="width:600px;",
                              outputXml=self.outputXml, internalId="PhaserRotationsTable")
        for title, select in (("#", "Number"), ("Ensemble", "Name"), ("Euler", "Euler"),
                              ("RF", "RF"), ("RFZ", "RFZ")):
            table.addData(title=title, select=select)
        self.drawTimeline(parent)
