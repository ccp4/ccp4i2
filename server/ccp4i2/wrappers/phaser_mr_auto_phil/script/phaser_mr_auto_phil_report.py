"""Report for phaser_mr_auto_phil: the recorded run, never an inference."""
from ccp4i2.wrappers.phaser_phil.script.phaser_report_base import PhaserReportBase


class phaser_mr_auto_phil_report(PhaserReportBase):
    TASKNAME = "phaser_mr_auto_phil"

    def drawVerdict(self, parent, running):
        super().drawVerdict(parent, running)
        if not running and not self.xmlnode.findtext("Verdict"):
            parent.addText(text="Phaser reported no verdict", style="color:red;")
        best = self.xmlnode.find("PhaserCurrentBestSolution")
        if running and best is not None and best.text:
            parent.addText(text="Current best solution: " + best.text.strip())

    def drawResults(self, parent, running):
        if self.xmlnode.find("Solutions") is not None:
            fold = parent.addFold(label="Solutions", brief="Solutions", initiallyOpen=True)
            self.drawSolutions(fold)
        fold = parent.addFold(label="Search strategy", brief="Strategy", initiallyOpen=not running)
        self.drawStrategy(fold, running)

    def drawSolutions(self, parent):
        solutions = self.xmlnode.find("Solutions")
        table = parent.addTable(xmlnode=solutions, select="Solution", style="width:700px;",
                                outputXml=self.outputXml, internalId="PhaserSolutionsTable")
        for title, select in (("#", "Number"), ("Space group", "spaceGroup"), ("LLG", "LLG"),
                              ("TFZ", "TFZ"), ("TFZ-equiv", "TFZeq"), ("R", "R"),
                              ("Clashes", "PAK"), ("Annotation", "Annotation")):
            table.addData(title=title, select=select)
        for i, solution in enumerate(solutions.findall("Solution")):
            placements = solution.find("Placements")
            if placements is None or len(placements) == 0:
                continue
            fold = parent.addFold(label=f"Solution {i + 1}: placements", initiallyOpen=(i == 0))
            table = fold.addTable(xmlnode=placements, select="Placement", style="width:600px;",
                                  outputXml=self.outputXml, internalId=f"PhaserPlacements{i}")
            for title, select in (("RFZ", "RFZ"), ("TFZ", "TFZ"), ("TFZ-equiv", "TFZeq"),
                                  ("Clashes", "PAK"), ("LLG", "LLG"), ("Note", "Note")):
                table.addData(title=title, select=select)
            history = solution.findtext("History")
            if history:
                fold.addText(text=f"History: {history}")
            unknown = solution.findtext("UnknownTokens")
            if unknown:
                fold.addText(text=f"Annotation tokens not understood: {unknown}",
                             style="color:orange;")

    def drawStrategy(self, parent, running):
        strategy = self.xmlnode.find("Strategy")
        attempts = strategy.findall("Attempt") if strategy is not None else \
            self.xmlnode.findall("Attempts/Attempt")
        if not attempts:
            parent.addText(text="No search attempt recorded yet" if running
                           else "No search attempt recorded")
        for attempt in attempts:
            component = attempt.get("component", "?")
            resolution = attempt.get("resolution")
            where = "" if resolution is None else (
                " at full resolution" if resolution == "full" else f" at {resolution} A")
            llg = attempt.get("llg")
            outcome = attempt.get("outcome", "?") + (f", LLG {llg}" if llg else "")
            parent.addText(text=f"Component {component}{where}: {outcome}")
            for warning in attempt.findall("Warning"):
                parent.addText(text="    " + (warning.text or ""), style="color:orange;")
            parent.addDiv(style="clear:both;")
        if strategy is not None and int(strategy.get("unparsed", "0")) > 0:
            parent.addText(text=f"{strategy.get('unparsed')} strategy block(s) in Phaser's "
                                "summary matched no known sentence -- see the summary below",
                           style="color:orange;")
        self.drawTimeline(parent)
