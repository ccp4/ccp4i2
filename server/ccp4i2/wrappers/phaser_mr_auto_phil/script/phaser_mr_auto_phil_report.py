"""Report for phaser_mr_auto_phil: the recorded run, never an inference."""
from ccp4i2.report import Report


class phaser_mr_auto_phil_report(Report):
    TASKNAME = "phaser_mr_auto_phil"
    RUNNING = True
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.outputXml = jobStatus is not None and jobStatus.lower() == "running"
        if jobStatus is None or jobStatus.lower() == "nooutput":
            return
        self.drawContent(jobStatus, self)

    def drawContent(self, jobStatus=None, parent=None):
        if parent is None:
            parent = self
        running = jobStatus is not None and jobStatus.lower() == "running"
        if running:
            self.outputXml = False
            self.drawProgress(parent)
        self.drawVerdict(parent, running)
        self.drawWarnings(parent)
        if self.xmlnode.find("Solutions") is not None:
            fold = parent.addFold(label="Solutions", brief="Solutions", initiallyOpen=True)
            self.drawSolutions(fold)
        fold = parent.addFold(label="Search strategy", brief="Strategy", initiallyOpen=not running)
        self.drawStrategy(fold, running)
        if self.xmlnode.find("Summaries") is not None:
            fold = parent.addFold(label="Phaser's summary, module by module", brief="Summary",
                                  initiallyOpen=False)
            for node in self.xmlnode.findall("Summaries/Summary"):
                sub = fold.addFold(label=node.get("module", "?"), initiallyOpen=False)
                sub.addPre(text=node.text or "")
        self.drawGraphs(parent)

    def drawProgress(self, parent):
        module = self.xmlnode.find("CurrentModule")
        if module is not None and module.text:
            parent.addText(text=f"Now: {module.text}", style="font-weight:bold;")
        label, maximum, value = ("No action in progress", 0, 0)
        node = self.xmlnode.find("CurrentActivity")
        if node is not None:
            label = node.findtext("label") or label
            maximum = node.findtext("max") or 0
            value = node.findtext("value") or 0
        parent.addProgress(style="width:500px;border:1px solid green;", value=value,
                           max=maximum, internalId="PhaserProgress", outputXml=self.outputXml,
                           label=label)

    def drawVerdict(self, parent, running):
        verdict = self.xmlnode.findtext("Verdict")
        if verdict:
            colour = "green" if "Single" in verdict else "orange"
            parent.addText(text=f"Phaser: {verdict}", style=f"color:{colour};font-weight:bold;")
        elif not running:
            parent.addText(text="Phaser reported no verdict", style="color:red;")
        best = self.xmlnode.find("PhaserCurrentBestSolution")
        if running and best is not None and best.text:
            parent.addText(text="Current best solution: " + best.text.strip())

    def drawWarnings(self, parent):
        warnings = [w.text for w in self.xmlnode.findall("PhaserWarnings/Warning") if w.text]
        if warnings:
            fold = parent.addFold(label=f"Phaser warnings ({len(warnings)})", initiallyOpen=False)
            for text in warnings:
                fold.addText(text=text, style="color:orange;")
                fold.addDiv(style="clear:both;")

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
        modules = self.xmlnode.findall("Modules/Module")
        if modules:
            fold = parent.addFold(label=f"Module timeline ({len(modules)})", initiallyOpen=False)
            fold.addPre(text="\n".join(f"{m.get('t', ''):>7}s  {m.get('name', '')}" for m in modules))

    def drawGraphs(self, parent):
        tables = self.xmlnode.findall(".//CCP4ApplicationOutput/CCP4Table")
        if not tables:
            return
        gallery = parent.addObjectGallery(style="float:left;", height="300px", width="700px",
                                          tableWidth="260px", contentWidth="420px")
        for i, table in enumerate(tables):
            graph = gallery.addFlotGraph(xmlnode=table, title=table.get("title"),
                                         internalId=f"PhaserGraph{i}", outputXml=self.outputXml,
                                         label=table.get("title"), style="width:420px;")
            graph.addPimpleData(xmlnode=table, usePlotly=False)
        parent.addDiv(style="clear:both;")
