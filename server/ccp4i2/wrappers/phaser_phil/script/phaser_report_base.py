"""The report sections every Phaser mode task shares: what the recorder
wrote -- progress, warnings, the module timeline, Phaser's summary block by
block, its loggraphs. The mode reports add their own results."""
from ccp4i2.report import Report


class PhaserReportBase(Report):
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
        self.drawResults(parent, running)
        self.drawSummaries(parent)
        self.drawGraphs(parent)

    def drawResults(self, parent, running):
        """The mode's own results; overridden."""

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

    def drawWarnings(self, parent):
        warnings = [w.text for w in self.xmlnode.findall("PhaserWarnings/Warning") if w.text]
        if warnings:
            fold = parent.addFold(label=f"Phaser warnings ({len(warnings)})", initiallyOpen=False)
            for text in warnings:
                fold.addText(text=text, style="color:orange;")
                fold.addDiv(style="clear:both;")

    def drawTimeline(self, parent):
        modules = self.xmlnode.findall("Modules/Module")
        if modules:
            fold = parent.addFold(label=f"Module timeline ({len(modules)})", initiallyOpen=False)
            fold.addPre(text="\n".join(f"{m.get('t', ''):>7}s  {m.get('name', '')}" for m in modules))

    def drawSummaries(self, parent):
        if self.xmlnode.find("Summaries") is None:
            return
        fold = parent.addFold(label="Phaser's summary, module by module", brief="Summary",
                              initiallyOpen=False)
        for node in self.xmlnode.findall("Summaries/Summary"):
            sub = fold.addFold(label=node.get("module", "?"), initiallyOpen=False)
            sub.addPre(text=node.text or "")

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
