import os

from ccp4i2.report import Report


class editbfac_report(Report):
    TASKNAME = 'editbfac'
    RUNNING = False

    def __init__(self, xmlnode=None,jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.addResults()
        # What was kept and how it was split, from the job's own summary. (It
        # said only "Edit B-factors finished." above the program's log.)
        model = self.xmlnode.find("Model") if self.xmlnode is not None else None
        if model is not None:
            domains = self.xmlnode.findall("Domain")
            # One text element: two in a row ran together with no space.
            text = "%s of %s residues kept: %s." % (
                model.get("residues"), self.xmlnode.findtext("InputResidues"),
                model.get("ranges"))
            if domains:
                text += " Split into %d region%s, each a file of its own:" % (
                    len(domains), "" if len(domains) == 1 else "s")
            parent.addText(text=text)
            if domains:
                table = parent.addTable()
                table.addData(title="Region", data=[d.get("chain") for d in domains])
                table.addData(title="Residue range", data=[d.get("ranges") for d in domains])
                table.addData(title="Residues", data=[d.get("residues") for d in domains])
        logf = os.path.join(self.getJobFolder(), "log.txt")
        if os.path.isfile(logf):
            fold = parent.addFold(label='Log from cctbx', initiallyOpen=model is None)
            with open(logf, encoding="utf-8", errors="replace") as f:
                fold.addPre(text=f.read())
