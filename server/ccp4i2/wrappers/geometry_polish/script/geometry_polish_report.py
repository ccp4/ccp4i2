"""Report for geometry_polish: what MolProbity said before and after, what each
stage did, every accepted change, and what is left for a human."""
from ccp4i2.report.CCP4ReportParser import Report

SUMMARY_ROWS = (
    ("Rotamer outliers", "rot_outliers"),
    ("Rotamer outliers (%)", "rot_pct"),
    ("Ramachandran outliers", "rama_outliers"),
    ("Ramachandran outliers (%)", "rama_pct"),
    ("Ramachandran favoured (%)", "rama_favored_pct"),
    ("Clashscore", "clashscore"),
    ("Clashes", "n_clashes"),
)
STAGE_LABELS = {"rotamers": "Rotamer outliers", "ramachandran": "Ramachandran outliers",
                "clashes": "Side-chain clashes"}


def _residue(node):
    """'SER 100A' from a change or review node."""
    return "{} {}{}".format(node.get("resname", ""), node.get("resno", ""), node.get("icode", "")).strip()


def _attr(node, name):
    if node is None:
        return ""
    return node.get(name, "")


class geometry_polish_report(Report):
    TASKNAME = "geometry_polish"
    RUNNING = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus in ("Running", "Running remotely"):
            self.addText(text="The polish is running; the tables below fill in as residues are decided.",
                         tone="emphasis")
        self.addSummary()
        self.addStages()
        self.addChanges()
        self.addReview()

    def addSummary(self):
        before = self.xmlnode.find(".//summary/before")
        after = self.xmlnode.find(".//summary/after")
        if before is None:
            return
        fold = self.addFold(label="MolProbity before and after", brief="Summary", initiallyOpen=True)
        if after is None:
            fold.addText(text="The 'after' column appears when the polish has finished.", tone="muted")
        table = fold.addTable()
        table.addData(title="", data=[label for label, _ in SUMMARY_ROWS])
        table.addData(title="Input", data=[_attr(before, key) for _, key in SUMMARY_ROWS])
        table.addData(title="Polished", data=[_attr(after, key) for _, key in SUMMARY_ROWS])

    def addStages(self):
        stages = self.xmlnode.findall(".//stages/stage")
        if not stages:
            return
        fold = self.addFold(label="What each stage did", brief="Stages", initiallyOpen=True)
        table = fold.addTable()
        table.addData(title="Stage", data=[STAGE_LABELS.get(s.get("name"), s.get("name")) for s in stages])
        table.addData(title="Targets", data=[s.get("targets", "") for s in stages])
        table.addData(title="Fixed", data=[s.get("fixed", "") for s in stages])
        table.addData(title="Left for review", data=[s.get("left", "") for s in stages])
        table.addData(title="Seconds", data=[s.get("seconds", "") for s in stages])

    def addChanges(self):
        changes = self.xmlnode.findall(".//changes/change")
        fold = self.addFold(label="Accepted changes ({})".format(len(changes)), brief="Changes",
                            initiallyOpen=False)
        if not changes:
            fold.addText(text="No change has been accepted yet.", tone="muted")
            return
        fold.addText(text="Density is the mean 2mFo-DFc value at the moved atoms, in sigma. "
                          "Every change listed kept its density within the guard and removed the "
                          "MolProbity flag without adding a contact.", tone="info")
        table = fold.addTable()
        table.addData(title="Chain", data=[c.get("chain", "") for c in changes])
        table.addData(title="Residue", data=[_residue(c) for c in changes])
        table.addData(title="Stage", data=[STAGE_LABELS.get(c.get("stage"), c.get("stage", "")) for c in changes])
        table.addData(title="Action", data=[c.get("action", "") for c in changes])
        table.addData(title="Density before", data=[c.get("density_before", "") for c in changes])
        table.addData(title="Density after", data=[c.get("density_after", "") for c in changes])
        table.addData(title="Note", data=[c.get("detail", "") for c in changes])

    def addReview(self):
        items = self.xmlnode.findall(".//review/item")
        fold = self.addFold(label="For a human to look at ({})".format(len(items)), brief="Review",
                            initiallyOpen=bool(items))
        if not items:
            fold.addText(text="Nothing was left unresolved.", tone="success")
            return
        fold.addText(text="These residues stay as they were: no candidate removed the MolProbity flag "
                          "without losing density or adding a clash. Outliers in strong density are "
                          "usually genuine and can be defended at deposition.", tone="warning")
        table = fold.addTable()
        table.addData(title="Kind", data=[i.get("kind", "") for i in items])
        table.addData(title="Chain", data=[i.get("chain", "") for i in items])
        table.addData(title="Residue", data=[_residue(i) for i in items])
        table.addData(title="Detail", data=[i.get("detail", "") for i in items])
