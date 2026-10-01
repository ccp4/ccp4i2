import os

from ccp4i2.report import Report


# Each sentence group in a div of its own: addText makes an inline span, and
# two in a row run together with no space between them.
class slicendice_report(Report):
    TASKNAME = 'slicendice'
    USEPROGRAMXML = True
    RUNNING = True
    CSS_VERSION = '0.1.0'

    def __init__(self, xmlnode=None, jobInfo={}, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, cssVersion=self.CSS_VERSION, **kw)
        results = self.addResults()
        results.addDiv().addText(text="SliceNDice prepares a predicted model for molecular replacement "
                        "(B-factors from its confidence, low-confidence residues removed), "
                        "slices it into rigid regions, and tries each set of slices in Phaser, "
                        "refining every placement with Refmac.")
        # The report is built from program.xml alone. (It read
        # slicendice_results.json at once, so the running report -- before
        # the file exists -- failed; and it pasted the log in unescaped.)
        best = self.xmlnode.find(".//RunInfo/Best") if self.xmlnode is not None else None
        if best is None:
            results.addDiv().addText(text="Running: results will appear here when the "
                            "placements have been refined.")
        else:
            self.summary(results, best)
        logf = os.path.join(jobInfo.get("fileroot", ""), "slicendice_0", "slicendice.log")
        if os.path.isfile(logf):
            fold = results.addFold(label="SliceNDice log", initiallyOpen=best is None)
            with open(logf, encoding="utf-8", errors="replace") as f:
                fold.addPre(text=f.read())

    def summary(self, parent, best):
        run = self.xmlnode.find(".//RunInfo")
        solved = best.findtext("Solved") == "True"
        rfree, r = float(best.findtext("RFree")), float(best.findtext("R"))
        n = best.findtext("bid")
        partial = best.findtext("Partial") == "True"
        if solved and partial:
            parent.addDiv().addText(text="Partly solved: the placement from %s split%s refined to R %.3f, "
                           "R-free %.3f, which passes SliceNDice's test (both below 0.45), but not "
                           "every search model was placed (below)." % (
                               n, "" if n == "1" else "s", r, rfree))
        elif solved:
            parent.addDiv().addText(text="Solved: the placement from %s split%s refined to R %.3f, "
                           "R-free %.3f." % (n, "" if n == "1" else "s", r, rfree))
        else:
            parent.addDiv().addText(text="No solution. The best placement (%s split%s) refined only to "
                           "R %.3f, R-free %.3f; SliceNDice counts a placement as a solution "
                           "when both are below 0.45." % (n, "" if n == "1" else "s", r, rfree))
        splits = {s.get("id"): [m.text for m in s.findall("Model")] for s in run.findall("Split")}
        sols = run.findall("Sol")
        table = parent.addTable()
        table.addData(title="Splits", data=[s.findtext("SolID") for s in sols])
        table.addData(title="Search models (residues)",
                      data=["; ".join(splits.get(s.findtext("SolID"), [])) for s in sols])
        # SliceNDice's LLG and TFZ per split are those of Phaser's last listed
        # solution, not its best: labelled as SliceNDice's, beside each piece's own.
        for title, tag in (("LLG (SliceNDice)", "llg"), ("TFZ (SliceNDice)", "tfz"),
                           ("R", "srf"), ("R-free", "sre")):
            table.addData(title=title, data=[s.findtext(tag) for s in sols])
        # Each search model's placement (Phaser's best solution, in the order
        # placed). SliceNDice's own TFZ is one number per split; a lobe placed
        # at TFZ 26 and one not placed at all (TFZ 6, LLG +18) read as 6.
        for sol in sols:
            comps = sol.findall("Component")
            if not comps:
                continue
            models = splits.get(sol.findtext("SolID"), [])
            parts, unsure, previous = [], [], 0.0
            for c in comps:
                i = int(c.get("cluster")) if c.get("cluster", "").isdigit() else None
                where = ("residues " + models[i]) if i is not None and i < len(models) else "model %s" % c.get("cluster")
                gain = float(c.get("llg") or 0) - previous
                previous = float(c.get("llg") or 0)
                clashes = c.get("clashes")
                parts.append("%s: TFZ %s, LLG +%.0f, %s clash%s" % (
                    where, c.get("tfz"), gain, clashes, "" if clashes == "1" else "es"))
                if c.get("tfz") and float(c.get("tfz")) < 8:
                    unsure.append(where)
            text = "%s split%s, each search model as Phaser placed it: %s." % (
                sol.findtext("SolID"), "" if sol.findtext("SolID") == "1" else "s", "; ".join(parts))
            if unsure:
                text += (" A TFZ below 8 is not a clear placement: %s probably did not find "
                         "its place, even where the refined R-free passes." % " and ".join(unsure))
            parent.addDiv().addText(text=text)
        tried = {s.findtext("SolID") for s in sols}
        untried = sorted(set(splits) - tried)
        if untried:
            parent.addDiv().addText(text="Not tried in molecular replacement: %s split%s (%s). "
                           "SliceNDice 0.1.3 runs MR on one of the splits it makes; to try "
                           "a particular number of splits, set the minimum and maximum "
                           "splits to it." % (
                               ", ".join(untried), "" if len(untried) == 1 else "s",
                               "; ".join("; ".join(splits[u]) for u in untried)))
