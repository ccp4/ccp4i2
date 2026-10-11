"""Report for phaser_mr_auto_phil: the recorded run, never an inference."""
from ccp4i2.wrappers.phaser_phil.script.phaser_report_base import PhaserReportBase


class phaser_mr_auto_phil_report(PhaserReportBase):
    TASKNAME = "phaser_mr_auto_phil"

    def drawVerdict(self, parent, running):
        super().drawVerdict(parent, running)
        if not running and not self.xmlnode.findtext("Verdict"):
            parent.addText(text="Phaser reported no verdict", tone="error")
        best = self.xmlnode.find("PhaserCurrentBestSolution")
        if running and best is not None and best.text:
            self.drawCurrentBest(parent, best.text.strip())

    def drawCurrentBest(self, parent, text):
        """Phaser's running best solution, read rather than quoted: one row
        per component placed so far (its RFZ, search TFZ, clashes and the
        LLG once it was in), the space group, and the refined TFZ of the
        whole; the annotation itself under a fold for anyone who wants it.
        It was printed raw: "SOLU SET RFZ=4.6 TFZ=6.2 PAK=0 LLG=40 TFZ==8.1
        SOLU SPAC P 21 21 21 SOLU 6DIM ENSE ...", one unbroken line."""
        from ccp4i2.wrappers.phaser_phil.script.phaser_run import annotation_tokens
        import re
        annotation = text.split("SOLU SET", 1)[1] if "SOLU SET" in text else text
        annotation = annotation.split("SOLU ", 1)[0]
        components, _unknown = annotation_tokens(annotation)
        space_group = re.search(r"SOLU SPAC\s+([^\n]*?)(?=\s+SOLU|\s*$)", text)
        refined = re.findall(r"TFZ==(-?[\d.]+)", annotation)
        parts = []
        if space_group:
            parts.append(f"space group {space_group.group(1).strip()}")
        if components:
            parts.append(f"{len(components)} component{'s' if len(components) != 1 else ''} placed")
        if refined:
            parts.append(f"refined TFZ {refined[-1]}")
        parent.addText(text="Current best solution: " + (", ".join(parts) if parts else text),
                       tone="emphasis")
        if components:
            def column(key):
                def show(v):
                    return str(int(v)) if isinstance(v, float) and v.is_integer() else str(v)
                return [show(c[key]) if c.get(key) is not None else "–" for c in components]
            table = parent.addTable(style="width:480px;", outputXml=self.outputXml,
                                    internalId="PhaserCurrentBestTable")
            table.addData(title="Component", data=[str(i + 1) for i in range(len(components))])
            table.addData(title="RFZ", data=column("RFZ"))
            table.addData(title="TFZ", data=column("TFZ"))
            table.addData(title="Clashes", data=column("PAK"))
            table.addData(title="LLG", data=column("LLG"))
            table.addData(title="Note", data=[", ".join(c.get("notes") or []) for c in components])
        fold = parent.addFold(label="Phaser's annotation", initiallyOpen=False)
        fold.addText(text=text, tone="code")

    def drawResults(self, parent, running):
        if self.xmlnode.find("Solutions") is not None:
            fold = parent.addFold(label="Solutions", brief="Solutions", initiallyOpen=True)
            self.drawSolutions(fold)
        fold = parent.addFold(label="Search strategy", brief="Strategy", initiallyOpen=not running)
        self.drawStrategy(fold, running)

    def drawSolutions(self, parent):
        solutions = self.xmlnode.find("Solutions")
        found = solutions.findall("Solution")
        # Phaser keeps every score on a solution, and one not calculated is
        # 0: a translation-function solution "has" LLG 0.00, R 0.00 and no
        # clashes. Its placements say what was calculated (the annotation
        # carries only those), so a score is shown when one of them has it.
        # R comes with refinement, as LLG does.
        computed = [{key for p in s.findall("Placements/Placement") for key in
                     ("TFZ", "TFZeq", "PAK", "LLG") if p.find(key) is not None}
                    for s in found]

        def column(field, needs):
            return [s.findtext(field) if needs in done and s.findtext(field) else "–"
                    for s, done in zip(found, computed)]

        table = parent.addTable(style="width:700px;", outputXml=self.outputXml,
                                internalId="PhaserSolutionsTable")
        table.addData(title="#", data=[s.findtext("Number") for s in found])
        table.addData(title="Space group", data=[s.findtext("spaceGroup") for s in found])
        table.addData(title="LLG", data=column("LLG", "LLG"))
        table.addData(title="TFZ", data=column("TFZ", "TFZ"))
        table.addData(title="TFZ-equiv", data=column("TFZeq", "TFZeq"))
        table.addData(title="R", data=column("R", "LLG"))
        table.addData(title="Clashes", data=column("PAK", "PAK"))
        table.addData(title="Annotation", data=[s.findtext("Annotation") or "" for s in found])
        for i, solution in enumerate(found):
            placements = solution.find("Placements")
            if placements is None or len(placements) == 0:
                continue
            fold = parent.addFold(label=f"Solution {i + 1}: placements", initiallyOpen=(i == 0))
            table = fold.addTable(style="width:600px;", outputXml=self.outputXml,
                                  internalId=f"PhaserPlacements{i}")
            rows = placements.findall("Placement")
            for title, select in (("RFZ", "RFZ"), ("TFZ", "TFZ"), ("TFZ-equiv", "TFZeq"),
                                  ("Clashes", "PAK"), ("LLG", "LLG")):
                table.addData(title=title, data=[p.findtext(select) or "–" for p in rows])
            table.addData(title="Note", data=[", ".join(n.text or "" for n in p.findall("Note"))
                                              or "" for p in rows])
            history = solution.findtext("History")
            if history:
                fold.addText(text=f"History: {history}")
            unknown = solution.findtext("UnknownTokens")
            if unknown:
                fold.addText(text=f"Annotation tokens not understood: {unknown}",
                             tone="warning")

    def drawStrategy(self, parent, running):
        strategy = self.xmlnode.find("Strategy")
        attempts = strategy.findall("Attempt") if strategy is not None else \
            self.xmlnode.findall("Attempts/Attempt")
        if not attempts and getattr(self, "SEARCHES", True):
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
                parent.addText(text="    " + (warning.text or ""), tone="warning")
            parent.addDiv()
        if strategy is not None and int(strategy.get("unparsed", "0")) > 0:
            parent.addText(text=f"{strategy.get('unparsed')} strategy block(s) in Phaser's "
                                "summary matched no known sentence -- see the summary below",
                           tone="warning")
        self.drawTimeline(parent)
