"""Report for phaser_ep_auto_phil: hands, sites, figures of merit, and the
substructure-completion cycles, all from the recorded run."""
from ccp4i2.wrappers.phaser_phil.script.phaser_report_base import PhaserReportBase


class phaser_ep_auto_phil_report(PhaserReportBase):
    TASKNAME = "phaser_ep_auto_phil"

    def drawVerdict(self, parent, running):
        overall = self.xmlnode.find("Overall")
        completion = self.xmlnode.find("Completion")
        if overall is not None:
            fom = overall.findtext("fom")
            llg = completion.get("llg") if completion is not None else None
            parent.addText(text=f"Overall figure of merit {fom}" + (f", log-likelihood {llg}" if llg else ""),
                           style="font-weight:bold;")
        if completion is not None and completion.get("converged") == "False":
            parent.addText(text="Substructure completion did not converge within the cycles allowed",
                           style="color:orange;")

    def drawResults(self, parent, running):
        hands = self.xmlnode.find("Hands")
        if hands is not None and len(hands):
            fold = parent.addFold(label="Hands", brief="Hands", initiallyOpen=True)
            table = fold.addTable(xmlnode=hands, select="Hand", style="width:500px;",
                                  outputXml=self.outputXml, internalId="PhaserHands")
            for title, select in (("Hand", "Number"), ("LLG", "LLG"), ("FOM", "FOM"),
                                  ("Coordinates", "PDB"), ("Phases", "MTZ")):
                table.addData(title=title, select=select)
            bins = hands.find("Hand/FOMByResolution")
            if bins is not None and len(bins):
                sub = fold.addFold(label="Figure of merit by resolution (hand 1)", initiallyOpen=False)
                table = sub.addTable(xmlnode=bins, select="Bin", style="width:400px;",
                                     outputXml=self.outputXml, internalId="PhaserFomBins")
                for title, select in (("Low", "LoRes"), ("High", "HiRes"), ("N", "Number"), ("FOM", "FOM")):
                    table.addData(title=title, select=select)
        sites = self.xmlnode.find("Sites")
        if sites is not None and len(sites):
            fold = parent.addFold(label=f"Substructure ({len(sites)} sites)", brief="Sites",
                                  initiallyOpen=True)
            table = fold.addTable(xmlnode=sites, select="Site", style="width:650px;",
                                  outputXml=self.outputXml, internalId="PhaserSites")
            for title, select in (("#", "Number"), ("Element", "Element"), ("Fractional x y z", "Frac"),
                                  ("Occupancy", "Occupancy"), ("B", "B"), ("History", "History")):
                table.addData(title=title, select=select)
        completion = self.xmlnode.find("Completion")
        fold = parent.addFold(label="Substructure completion", brief="Completion",
                              initiallyOpen=not running)
        cycles = completion.findall("Cycle") if completion is not None else []
        if not cycles:
            fold.addText(text="No completion cycle recorded yet" if running
                         else "No completion cycle recorded")
        for c in cycles:
            deleted = c.get("deleted")
            text = (f"Cycle {c.get('cycle')}: LLG {c.get('llg', '?')}, FOM {c.get('fom', '?')}, "
                    f"{c.get('added', '0')} site(s) found" + (f", deleted #{deleted}" if deleted else ""))
            if c.get("converged") == "True":
                text += " -- converged"
            fold.addText(text=text)
            fold.addDiv(style="clear:both;")
        self.drawTimeline(fold)
