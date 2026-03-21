import os
import xml.etree.ElementTree as etree

from ccp4i2.report import Report


class matthews_report(Report):
    TASKNAME = 'matthews'
    RUNNING = False
    USEPROGRAMXML = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if xmlnode is None:
            return

        fold = self.addFold(label="Matthews coefficient analysis", initiallyOpen=True)

        # Cell volume
        cellVolume = xmlnode.findall(".//cellVolume")
        if cellVolume:
            fold.append('<p>Cell volume: {0:.1f} &#197;<sup>3</sup></p>'.format(
                float(cellVolume[0].text)))

        # Probability table
        compositions = xmlnode.findall(".//matthewsCompositions/composition")
        if compositions:
            table = fold.addTable(
                select=".//matthewsCompositions/composition",
                xmlnode=xmlnode,
            )
            for title, select in [
                ["Molecules in ASU", "nMolecules"],
                ["Solvent %", "solventPercentage"],
                ["Matthews coefficient", "matthewsCoeff"],
                ["Probability", "matthewsProbability"],
            ]:
                table.addData(title=title, select=select)

        # Summary
        summary = xmlnode.find(".//summary")
        if summary is not None:
            nmol = summary.findtext("mostLikelyNmol", "?")
            solvent = summary.findtext("solventContent", "?")
            fold.append(
                '<p><b>Most likely: {0} molecule(s) in ASU, '
                '{1}% solvent</b></p>'.format(nmol, solvent)
            )
