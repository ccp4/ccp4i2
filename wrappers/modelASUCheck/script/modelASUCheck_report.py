"""
    modelASUCheck_report.py: CCP4 GUI Project

    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.

    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

from report.CCP4ReportParser import Report


class modelASUCheck_report(Report):
    TASKNAME = 'modelASUCheck'

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super(). __init__(xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addAlignReport(xmlnode, initiallyOpen=True)

    def addAlignReport(self, xml, initiallyOpen=False, parent=None):
        if parent is None:
            parent = self

        xmlNodes = xml.findall('.//SequenceAlignment/Alignment')
        if len(xmlNodes) > 0:
            parent.addDiv(style="clear:both;")
            fold = parent.addFold(label='Sequence alignment information', brief='Seq. Align. Info.', initiallyOpen=initiallyOpen)
            lineLength = 60
            for alignNode in xmlNodes:
                chainID = alignNode.findall("ChainID")[0].text
                matchCount = alignNode.findall("match_count")[0].text
                identity = alignNode.findall("identity")[0].text
                cigar = alignNode.findall("CIGAR")[0].text
                align1 = alignNode.findall("align_1")[0].text
                align2 = alignNode.findall("align_2")[0].text
                alignMatch = alignNode.findall("align_match")[0].text
                labelDiv = fold.addDiv(style='border:0px solid black; width:700px; overflow:auto;')
                labelDiv.append("<b>Chain: "+chainID+"</b>")

                tableText = "<table>\n"
                tableText += f"<tr><th>Match count</th><td>{matchCount}</td></tr>"
                tableText += f"<tr><th>Identity</th><td>{identity}</td></tr>"
                tableText += f"<tr><th>CIGAR</th><td>{cigar}</td></tr>"
                tableText += "</table>\n"
                fold.append(tableText)

                preText = ""
                for i in range(0, len(alignMatch), lineLength):
                    preText += str(i + 1)
                    preText += (60 - len(str(i) + str(i + lineLength))) * " "
                    preText += str(min(i + lineLength, len(alignMatch)))
                    preText += "\n"
                    preText += align1[i : i + lineLength] + "\n"
                    preText += alignMatch[i : i + lineLength] + "\n"
                    preText += align2[i : i + lineLength] + "\n"
                    preText += "\n"
                fold.append('<pre>'+preText+'</pre>')
