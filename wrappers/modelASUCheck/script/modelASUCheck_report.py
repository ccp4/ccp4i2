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
        nonMatchedXmlNodes = xml.findall('.//SequenceAlignment/NonAlignedModelChains/Alignment')

        if len(xmlNodes) > 0 or len(nonMatchedXmlNodes)>0:
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

                if alignMatch is None:
                    fold.append('<span style="color:orange">No chain in the model is matched to chain: '+chainID+' defined in the AU contents. This might mean that the chain defined in the AU contents has not been modelled or that a model chain has not had its residue types defined.</span>')
                    continue

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

        if len(nonMatchedXmlNodes)>0:
            labelDiv = fold.addDiv(style="clear:both;")
            labelDiv.append("<b>Non matched chains from structure</b>")
            labelDiv = fold.addDiv(style="clear:both;")
            labelDiv.append('<span>The following chains in the structure have not been matched to anything in the AU contents. This could be because they are incorrectly or partially modelled, or have undefined residue types. It is possible that on or more of the chains below should be merged with matched chains above to improve the alignment. For each unmodelled chain, the alignment with each unique sequence in the AU contents is shown below. </span><br/>')
            labelDiv.addDiv(style="clear:both;")
            for nonAlignNode in nonMatchedXmlNodes:
                chainID = nonAlignNode.findall("ChainID")[0].text
                iseq = nonAlignNode.findall("SeqAUNo")[0].text
                labelDiv.addDiv(style="clear:both;")
                labelDiv.append('<span><b>Chain '+chainID+'</b> match with sequence number '+iseq+' in AU contents</span>')
                labelDiv = fold.addDiv(style="clear:both;")
                labelDiv.addDiv(style="clear:both;")
                labelDiv.append('<br/>')

                matchCount = nonAlignNode.findall("match_count")[0].text
                identity = nonAlignNode.findall("identity")[0].text
                cigar = nonAlignNode.findall("CIGAR")[0].text
                align1 = nonAlignNode.findall("align_1")[0].text
                align2 = nonAlignNode.findall("align_2")[0].text
                alignMatch = nonAlignNode.findall("align_match")[0].text
                labelDiv = fold.addDiv(style='border:0px solid black; width:700px; overflow:auto;')

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
