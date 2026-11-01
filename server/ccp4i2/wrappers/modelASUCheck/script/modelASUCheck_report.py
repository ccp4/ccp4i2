from ccp4i2.report.CCP4ReportParser import Report
from ccp4i2.report.svg import inline_svg


def coverage_from_alignment(align_seq, align_model, align_match):
    """Per-residue coverage of a construct by a model chain, from the three
    gapped alignment strings the sequence check writes.

    Returns ``(length, modelled, mismatches, unmodelled_runs)``: the
    construct length, the set of construct positions the model covers with
    the same residue, the set it covers with a different residue, and the
    runs ``[(first, last), ...]`` of construct positions no model residue
    covers. Positions are 1-based on the construct.
    """
    modelled, mismatches = set(), set()
    pos = 0
    for s, m, k in zip(align_seq or "", align_model or "", align_match or ""):
        if s == "-":
            continue  # a model residue the construct lacks: not on this axis
        pos += 1
        if m == "-":
            continue
        if k == "|":
            modelled.add(pos)
        else:
            mismatches.add(pos)
    length = pos
    runs = []
    start = None
    for p in range(1, length + 1):
        covered = p in modelled or p in mismatches
        if not covered and start is None:
            start = p
        if covered and start is not None:
            runs.append((start, p - 1))
            start = None
    if start is not None:
        runs.append((start, length))
    return length, modelled, mismatches, runs


def describe_runs(runs):
    """``1-10, 275-300, 483-484`` (a single position is written once)."""
    return ", ".join("%d-%d" % r if r[0] != r[1] else "%d" % r[0] for r in runs)


def coverage_track_svg(length, modelled, mismatches, width=600, height=14):
    """An SVG strip along the construct: covered positions filled, mismatches
    in a distinct colour, gaps left empty. Drawn as runs, so a 500-residue
    chain is a handful of rectangles."""
    if length <= 0:
        return ""
    scale = float(width) / length
    parts = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" '
             'viewBox="0 0 %d %d">' % (width, height, width, height),
             '<rect x="0" y="2" width="%d" height="%d" fill="#e0e0e0" stroke="#9e9e9e" '
             'stroke-width="0.5"/>' % (width, height - 4)]
    for positions, colour in ((modelled, "#4caf50"), (mismatches, "#ff9800")):
        for first, last in _runs(sorted(positions)):
            x = (first - 1) * scale
            w = max((last - first + 1) * scale, 0.5)
            parts.append('<rect x="%.2f" y="2" width="%.2f" height="%d" fill="%s"/>'
                         % (x, w, height - 4, colour))
    parts.append("</svg>")
    return "".join(parts)


def _runs(positions):
    runs = []
    for p in positions:
        if runs and runs[-1][1] == p - 1:
            runs[-1][1] = p
        else:
            runs.append([p, p])
    return [tuple(r) for r in runs]


class modelASUCheck_report(Report):
    TASKNAME = 'modelASUCheck'

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super(). __init__(xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return
        self.addAlignSummary(xmlnode, parent=self)
        self.addAlignReport(xmlnode, initiallyOpen=False, parent=self)

    def addAlignSummary(self, xml, parent=None):
        """Per chain, a coverage track along the construct and a caption
        saying which residues are unmodelled.

        A progress bar of matches over length was what this used to draw,
        and it read as "the missing residues are at the C-terminus": for a
        chain missing its tag, an internal loop and two C-terminal residues
        the bar showed one block missing at the end. The track puts every
        gap where it is.
        """
        if parent is None:
            parent = self
        alignNodes = xml.findall('.//SequenceAlignment/Alignment')
        nonMatchedNodes = xml.findall(
            './/SequenceAlignment/NonAlignedModelChains/Alignment')
        if not alignNodes and not nonMatchedNodes:
            return
        summaryFold = parent.addFold(
            label='Sequence alignment summary', brief='Summary',
            initiallyOpen=True)
        for node in alignNodes:
            chain = node.findtext('ChainID', '')
            identity = node.findtext('identity', '')
            length, modelled, mismatches, runs = coverage_from_alignment(
                node.findtext('align_1', ''), node.findtext('align_2', ''),
                node.findtext('align_match', ''))
            if length == 0:
                summaryFold.addText(text='Chain %s: no alignment' % chain, tone='warning')
                continue
            covered = len(modelled) + len(mismatches)
            caption = 'Chain %s: %d of %d residues modelled' % (chain, covered, length)
            if identity:
                caption += ' (%s%% identity)' % identity
            if runs:
                caption += '; unmodelled %s' % describe_runs(runs)
            if mismatches:
                caption += '; %d residue(s) differ from the construct' % len(mismatches)
            summaryFold.addText(text=caption, tone='emphasis' if not runs else 'info')
            summaryFold.append(inline_svg(coverage_track_svg(length, modelled, mismatches)))
        if nonMatchedNodes:
            chains = sorted(set(
                n.findtext('ChainID', '?') for n in nonMatchedNodes))
            summaryFold.addText(text='Unmatched model chains: ' + ', '.join(chains),
                                tone='warning')

    def addAlignReport(self, xml, initiallyOpen=False, parent=None):
        if parent is None:
            parent = self

        xmlNodes = xml.findall('.//SequenceAlignment/Alignment')
        nonMatchedXmlNodes = xml.findall('.//SequenceAlignment/NonAlignedModelChains/Alignment')

        if len(xmlNodes) > 0 or len(nonMatchedXmlNodes)>0:
            parent.addDiv()
            fold = parent.addFold(label='Sequence alignment detail', brief='Alignment Detail', initiallyOpen=initiallyOpen)
            lineLength = 60
            for alignNode in xmlNodes:
                chainID = alignNode.findall("ChainID")[0].text
                matchCount = alignNode.findall("match_count")[0].text
                identity = alignNode.findall("identity")[0].text
                cigar = alignNode.findall("CIGAR")[0].text
                align1 = alignNode.findall("align_1")[0].text
                align2 = alignNode.findall("align_2")[0].text
                alignMatch = alignNode.findall("align_match")[0].text
                labelDiv = fold.addDiv()
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
            labelDiv = fold.addDiv()
            labelDiv.append("<b>Non matched chains from structure</b>")
            labelDiv = fold.addDiv()
            labelDiv.append('<span>The following chains in the structure have not been matched to anything in the AU contents. This could be because they are incorrectly or partially modelled, or have undefined residue types. It is possible that on or more of the chains below should be merged with matched chains above to improve the alignment. For each unmodelled chain, the alignment with each unique sequence in the AU contents is shown below. </span><br/>')
            labelDiv.addDiv()
            for nonAlignNode in nonMatchedXmlNodes:
                chainID = nonAlignNode.findall("ChainID")[0].text
                iseq = nonAlignNode.findall("SeqAUNo")[0].text
                labelDiv.addDiv()
                labelDiv.append('<span><b>Chain '+chainID+'</b> match with sequence number '+iseq+' in AU contents</span>')
                labelDiv = fold.addDiv()
                labelDiv.addDiv()

                matchCount = nonAlignNode.findall("match_count")[0].text
                identity = nonAlignNode.findall("identity")[0].text
                cigar = nonAlignNode.findall("CIGAR")[0].text
                align1 = nonAlignNode.findall("align_1")[0].text
                align2 = nonAlignNode.findall("align_2")[0].text
                alignMatch = nonAlignNode.findall("align_match")[0].text
                labelDiv = fold.addDiv()

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
