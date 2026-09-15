from ccp4i2.report import Report


class molrep_map_report(Report):
    TASKNAME = 'molrep_map'
    RUNNING = False

    # molrep.doc peak-table columns, as scraped (tag name -> display title).
    _PEAK_COLUMNS = [
        ("RF", "RF"), ("TF", "TF"), ("TF_sig", "TF/sig"),
        ("TFcntrst", "TF contrast"), ("PFind", "PFind"), ("PF", "PF"),
        ("PFmod", "PFmod"), ("wRfac", "wRfac"), ("Score", "Score"),
        ("Cntrst", "Contrast"), ("for", "FOM"),
    ]

    _CONFIDENCE_TEXT = {
        "confident": ("A clear winner: the recommended hand's model fits the map "
                      "well and the two hands are well separated."),
        "ambiguous": ("The two hands are hard to separate. The recommendation is a "
                      "best guess -- inspect both before committing."),
        "weak": ("Both hands fit the map poorly, so the hand may not be "
                 "determinable from this data. Inspect both."),
        "single": ("Only one hand produced a placement."),
        "none": ("No placement was scored."),
    }

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        results = self.addResults()
        self._add_recommendation(results)

        for hand in ['Original', 'Flipped']:
            fold = results.addFold(label=f'{hand} hand')
            node = self.xmlnode.findall(f'./{hand}')
            if node and node[0].get('placed') == 'false':
                fold.addText(text=f'No model was placed in the {hand.lower()} map.')
                continue
            table = fold.addTable(
                select=f"./{hand}/MolrepResult/RFpeaks/RFpeak")
            for tag, title in self._PEAK_COLUMNS:
                table.addData(title=title, select=tag)

            graph = fold.addFlotGraph(
                title="Best TF peak vs RF peak No",
                select=f"./{hand}/MolrepResult/RFpeaks/RFpeak")
            graph.addData(title="RF_peak_No", select="RF")
            graph.addData(title="Score", select="Score")
            graph.addData(title="TF_sig", select="TF_sig")
            graph.addPlot(plot='''<plot>
        <title>Best TF peak score vs RF peak No</title>
        <plottype>xy</plottype>
        <plotline xcol="1" ycol="2">
        <linestyle>.</linestyle>
        <markeredgewidth>0</markeredgewidth>
        <colour>blue</colour>
        </plotline>
        </plot>''')
            graph.addPlot(plot='''<plot>
        <title>TF/sig(TF) vs RF peak No</title>
        <plottype>xy</plottype>
        <plotline xcol="1" ycol="3">
        <linestyle>.</linestyle>
        <markeredgewidth>0</markeredgewidth>
        <colour>red</colour>
        </plotline>
        </plot>''')

        self.addTaskReferences()

    def _add_recommendation(self, results):
        rec = self.xmlnode.findall('.//recommendation')
        if not rec:
            return
        rec = rec[0]
        hand = rec.get('hand', 'Original')
        confidence = rec.get('confidence', 'none')
        label = 'original' if hand == 'Original' else 'inverted'

        cc_o = rec.get('cc_original')
        cc_f = rec.get('cc_flipped')
        cc_bits = []
        if cc_o is not None:
            cc_bits.append(f'Original {cc_o}')
        if cc_f is not None:
            cc_bits.append(f'Flipped {cc_f}')
        cc_line = ('Real-space map-model CC: ' + ', '.join(cc_bits) + '. ') if cc_bits else ''

        # addText renders as escaped plain text (the frontend shows tags/entities
        # literally), so keep this plain -- no HTML markup or entities.
        note = self._CONFIDENCE_TEXT.get(confidence, '')
        results.addText(text=(
            f'Recommended hand: {hand} ({label} map). Confidence: {confidence}. '
            f'{cc_line}{note} Both hands are provided below; the recommended hand '
            'carries any half maps through for cross-validated refinement.'))

    def addTaskReferences(self):
        try:
            super().addTaskReferences()
        except Exception:
            pass
