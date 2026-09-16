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

    # Confidence -> colour for the verdict badge (table cells render inline HTML).
    _CONFIDENCE_COLOUR = {
        "confident": "#2e7d32",   # green
        "ambiguous": "#e65100",   # orange
        "weak": "#c62828",        # red
        "single": "#616161",      # grey
        "none": "#616161",
    }

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        results = self.addResults()
        self._add_recommendation(results)

        for hand in ['Original', 'Inverted']:
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
        """The verdict: a coloured hand-comparison table (map-model CC is the
        headline; molrep score alongside), then the confidence in words.

        Note the asymmetry with the peak tables below: ``addText`` escapes, but
        table *cells* render inline HTML (``_set_cell_content``), so the emphasis
        and colour live in the table, not in free text.
        """
        rec = self.xmlnode.findall('.//recommendation')
        if not rec:
            return
        rec = rec[0]
        hand = rec.get('hand', 'Original')
        confidence = rec.get('confidence', 'none')
        colour = self._CONFIDENCE_COLOUR.get(confidence, '#616161')

        results.addText(text=(
            'Cryo-EM hand assignment: the model is placed into the map and into '
            'its mirror image, and the better real-space fit (map-model '
            'correlation) wins.'))

        hands = ['Original', 'Inverted']
        cc = {'Original': rec.get('cc_original'), 'Inverted': rec.get('cc_inverted')}

        def emphasise(h, s):
            return f'<b>{s}</b>' if h == hand else s

        table = results.addTable(title='Hand assignment')
        table.addData(title='Hand', data=[
            emphasise(h, 'Original (as given)' if h == 'Original'
                      else 'Inverted (mirror image)') for h in hands])
        table.addData(title='Placed', data=[
            self._hand_placed(h) for h in hands])
        table.addData(title='Map-model CC', data=[
            emphasise(h, cc[h] if cc[h] is not None else '-') for h in hands])
        table.addData(title='molrep score', data=[
            self._hand_score(h) for h in hands])
        # Cell HTML is parsed as XML (_set_cell_content), so use numeric char
        # refs only -- named HTML entities like &ndash; fail the parse and the
        # whole cell falls back to escaped text.
        table.addData(title='Verdict', data=[
            (f'<span style="color:{colour};font-weight:bold">&#10003; '
             f'recommended &#8211; {confidence}</span>') if h == hand else ''
            for h in hands])

        results.addText(text=self._CONFIDENCE_TEXT.get(confidence, '')
                        + ' Both hands are provided below; the recommended hand '
                        'carries any half maps through for cross-validated '
                        'refinement.')

    def _hand_placed(self, hand):
        el = self.xmlnode.find(f'./{hand}')
        if el is None:
            return '-'
        if el.get('placed') != 'true':
            return 'no'
        return 'timed out' if el.get('timed_out') == 'true' else 'yes'

    def _hand_score(self, hand):
        el = self.xmlnode.find(f'./{hand}')
        return el.get('score', '-') if el is not None else '-'

    def addTaskReferences(self):
        try:
            super().addTaskReferences()
        except Exception:
            pass
