from ccp4i2.report import Report


class molrep_map_report(Report):
    TASKNAME = 'molrep_map'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        results = self.addResults()

        rec = self.xmlnode.findall('.//recommendation')
        if rec:
            hand = rec[0].get('hand', 'Original')
            label = 'original' if hand == 'Original' else 'inverted'
            results.addText(text=(
                f'Recommended hand: <b>{hand}</b> ({label} map). '
                'Both hands are provided below for inspection; the recommended '
                'hand carries any half maps through for cross-validated refinement.'))

        for hand in ['Original', 'Flipped']:
            fold = results.addFold(label=f'{hand} hand')
            node = self.xmlnode.findall(f'./{hand}')
            if node and node[0].get('placed') == 'false':
                fold.addText(text=f'No model was placed in the {hand.lower()} map.')
                continue
            table = fold.addTable(
                select=f"./{hand}/MolrepResult/RFpeaks/RFpeak")
            for title, select in [["RF", "RF"], ["TF", "TF"],
                                  ["Tf_sig", "Tf_sig"], ["TFcntrst", "TFcntrst"],
                                  ["PFind", "PFind"], ["PF", "PF"],
                                  ["PFmin", "PFmin"], ["wRfac", "wRfac"],
                                  ["Score", "Score"], ["Cntrst", "Cntrst"],
                                  ["For", "for"]]:
                table.addData(title=title, select=select)

            graph = fold.addGraph(
                title="Best TF peak vs RF peak No",
                select=f"./{hand}/MolrepResult/RFpeaks/RFpeak")
            graph.addData(title="RF_peak_No", select="RF")
            graph.addData(title="Score", select="Score")
            graph.addData(title="Tf_sig", select="Tf_sig")
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

    def addTaskReferences(self):
        try:
            super().addTaskReferences()
        except Exception:
            pass
