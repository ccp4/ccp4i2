from ccp4i2.report import Report


class pandda_campaign_report(Report):
    TASKNAME = 'pandda_campaign'
    # program.xml is rewritten as PanDDA reports progress, so the page is
    # worth reading while the run is going.
    RUNNING = True

    _STATE_TEXT = {
        'staged': 'Input tree staged; PanDDA is about to start.',
        'running': 'PanDDA is running.',
        'finished': 'PanDDA finished and wrote its events table.',
        'partial': 'PanDDA wrote processed datasets but no events table: a partial run. '
                   'Fan-out can still take what is there.',
        'failed': 'PanDDA failed.',
        'stage_only': 'Stage-only run. The input tree is staged; run PanDDA elsewhere with the '
                      'command below, then fan out from its pandda2_out with this job\'s manifest.',
    }

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style="clear:both;")
        if xmlnode is None:
            return
        text = lambda path, default='': (
            xmlnode.findtext(path) if xmlnode.findtext(path) is not None else default)
        state = text('state', 'unknown')
        self.addText(text=self._STATE_TEXT.get(state, state),
                     style='font-weight:bold;' + ('color:#c62828;' if state == 'failed' else ''))
        if state == 'failed' and text('failure'):
            self.addText(text=f"Classified as: {text('failure')} (see the job's diagnostics for what to do).")

        progress = xmlnode.find('progress')
        if state == 'running':
            if progress is not None:
                self.addText(text=f"Dataset {progress.get('done')} of {progress.get('total')}.")
            else:
                self.addText(text='Progress unknown: this PanDDA build does not report it.')

        summary = self.addFold(label='Run', brief='Run', initiallyOpen=True)
        table = summary.addTable(transpose=True)
        table.addData(title='Datasets staged', data=[text('n_datasets')])
        table.addData(title='Datasets processed', data=[text('n_processed', '-')])
        table.addData(title='Events', data=[text('n_events', '-')])
        table.addData(title='Wall time (s)', data=[text('wall_seconds', '-')])
        table.addData(title='Staged tree', data=[text('staging_dir')])
        table.addData(title='Output tree', data=[text('out_dir')])
        table.addData(title='Executable', data=[text('executable', '-')])
        table.addData(title='Contract', data=[text('contract')])
        if text('probe'):
            summary.addPre(text=text('probe'))

        datasets = xmlnode.findall('datasets/dataset')
        if datasets:
            fold = self.addFold(label='Datasets', brief='Datasets', initiallyOpen=False)
            table = fold.addTable()
            table.addData(title='xtal', data=[d.get('xtal') for d in datasets])
            table.addData(title='Label', data=[d.get('label') for d in datasets])
            table.addData(title='Dictionary', data=[d.get('dict') for d in datasets])
