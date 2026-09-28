from ccp4i2.report import Report


class pandda_campaign_report(Report):
    TASKNAME = 'pandda_campaign'
    # program.xml is rewritten as PanDDA reports progress, so the page is
    # worth reading while the run is going.
    RUNNING = True

    _STATE_TEXT = {
        'staged': 'Input tree staged; PanDDA is about to start.',
        'running': 'PanDDA is running.',
        'dispatched': 'PanDDA was handed to a run target and is running elsewhere; the job waits '
                      'here until the run is reconciled (Check remote run on the job, or '
                      'manage.py reconcile_dispatch).',
        'finished': 'PanDDA finished and wrote its events table.',
        'partial': 'PanDDA wrote processed datasets but no events table: a partial run. '
                   'Fan-out can still take what is there.',
        'failed': 'PanDDA failed.',
        'empty': 'PanDDA ran to the end but analysed no dataset. Its reasons are below; '
                 'too few comparators is the usual one, and the minimum is a parameter.',
        'stage_only': 'Stage-only run. The input tree is staged; run PanDDA elsewhere with the '
                      'command below, then fan out from its pandda2_out with this job\'s manifest.',
    }

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style="clear:both;")
        if xmlnode is None:
            return
        def text(path, default=''):
            value = xmlnode.findtext(path)
            return value if value is not None else default
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
        table.addData(title='Datasets loaded', data=[text('n_processed', '-')])
        table.addData(title='Datasets analysed', data=[text('n_analysed', '-')])
        table.addData(title='Events', data=[text('n_events', '-')])
        table.addData(title='Wall time (s)', data=[text('wall_seconds', '-')])
        table.addData(title='Staged tree', data=[text('staging_dir')])
        table.addData(title='Output tree', data=[text('out_dir')])
        table.addData(title='Executable', data=[text('executable', '-')])
        table.addData(title='Contract', data=[text('contract')])
        reasons = [r.text for r in xmlnode.findall('reasons/reason') if r.text]
        if reasons:
            summary.addText(text='Why datasets were left unanalysed, in PanDDA\'s words:')
            summary.addPre(text='\n'.join(reasons))
        if text('probe'):
            summary.addPre(text=text('probe'))

        analysis = xmlnode.find('analysis')
        if analysis is not None:
            self._add_analysis(analysis, self._labels(xmlnode))

        datasets = xmlnode.findall('datasets/dataset')
        if datasets:
            fold = self.addFold(label='Datasets', brief='Datasets', initiallyOpen=False)
            table = fold.addTable()
            table.addData(title='xtal', data=[d.get('xtal') for d in datasets])
            table.addData(title='Label', data=[d.get('label') for d in datasets])
            table.addData(title='Dictionary', data=[d.get('dict') for d in datasets])

    # -- the run's own analysis, recapitulated from its tables ---------------

    @staticmethod
    def _labels(xmlnode):
        """xtal -> the dataset's label (usually the member project name)."""
        return {d.get('xtal'): d.get('label') for d in xmlnode.findall('datasets/dataset')}

    def _add_analysis(self, node, labels):
        stats_node = node.find('stats')
        stats = {child.tag: child.text or '' for child in stats_node} if stats_node is not None else {}
        events = node.findall('events/event')
        sites = node.findall('sites/site')
        datasets = node.findall('datasets/dataset')

        fold = self.addFold(label='Analysis', brief='Analysis', initiallyOpen=True)
        # The events table is written at the end of a run, so mid-run there is
        # nothing in `events` however many the datasets have already found.
        # Saying "PanDDA found no events" then is simply untrue, and it sat
        # above a per-dataset column showing several.
        pending = stats.get('events_table') == 'False'
        counted = self._n(stats.get('n_events')) or 0
        if pending:
            fold.addText(text=(
                f"The run is still going. Its datasets report {int(counted)} event(s) so far, counted "
                f"from each dataset as it finishes. PanDDA writes the events table, and with it the "
                f"sites, the scores and which events it judges interesting, only when the whole run "
                f"ends -- so those are empty below until then."))
        elif not events:
            fold.addText(text='PanDDA found no events. The per-dataset table below shows how far each '
                              'dataset got: how many comparators it had, how many candidate events '
                              'survived the size and score filters.')
        else:
            fold.addText(text=f"{stats.get('n_events')} event(s) in {stats.get('n_datasets_with_events')} "
                              f"dataset(s), clustering into {stats.get('n_sites')} site(s). The best "
                              f"hit-in-site probability is {stats.get('best_hit_probability') or '-'} and "
                              f"the best event score {stats.get('best_score') or '-'}. Each event is "
                              "delivered to its dataset's project by fan-out; look at it there, in the "
                              "scene the receipt job writes.")
        # Said before the numbers, because a run processed far below its own
        # data is not a result to read -- it is a run to do again. PanDDA
        # processes a dataset at the worst resolution among its comparators,
        # so one bad crystal can set the resolution for a whole campaign.
        # Numbers, not the text they arrive as: `stats` is read straight off
        # the XML, so every value here is a string and arithmetic on one
        # raised rather than rendering (the whole Analysis section became a
        # REPORT_GENERATION_FAILED the moment a run had both resolutions).
        processing = self._n(stats.get('median_resolution'))
        best_input = self._n(stats.get('best_input_resolution'))
        dragged = int(self._n(stats.get('n_input_better_than_processing')) or 0)
        if processing and best_input and processing > best_input + 0.5 and dragged:
            fold.addText(text=(
                f"WARNING: processed at {processing:.2f} A, though {dragged} of its datasets are "
                f"better than that and the best is {best_input:.2f} A. PanDDA processes a dataset at "
                f"the worst resolution among its comparators, and takes comparators regardless of "
                f"resolution until it has max_shell_datasets of them -- so on a campaign smaller than "
                f"that, one bad crystal sets the resolution for every dataset. Events found at "
                f"{processing:.2f} A say little about fragments. Remove the low-resolution datasets, "
                f"or lower MAX_SHELL_DATASETS below the number of good ones, and run it again."),
                style="color:#b71c1c; font-weight:600;")

        headline = fold.addTable(transpose=True)
        headline.addData(title='Datasets analysed', data=[stats.get('n_analysed', '-')])
        headline.addData(title='Median processing resolution (A)', data=[stats.get('median_resolution') or '-'])
        headline.addData(title='Best dataset resolution (A)', data=[stats.get('best_input_resolution') or '-'])
        headline.addData(title='Events', data=[stats.get('n_events', '-')])
        headline.addData(title='Datasets with events', data=[stats.get('n_datasets_with_events', '-')])
        headline.addData(title='Sites', data=[stats.get('n_sites', '-')])
        headline.addData(title='Events PanDDA flagged interesting', data=[stats.get('n_interesting', '-')])
        headline.addData(title='Best hit-in-site probability', data=[stats.get('best_hit_probability') or '-'])
        headline.addData(title='Best event score', data=[stats.get('best_score') or '-'])

        if sites:
            self._add_site_chart(fold, sites)
        self._add_histograms(fold, node)

        if events:
            table = fold.addTable()
            table.addData(title='Dataset', data=[labels.get(e.get('dtag'), e.get('dtag')) for e in events])
            table.addData(title='xtal', data=[e.get('dtag') for e in events])
            table.addData(title='Event', data=[e.get('event_idx') for e in events])
            table.addData(title='Site', data=[e.get('site_idx') for e in events])
            table.addData(title='1-BDC', data=[self._f(e.get('event_fraction')) for e in events])
            table.addData(title='Z peak', data=[self._f(e.get('z_peak')) for e in events])
            table.addData(title='Cluster size', data=[self._i(e.get('cluster_size')) for e in events])
            table.addData(title='Event score', data=[self._f(e.get('score')) for e in events])
            table.addData(title='Hit prob.', data=[self._f(e.get('hit_probability')) for e in events])
            table.addData(title='Build score', data=[self._f(e.get('build_score')) for e in events])
            table.addData(title='RSCC', data=[self._f(e.get('rscc')) for e in events])
            table.addData(title='Resolution (A)', data=[self._f(e.get('resolution')) for e in events])
            table.addData(title='R-free', data=[self._f(e.get('r_free')) for e in events])
            table.addData(title='Interesting', data=[e.get('interesting') or '-' for e in events])
            table.addData(title='Centroid', data=[self._centroid(e) for e in events])

        if sites:
            fold.addText(text='Sites: where the events cluster. Several datasets at one site is the '
                              'classic fragment-screen signal; one site with many events from one '
                              'dataset is more often a modelling artefact.')
            table = fold.addTable()
            table.addData(title='Site', data=[s.get('site_idx') for s in sites])
            table.addData(title='Events', data=[s.get('n_events') for s in sites])
            table.addData(title='Datasets', data=[s.get('n_datasets') for s in sites])
            table.addData(title='Interesting', data=[s.get('n_interesting') for s in sites])
            table.addData(title='Best event score', data=[self._f(s.get('best_score')) for s in sites])
            table.addData(title='Best hit prob.', data=[self._f(s.get('best_hit_probability')) for s in sites])
            table.addData(title='Centroid', data=[s.get('centroid') or '-' for s in sites])

        if datasets:
            per = self.addFold(label='Per dataset', brief='Per dataset', initiallyOpen=False)
            per.addText(text='How far each dataset got. Candidate events are counted before and after '
                             'PanDDA\'s size and score filters; a dataset with candidates but no events '
                             'lost them to a filter, not to the search.')
            table = per.addTable()
            table.addData(title='Dataset', data=[labels.get(d.get('dtag'), d.get('dtag')) for d in datasets])
            table.addData(title='xtal', data=[d.get('dtag') for d in datasets])
            table.addData(title='Analysed', data=[self._yes(d.get('analysed')) for d in datasets])
            table.addData(title='Resolution (A)', data=[self._f(d.get('resolution')) for d in datasets])
            table.addData(title='Comparators', data=[d.get('n_comparators') or '-' for d in datasets])
            table.addData(title='Models', data=[d.get('n_models') or '-' for d in datasets])
            table.addData(title='Selected model', data=[d.get('selected_model') or '-' for d in datasets])
            table.addData(title='Candidates', data=[d.get('n_initial_events') or '-' for d in datasets])
            table.addData(title='After size filter', data=[d.get('n_size_filtered_events') or '-' for d in datasets])
            table.addData(title='After score filter', data=[d.get('n_score_filtered_events') or '-' for d in datasets])
            table.addData(title='Events', data=[d.get('n_events') or '-' for d in datasets])
            table.addData(title='R-work', data=[self._f(d.get('r_work')) for d in datasets])
            table.addData(title='R-free', data=[self._f(d.get('r_free')) for d in datasets])
            table.addData(title='Best hit prob.', data=[self._f(d.get('best_hit_probability')) for d in datasets])

    def _add_site_chart(self, parent, sites):
        graph = parent.addFlotGraph(title='Events per site', xmlnode=self.xmlnode,
                                    style='height:220px;width:45%;float:left;')
        graph.addData(title='Site', data=[s.get('site_idx') for s in sites])
        graph.addData(title='Events', data=[s.get('n_events') for s in sites])
        graph.addData(title='Datasets', data=[s.get('n_datasets') for s in sites])
        for title, col, colour in (('Events per site', 2, '#1976d2'), ('Datasets per site', 3, '#388e3c')):
            plot = graph.addPlotObject()
            plot.append('title', title)
            plot.append('plottype', 'xy')
            plot.append('xlabel', 'Site')
            plot.append('ylabel', 'Count')
            plot.append('xintegral', 'true')
            bar = plot.append('barchart', col=1, tcol=col)
            bar.append('colour', colour)

    _HISTOGRAMS = (
        ('event_fraction', 'Event fraction (1-BDC)', 'Events'),
        ('hit_probability', 'Hit-in-site probability', 'Events'),
        ('resolution', 'Processing resolution (A)', 'Datasets'),
        ('r_free', 'R-free', 'Datasets'),
    )

    def _add_histograms(self, parent, node):
        """Reinspect's summary histograms, binned server-side because the
        report viewer draws bar charts but not histograms. One graph per
        histogram: the viewer pads a shared data grid with '-' where columns
        differ in length, and a bar at x='-' is not a bar."""
        for name, title, unit in self._HISTOGRAMS:
            bins = node.findall(f"histograms/histogram[@name='{name}']/bin")
            if not bins:
                continue
            graph = parent.addFlotGraph(title=title, xmlnode=self.xmlnode,
                                        style='height:220px;width:45%;float:left;')
            graph.addData(title=title, data=[b.get('centre') for b in bins])
            graph.addData(title=unit, data=[b.get('count') for b in bins])
            plot = graph.addPlotObject()
            plot.append('title', title)
            plot.append('plottype', 'xy')
            plot.append('xlabel', title)
            plot.append('ylabel', unit)
            bar = plot.append('barchart', col=1, tcol=2)
            bar.append('colour', '#1976d2')
        parent.addDiv(style='clear:both;')

    @staticmethod
    def _n(value):
        """``value`` as a float, or None when it is not a number.

        The report reads its stats as XML text; anything doing arithmetic or
        ``:.2f`` formatting with one needs it as a number first."""
        try:
            return float(value)
        except (TypeError, ValueError):
            return None

    @staticmethod
    def _f(value):
        try:
            return f'{float(value):.3f}'
        except (TypeError, ValueError):
            return '-'

    @staticmethod
    def _i(value):
        try:
            return str(int(float(value)))
        except (TypeError, ValueError):
            return '-'

    @staticmethod
    def _yes(value):
        return {'true': 'yes', 'false': 'no'}.get(value or '', '-')

    @classmethod
    def _centroid(cls, event):
        parts = [event.get(k) for k in ('x', 'y', 'z')]
        if not all(parts):
            return '-'
        return ' '.join(f'{float(p):.1f}' for p in parts)
