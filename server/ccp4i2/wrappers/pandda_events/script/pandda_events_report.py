from ccp4i2.report import Report


class pandda_events_report(Report):
    TASKNAME = 'pandda_events'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style="clear:both;")
        if xmlnode is None:
            return

        text = lambda path, default='': (
            xmlnode.findtext(path) if xmlnode.findtext(path) is not None else default)
        dtag = text('dtag')
        shortfalls = [n.text for n in xmlnode.findall('shortfalls/item')]
        incomplete = text('run_incomplete') == 'True'

        if shortfalls:
            self.addText(text=(
                f'{dtag}: receipt is SHORT. PanDDA declared {len(shortfalls)} '
                f'output(s) that are not in the tree: ' + '; '.join(shortfalls) + '. '
                'Everything that did arrive is available below.'),
                style='color:#c62828;font-weight:bold;')
        if incomplete:
            self.addText(text=('This receipt was taken from a run that did not finish; '
                               'the events here are what that run got to.'),
                         style='color:#e65100;')
        if text('events_table_present') == 'False':
            self.addText(text=('The run-level events table was not written, so site '
                               'numbers and hit probabilities are absent.'))

        summary = self.addFold(label='Summary', brief='Summary', initiallyOpen=True)
        counts = summary.addTable(transpose=True)
        counts.addData(title='Dataset', data=[dtag])
        counts.addData(title='Apo model (model of record)', data=[text('apo_model')])
        counts.addData(title='Z-map', data=[text('zmap')])
        counts.addData(title='Merged PanDDA model', data=[text('pandda_model')])
        counts.addData(title='Events expected / event maps delivered',
                       data=[f"{text('counts/events_expected')} / {text('counts/event_maps_delivered')}"])
        counts.addData(title='Poses expected / delivered',
                       data=[f"{text('counts/poses_expected')} / {text('counts/poses_delivered')}"])

        events = xmlnode.findall('events/event')
        if not events:
            summary.addText(text='PanDDA recorded no events for this dataset.')
            return

        fold = self.addFold(label='Events', brief='Events', initiallyOpen=True)
        fold.addText(text=('Every pose is a candidate merged onto the apo model, to be '
                           'judged. Optimal contour is in absolute map units, not sigma.'))
        table = fold.addTable()
        get = lambda node, tag: node.findtext(tag) if node.findtext(tag) is not None else '-'
        table.addData(title='Event', data=[n.get('idx') for n in events])
        table.addData(title='Site', data=[get(n, 'site_idx') for n in events])
        table.addData(title='BDC', data=[self._f(get(n, 'bdc')) for n in events])
        table.addData(title='Event score', data=[self._f(get(n, 'score')) for n in events])
        table.addData(title='Hit prob.', data=[self._f(get(n, 'hit_probability')) for n in events])
        table.addData(title='Build score', data=[self._f(get(n, 'build_score')) for n in events])
        table.addData(title='RSCC', data=[self._f(get(n, 'rscc')) for n in events])
        table.addData(title='Optimal contour', data=[self._f(get(n, 'optimal_contour')) for n in events])
        table.addData(title='Centroid', data=[get(n, 'centroid') for n in events])
        table.addData(title='Event map', data=[get(n, 'event_map') for n in events])
        table.addData(title='Pose', data=[get(n, 'pose') for n in events])

    @staticmethod
    def _f(value):
        try:
            return f'{float(value):.2f}'
        except (TypeError, ValueError):
            return value
