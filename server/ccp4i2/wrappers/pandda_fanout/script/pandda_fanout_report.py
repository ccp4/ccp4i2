from ccp4i2.report import Report


class pandda_fanout_report(Report):
    TASKNAME = 'pandda_fanout'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style="clear:both;")
        if xmlnode is None:
            return
        text = lambda path, default='': (
            xmlnode.findtext(path) if xmlnode.findtext(path) is not None else default)
        if text('dry_run') == 'True':
            self.addText(text='Preview only: nothing was created.', style='font-weight:bold;')
        if text('incomplete') == 'True':
            self.addText(text=('The run did not write its events table: a partial tree. Receipts were '
                               'taken for what is there and told the run did not finish.'),
                         style='color:#e65100;')
        self.addText(text=('This job records what the fan-out created at the time. Each receipt is '
                           'its own job in its own project, and that is where its state lives.'))

        summary = self.addFold(label='Summary', brief='Summary', initiallyOpen=True)
        table = summary.addTable(transpose=True)
        table.addData(title='Output tree', data=[text('tree')])
        table.addData(title='Run job', data=[text('run_job_uuid') or '(none: external tree)'])
        for action, title in (('created', 'Receipts created'), ('skipped', 'Already had a receipt'),
                              ('absent', 'Not in the tree'), ('failed', 'Failed'),
                              ('no_project', 'Project not in this database')):
            table.addData(title=title, data=[text(f'n_{action}', '0')])

        rows = xmlnode.findall('datasets/dataset')
        if rows:
            fold = self.addFold(label='Datasets', brief='Datasets', initiallyOpen=True)
            table = fold.addTable()
            table.addData(title='xtal', data=[r.get('xtal') for r in rows])
            table.addData(title='Label', data=[r.get('label') for r in rows])
            table.addData(title='Project', data=[r.get('project') or '-' for r in rows])
            table.addData(title='Outcome', data=[r.get('action') for r in rows])
            table.addData(title='Detail', data=[r.text or '' for r in rows])
