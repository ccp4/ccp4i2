from ......report.CCP4ReportParser import Report


class pisa_analyse_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'pisa_analyse'
    USEPROGRAMXML = False

    def __init__(self, *args, **kws):
        Report.__init__(self, *args, **kws)
        self.addText(text='There is no meaningful report from the first run of Pisa - look at the second run')
