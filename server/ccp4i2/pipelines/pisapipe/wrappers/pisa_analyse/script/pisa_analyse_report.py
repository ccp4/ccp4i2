from ccp4i2.report import Report


class pisa_analyse_report(Report):
    TASKNAME = 'pisa_analyse'
    USEPROGRAMXML = False

    def __init__(self,*args,**kws):
        Report.__init__(self, *args, **kws)
        self.addText(text='There is no meaningful report from the first run of Pisa - look at the second run')
