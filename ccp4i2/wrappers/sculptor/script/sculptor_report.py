from ....report.CCP4ReportParser import Report


class sculptor_report(Report):
    TASKNAME = 'sculptor'
    RUNNING = False
    USEPROGRAMXML = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        self.addText(text='Finished: ')
        pdbsWrittenPath = './/sculptor/number_output_files'
        pdbsWrittenStringS = xmlnode.findall(pdbsWrittenPath)
        if len(pdbsWrittenStringS) > 0:
            pdbsWrittenString = xmlnode.findall(pdbsWrittenPath)[0].text
            self.addText(text=pdbsWrittenString+' edited model(s) created')
