from ....report.CCP4ReportParser import Report


class chainsaw_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'chainsaw'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus, **kw)
        
        if jobStatus is not None and jobStatus.lower() == 'nooutput': return
        
        self.defaultReport(parent=self)
    
    def defaultReport(self, parent= None):
        if parent is None: parent= self
        parent.addDiv(style='clear:both;')

        for outputNode in self.xmlnode.findall('.//outputData'):
            newTable = parent.addTable(xmlnode=self.xmlnode)
            newTable.addData(title="Deleted", select=".//outputData/NO_DELETED")
            newTable.addData(title="Conserved",   select=".//outputData/NO_CONSERVED")
            newTable.addData(title="Mutated",   select=".//outputData/NO_MUTATED")
