from ....report.CCP4ReportParser import Report


class ZZPipelineNameZZ_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ZZPipelineNameZZ'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style="clear:both;")
        self.addDefaultReport(self)
        
    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        for ZZFirstPluginNameZZNode in self.xmlnode.findall(".//ZZFirstPluginNameZZ"):
            try:
                cycleNode = ZZFirstPluginNameZZNode.findall("Cycle")[0]
            except:
                print("Missing cycle")
            try:
                newFold = parent.addFold(label="Log text for iteration "+cycleNode.text, initiallyOpen=True)
            except:
                print("Unable to make fold")
