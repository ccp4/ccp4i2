from ccp4i2.report.CCP4ReportParser import Report

class ZZPluginNameZZ_report(Report):
    TASKNAME = 'ZZPluginNameZZ'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv(style="clear:both;")
        self.addDefaultReport(self)
        
    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        if len(self.xmlnode.findall("LogText")) > 0:
            newFold = parent.addFold(label="Log text", initiallyOpen=True)
            newFold.addPre(text = self.xmlnode.findall("LogText")[0].text)
