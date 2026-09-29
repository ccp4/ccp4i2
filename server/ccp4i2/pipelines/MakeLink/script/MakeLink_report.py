from ccp4i2.report.CCP4ReportParser import Report


class MakeLink_report(Report):
    TASKNAME = 'MakeLink'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv(style="clear:both;")
        self.addDefaultReport(self)
        clearingDiv = self.addDiv(style="clear:both;")

    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        for AcedrgLinkNode in self.xmlnode.findall("."):
            try:
                cycleNode = AcedrgLinkNode.findall("Cycle")[0]
            except:
                print("Missing cycle")
            try:
                logTextNode = AcedrgLinkNode.findall("LogText")[0]
            except:
                print("Missing logText")
            try:
                newFold = parent.addFold(label="Log text for iteration "+cycleNode.text, initiallyOpen=True)
            except:
                print("Unable to make fold")
            try:
                newFold.addPre(text = logTextNode.text)
            except:
                print("Unable to add Pre")
