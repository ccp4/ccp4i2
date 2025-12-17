from ccp4i2.report.CCP4ReportParser import Report
import sys

class phaser_EP_LLG_report(Report):
    TASKNAME = 'phaser_EP_LLG'
    RUNNING = True
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        self.addDiv(style='clear:both;')
        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.drawContent(jobStatus, self)
        
    def drawContent(self, jobStatus=None, parent=None):
        if parent is None: parent = self
        self.drawWarnings(parent=parent)
        self.addSmartieGraphs(parent=parent)

    def drawWarnings(self, parent=None):
        if parent is None: parent = self
        warnings = self.xmlnode.findall('.//PhaserWarnings/Warning')
        warningsFolder = parent.addFold(label='PHASER warnings', initiallyOpen=True)
        if len(warnings)>0:
            for warning in warnings:
                warningsFolder.addText(text=warning.text,style='color:red;')
                warningsFolder.append('<br/>')
        else:
            warningsFolder.addText(text='No warnings from PHASER')
        return

    def addSmartieGraphs(self, parent=None):
        if parent is None: parent=self
        reportFold = parent.addFold(label='Plots from PHASER output',initiallyOpen=False)
        graphTableList = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table')
        graphgroup = reportFold.addFlotGraphGroup(style="width:500px;  height:300px;")
        for graphTableNode in graphTableList:
            graph = graphgroup.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title") )
            graph = graph.addPimpleData(xmlnode=graphTableNode)
        parent.addDiv(style='clear:both;')



