from ......report.CCP4ReportParser import Report


class phaser_EP_AUTO_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_EP_AUTO'
    RUNNING = True
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        self.addDiv(style='clear:both;')
        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.drawContent(jobStatus, self)
        
    def drawContent(self, jobStatus=None, parent=None):
        if parent is None: parent = self
        self.drawWarnings(parent=parent)
        atomSetNodes = self.xmlnode.findall('AtomSet')
        for atomSetNode in atomSetNodes:
            self.addText(text=atomSetNode.text,style="font-size:125%;")
            solNodes = atomSetNode.findall('AtomList')
            if len(solNodes)>0:
                self.addPre(text=solNodes[0].text)
            rfacNodes = atomSetNode.findall('R-factor')
            if len(rfacNodes)>0:
                self.addPre(text=rfacNodes[0].text)
            llNodes = atomSetNode.findall('Log-Likelihood')
            if len(llNodes)>0:
                self.addPre(text=llNodes[0].text)
            fomNodes = atomSetNode.findall('FiguresOfMerit')
            if len(fomNodes)>0:
                self.addPre(text=fomNodes[0].text)
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
        reportFold = parent.addFold(label='Plots from PHASER output',initiallyOpen=True)
        graphTableList = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table')

        gallery = reportFold.addObjectGallery(height='270px',contentWidth='500px',tableWidth='300px',style='float:left;width:820px;')
        for iGraph, graphTableNode in enumerate(graphTableList):
            graph = gallery.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title"), internalId='SmartiePlot'+str(iGraph), outputXml=False, label=graphTableNode.get("title"),style="width:500px;" )
            graph = graph.addPimpleData(xmlnode=graphTableNode)
        parent.addDiv(style='clear:both;')
        return
        
        graphgroup = reportFold.addFlotGraphGroup(style="width:500px;  height:300px;",outputXml=False, internalId='SmartiePlot')
        for iGraph, graphTableNode in enumerate(graphTableList):
            graph = graphgroup.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title"),internalId='SmartiePlot'+str(iGraph), outputXml=False )
            graph = graph.addPimpleData(xmlnode=graphTableNode)
        parent.addDiv(style='clear:both;')
        return



