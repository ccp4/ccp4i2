import math

from ccp4i2.report import Report


class ShelxCD_report(Report):
    TASKNAME = 'ShelxCD'
    RUNNING = True
    SEPARATEDATA=True

    def __init__(self, *args, **kws):
        Report.__init__(self, *args, **kws)
        self.outputXml = False
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None: parent = self
        if len(self.xmlnode.findall('.//Shelxd'))==0:
            self.shelXCReport(parent, initiallyOpen=True )
        else:
            self.shelXCReport(parent, initiallyOpen=False )
            self.shelXDReport(parent, initiallyOpen=True)

    def shelXCReport(self,parent=None,initiallyOpen = False):
        if parent is None: parent = self

        tagsAndKeys = {'N(data)':'NData','<I/sig>':'IOverSig','%Complete':'Completeness','<d"/sig>':'AnomalousSignal'}
        keysAndTags = {}
        for tag in tagsAndKeys: keysAndTags[tagsAndKeys[tag]] = tag

        shelxcFolder = parent.addFold(label='Shelxc data analysis', initiallyOpen = initiallyOpen)
        shelxcFolder.addText(text='Analysis of the different datasets provided')
        datasetNodes = self.xmlnode.findall('.//Dataset')
        if len(datasetNodes) > 0:
            galleryObject = shelxcFolder.addObjectGallery(style='width:100%;')
        for iDataset, datasetNode in enumerate(datasetNodes):
            if len(datasetNode.findall('DataAnalysis'))>0:
                try: title=datasetNode.findall('Name')[0].text
                except: title = 'Dataset'
                datasetDiv = galleryObject.addDiv(label=title)
                leftDiv, rightDiv = datasetDiv.addTwoColumnLayout(left_span=4, right_span=8, spacing=2)
                table = leftDiv.addTable(transpose=True,xmlnode=datasetNode,internalId='ShelxCDataset_Table_'+str(iDataset), outputXml=False)
                table.addData(title='<|E^2-1|>',select='eSquaredMinus1')
                table.addData(title='Refl. read',select='ReflectionsRead')
                table.addData(title='Unique refl.',select='UniqueReflections')
                table.addData(title='Highest res.',select='HighestResolution')
                graph = rightDiv.addFlotGraph(xmlnode=datasetNode, select = 'DataAnalysis/Bin',title=title,internalId='ShelxCDataset_Graph_'+str(iDataset), outputXml=False)
                resNodes = datasetNode.findall('DataAnalysis/Bin/HighRes')
                resData = [1./math.pow(float(resNode.text),2.0) for resNode in resNodes]
                graph.addData(title='Resolution',data=resData)
                iKey = 2
                for key in keysAndTags:
                    graph.addData(select = key, title=keysAndTags[key])
                    plot = graph.addPlotObject()
                    plot.append('title',keysAndTags[key])
                    plot.append('plottype','xy')
                    plot.append('xintegral','false')
                    plot.append('xscale','oneoversqrt')
                    plotLine = plot.append('plotline',xcol=1,ycol=iKey)
                    iKey += 1

    def shelXDReport(self, parent=None, initiallyOpen=False):
        if parent is None: parent = self
        shelXDNode = self.xmlnode.find('.//Shelxd')
        if shelXDNode is not None:
            shelxdFold = parent.addFold(label='ShelxD solutions',initiallyOpen=initiallyOpen)
            shelxdFold.addText(text='Quality metrics of different solutions')
            leftDiv, rightDiv = shelxdFold.addTwoColumnLayout(left_span=6, right_span=6, spacing=2)
            graph = leftDiv.addFlotGraph(xmlnode=shelXDNode,select='Try',title='Figures of merit',outputXml=False,internalId="FiguresOfMerit")
            plot = graph.addPlotObject()
            plot.append('title','Scores for each solution')
            iProperty = 1
            for property in ['Number','CCAll','CCWeak','CFOM','PATFOM']:
                graph.addData(select=property, title=property)
                if iProperty > 1:
                    plot.append('plottype','xy')
                    plot.append('xintegral','true')
                    plotLine = plot.append('plotline',xcol=1,ycol=iProperty)
                iProperty += 1
            graph = rightDiv.addFlotGraph(xmlnode=shelXDNode,select='Try',title='CCall vs CCweak',outputXml=False,internalId="CCAllVsCCWeak")
            plot = graph.addPlotObject()
            plot.append('title','CCall vs CCweak')
            iProperty = 1
            for property in ['CCAll','CCWeak']:
                graph.addData(select=property, title=property)
                if iProperty > 1:
                    plot.append('title',property)
                    plot.append('plottype','xy')
                    plot.append('xintegral','true')
                    plotLine = plot.append('plotline',xcol=1,ycol=iProperty)
                    plotLine.append('linestyle','.')
                iProperty += 1

            lstNodes = shelXDNode.findall('LstText')
            if len(lstNodes) > 0:
                lstFold = shelxdFold.addFold(label='Text of result_fa.lst file', initiallyOpen=False)
                for lstNode in lstNodes:
                    lstFold.addPre(outputXml=self.outputXml, internalId="ShelxDLstText", text = lstNode.text)


    def shelXEReport(self, parent=None, initiallyOpen=False, idRoot=""):
        if parent is None: parent = self
        globalTraceNodes = self.xmlnode.findall('Shelxe/GlobalAutotracingCycle')
        if len(globalTraceNodes) >0:
            shelxeFold = parent.addFold(label='ShelxE progress',initiallyOpen=initiallyOpen)
            shelxeFold.addText(text='Progress of density modification for each cycle of autobuild')
            galleryObject = shelxeFold.addObjectGallery(style='width:100%;')
            for globalTraceNode in reversed(globalTraceNodes):
                cycleNumber = globalTraceNode.findall('Number')[0].text
                cycleDiv = galleryObject.addDiv(label=cycleNumber)
                summaryText = ''
                fomNodes = globalTraceNode.findall('OverallFOMs/FOM')
                if len(fomNodes) > 0:
                    summaryText += ('FOM: ' + fomNodes[-1].text)
                ccNodes = globalTraceNode.findall('OverallFOMs/PseudoFreeCC')
                if len(ccNodes) > 0:
                    summaryText += (', Pseudo Free-CC: ' + ccNodes[-1].text)
                contrastNodes = globalTraceNode.findall('Cycle/Contrast')
                if len(contrastNodes) > 0:
                    summaryText += (', Contrast: ' + contrastNodes[-1].text)
                connectNodes = globalTraceNode.findall('Cycle/Connect')
                if len(connectNodes) > 0:
                    summaryText += (', Connect: ' + connectNodes[-1].text)
                cycleDiv.append('<span style="font-size:110%;">'+summaryText+'</span>')
                if len(connectNodes) > 0 or len(contrastNodes) > 0:
                    graph = cycleDiv.addFlotGraph(xmlnode=globalTraceNode,select='Cycle',outputXml=False, title='ShelXE', internalId=idRoot+'ShelXEStats_BuildCycle_'+cycleNumber)
                    plot = graph.addPlotObject()
                    iProperty = 1
                    for property in ['Number','Wt','Contrast','Connect']:
                        graph.addData(select=property, title=property)
                        if iProperty > 1:
                            plot.append('title',property)
                            plot.append('plottype','xy')
                            plot.append('xintegral','true')
                            plotLine = plot.append('plotline',xcol=1,ycol=iProperty)
                        iProperty += 1
            lstNodes = self.xmlnode.findall('Shelxe/LstText')
            if len(lstNodes) > 0:
                lstFold = shelxeFold.addFold(label='Text of ".lst" file', initiallyOpen=False)
                for lstNode in lstNodes:
                    lstFold.addPre(outputXml=self.outputXml, internalId=idRoot+"ShelxELstText", text = lstNode.text)
