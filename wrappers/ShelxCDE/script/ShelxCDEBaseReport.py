from report.CCP4ReportParser import *
import sys
from lxml import etree
import math

class ShelxCDEBaseReport(Report):
    
    def __init__(self, *args, **kws):
        Report.__init__(self, *args, **kws)
        #self.outputXml = self.jobStatus is not None and self.jobStatus.lower() == 'running'
        self.outputXml = False
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.addDiv(style='clear:both;')
        self.defaultReport()

    def shelXCReport(self,parent=None,initiallyOpen = False):
        if parent is None: parent = self
        
        tagsAndKeys = {'N(data)':'NData','<I/sig>':'IOverSig','%Complete':'Completeness','<d"/sig>':'AnomalousSignal'}
        keysAndTags = {}
        for tag in tagsAndKeys: keysAndTags[tagsAndKeys[tag]] = tag
        
        shelxcFolder = parent.addFold(label='Shelxc data analysis', initiallyOpen = initiallyOpen)
        shelxcFolder.addText(text='Analysis of the different datasets provided')
        datasetNodes = self.xmlnode.xpath('//Dataset')
        if len(datasetNodes) > 0:
            galleryObject = shelxcFolder.addObjectGallery(contentWidth='600px',height='214px',tableWidth='50px',style='border:1px solid black;float:left;width:660px;height:224px;')
        for iDataset, datasetNode in enumerate(datasetNodes):
            if len(datasetNode.xpath('DataAnalysis'))>0:
                try: title=datasetNode.xpath('Name')[0].text
                except: title = 'Dataset'
                datasetDiv = galleryObject.addDiv(style='border:0px solid white;width:590px;height:210px;float:left;',label=title)
                tableDiv = datasetDiv.addDiv(style='border:0px solid white;width:200px;height:210px;float:left;')
                table = tableDiv.addTable(transpose=True,xmlnode=datasetNode,internalId='ShelxCDataset_Table_'+str(iDataset), outputXml=False)
                table.addData(title='<|E^2-1|>',select='eSquaredMinus1')
                table.addData(title='Refl. read',select='ReflectionsRead')
                table.addData(title='Unique refl.',select='UniqueReflections')
                table.addData(title='Highest res.',select='HighestResolution')
                graph = datasetDiv.addFlotGraph(xmlnode=datasetNode, select = 'DataAnalysis/Bin',style='width:300px;height:200px;border:0px solid white;float:left;',title=title,internalId='ShelxCDataset_Graph_'+str(iDataset), outputXml=False)
                resNodes = datasetNode.xpath('DataAnalysis/Bin/HighRes')
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
        parent.addDiv(style='clear:both;')

    def shelXDReport(self, parent=None, initiallyOpen=False):
        if parent is None: parent = self
        shelXDNode = self.xmlnode.xpath0('//Shelxd')
        if shelXDNode is not None:
            shelxdFold = parent.addFold(label='ShelxD solutions',initiallyOpen=initiallyOpen)
            shelxdFold.addText(text='Quality metrics of different solutions')
            graphsDiv = shelxdFold.addDiv(style='width:620px; height:350px; border:1px solid black;float:left;')
            graph = graphsDiv.addFlotGraph(xmlnode=shelXDNode,select='Try',style='width:300px;height:300px;border:0px solid white;float:left;',title='Figures of merit',outputXml=False,internalId="FiguresOfMerit")
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
            graph = graphsDiv.addFlotGraph(xmlnode=shelXDNode,select='Try',style='width:300px;height:300px;border:0px solid white;float:left;',title='CCall vs CCweak',outputXml=False,internalId="CCAllVsCCWeak")
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
            shelxdFold.addDiv(style='clear:both;')
            
            lstNodes = shelXDNode.xpath('LstText')
            if len(lstNodes) > 0:
                lstFold = shelxdFold.addFold(label='Text of result_fa.lst file', initiallyOpen=False)
                for lstNode in lstNodes:
                    lstFold.addPre(outputXml=self.outputXml, internalId="ShelxDLstText", text = lstNode.text)
            

    def shelXEReport(self, parent=None, initiallyOpen=False, idRoot=""):
        if parent is None: parent = self
        globalTraceNodes = self.xmlnode.xpath('Shelxe/GlobalAutotracingCycle')
        if len(globalTraceNodes) >0:
            shelxeFold = parent.addFold(label='ShelxE progress',initiallyOpen=initiallyOpen)
            shelxeFold.addText(text='Progress of density modification for each cycle of autobuild')
            galleryObject = shelxeFold.addObjectGallery(contentWidth='600px',height='230px',tableWidth='50px',style='border:1px solid black;float:left;width:660px;height:230px;')
            for globalTraceNode in reversed(globalTraceNodes):
                cycleNumber = globalTraceNode.xpath('Number')[0].text
                cycleDiv = galleryObject.addDiv(style='border:1px solid green;height:230px; width:600px;',label=cycleNumber)
                summaryText = ''
                fomNodes = globalTraceNode.xpath('OverallFOMs/FOM')
                if len(fomNodes) > 0:
                    summaryText += ('FOM: ' + fomNodes[-1].text)
                ccNodes = globalTraceNode.xpath('OverallFOMs/PseudoFreeCC')
                if len(ccNodes) > 0:
                    summaryText += (', Pseudo Free-CC: ' + ccNodes[-1].text)
                contrastNodes = globalTraceNode.xpath('Cycle/Contrast')
                if len(contrastNodes) > 0:
                    summaryText += (', Contrast: ' + contrastNodes[-1].text)
                connectNodes = globalTraceNode.xpath('Cycle/Connect')
                if len(connectNodes) > 0:
                    summaryText += (', Connect: ' + connectNodes[-1].text)
                cycleDiv.append('<span style="font-size:110%;">'+summaryText+'</span>')
                if len(connectNodes) > 0 or len(contrastNodes) > 0:
                    graph = cycleDiv.addFlotGraph(xmlnode=globalTraceNode,select='Cycle',style='width:450px;height:200px;border:0px solid white;',outputXml=False, title='ShelXE', internalId=idRoot+'ShelXEStats_BuildCycle_'+cycleNumber)
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
            lstNodes = self.xmlnode.xpath('Shelxe/LstText')
            if len(lstNodes) > 0:
                shelxeFold.addDiv(style='clear:both;')
                lstFold = shelxeFold.addFold(label='Text of ".lst" file', initiallyOpen=False)
                for lstNode in lstNodes:
                    lstFold.addPre(outputXml=self.outputXml, internalId=idRoot+"ShelxELstText", text = lstNode.text)

