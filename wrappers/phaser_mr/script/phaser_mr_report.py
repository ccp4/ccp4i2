from report.CCP4ReportParser import *
import sys
import xml.etree.ElementTree as etree
import base64

class pimpleGraph():
    def __init__(self,title=None, separator =' ', xmlnode=None):
        import copy
        if xmlnode is not None:
            self.xmlnode = copy.deepcopy(xmlnode)
            if title is not None: self.xmlnode.set('title',title)
            self.dataNode = self.xmlnode.findall('./data')[0]
            self.data = [dataLine.split(' ') for dataLine in self.dataNode.text.split('\n')]
            #The above will append a row of data with a single empty element if the last row in the XML terminates with a newline
            if len(self.data)>0 and len(self.data[-1]) == 1 and self.data[-1][0] == '':
                self.data = self.data[0:-1]
            self.headersNode = self.xmlnode.findall('./headers')[0]
            self.separator = self.headersNode.get('separator')
            self.headers = self.headersNode.text.strip().split(self.separator)
            return
        self.xmlnode = etree.Element('CCP4Table')
        self.xmlnode.set('title',title)
        self.dataNode = etree.SubElement(self.xmlnode,'data')
        self.data = [[]]
        self.dataNode.text = ''
        self.separator = separator
        self.headersNode = etree.SubElement(self.xmlnode,'headers',separator=self.separator)
        self.headersNode.text = ''
        self.headers = []
        return
    
    def dataAsText(self):
        import string
        textLines = [string.join(dataRow,' ') for dataRow in self.data]
        text = string.join(textLines,'\n')
        return text
    
    def nRows(self):
        return len(self.data)
    
    def nColumns(self):
        if self.nRows() == 0: return 0
        return len(self.data[0])
    
    def dataRow(self,row):
        return self.data[row]
    
    def appendPimpleGraph(self, otherPimpleGraph = None):
        if otherPimpleGraph is None: return
        import string
        import copy
        self.headers += otherPimpleGraph.headers
        self.headersNode.text = string.join(self.headers,self.separator)
        originalNRows = self.nRows()
        originalNColumns = self.nColumns()
        finalNRows = max(originalNRows, otherPimpleGraph.nRows())
        for row in range(finalNRows):
            if row>=originalNRows:
                self.data.append([])
                self.data[-1] += (['-' for i in range(originalNColumns)])
            if row>=otherPimpleGraph.nRows():
                self.data[row] += (['-' for i in range(otherPimpleGraph.nColumns())])
            else:
                self.data[row] += otherPimpleGraph.dataRow(row)
        self.dataNode.text = self.dataAsText()
        for plot in otherPimpleGraph.xmlnode.findall('./plot'):
            copiedPlot = copy.deepcopy(plot)
            for plotline in copiedPlot.findall('plotline'):
                intXcol = int(plotline.get('xcol'))
                plotline.set('xcol',str(intXcol+originalNColumns))
                intYcol = int(plotline.get('ycol'))
                plotline.set('ycol',str(intYcol+originalNColumns))
            self.xmlnode.append(copiedPlot)
        return

    def columnWithHeader(self,header):
        if not header in self.headers: return None
        iColumn = self.headers.index(header)
        return [dataRow[iColumn] for dataRow in self.data]

class phaser_mr_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_mr'
    RUNNING = True
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        if jobStatus.lower() == 'nooutput': return

        solutionNodes = xmlnode.findall(".//Solution")
        if len(solutionNodes)>0:
            self.append('<span style = "font-size:150%;"> Congratulations: solutions have been found</span>')
            solutionsFold = self.addFold(label='Solutions found', initiallyOpen=True,brief='Solutions')
            sideBySideDiv = solutionsFold.addDiv()
            lhs = sideBySideDiv.addDiv(style="float:left;width:250px;")
            self.addSolutionsTable(parent = lhs)
            rhs = sideBySideDiv.addDiv(style="float:left;")
            self.addSolutionsGraph(parent = rhs)
            clearingDiv = self.addDiv(style="clear:both;")
        elif jobStatus.lower() != 'running':
            self.append('<span style = "font-size:150%;"> Commiserations: no solutions have been found</span>')
            self.append('<span style = "font-size:100%;"> The last summary segment from Phaser follows:</span>')
            summaryNodes = self.xmlnode.findall('.//CCP4Summary')
            if len(summaryNodes)>1: self.addPre(text = base64.b64decode(summaryNodes[-2].text))
            if len(summaryNodes)>0: self.addPre(text = base64.b64decode(summaryNodes[-1].text))
        
        cellContentAnalysisFold = self.addFold(label = 'Cell Content analysis', initiallyOpen = False,brief='Cell')
        sideBySideDiv = cellContentAnalysisFold.addDiv()
        lhs = sideBySideDiv.addDiv(style="float:left;width:250px;")
        self.addAnalysisOfProblem(parent=lhs)
        rhs = sideBySideDiv.addDiv(style="float:left;")
        self.addPhaserCompositionAnalysis(parent = rhs)
        clearingDiv = self.addDiv(style="clear:both;")

        twinningAnalysisFold = self.addFold(label = 'Twinning analysis', initiallyOpen = False,brief='Twinning')
        self.addTwinDetection(parent = twinningAnalysisFold)

        searchResultsFold = self.addFold(label='Search results', initiallyOpen=True,brief='Search')
        self.addSearchResults(parent = searchResultsFold)

        summariesFold = self.addFold(label='Summaries output by Phaser', initiallyOpen=False,brief='Summaries')
        self.addSummaries(parent = summariesFold)

    def addAnalysisOfProblem(self, parent=None):
        if parent is None: parent = self
        parent.append('<span> Phaser was told that the asymmetric unit contains %s component(s) of %s type(s).  The neighbouring graph analyses the prior probability of how many copies of the composition provided to Phaser are likely to be in the asymmetric unit</span>'%(self.xmlnode.findall(".//Target/TotalComps")[0].text, self.xmlnode.findall(".//Target/CompTypes")[0].text))
    

    def addPhaserCompositionAnalysis(self,parent=None):
        if parent is None: parent = self
        cellContentAnalysisTableNode = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table[@title="Cell Content Analysis"]')[0].text
        if cellContentAnalysisTableNode is not None:
            cellContentAnalysisTableGraph = parent.addFlotGraph( xmlnode=cellContentAnalysisTableNode, title=cellContentAnalysisTableNode.get("title") , style = "width:300px; height:200px;")
            cellContentAnalysisTableData = cellContentAnalysisTableGraph.addPimpleData(xmlnode=cellContentAnalysisTableNode)

    def addSolutionsTable(self, parent = None):
        if parent is None: parent = self
        solutionNodes = self.xmlnode.findall('//Solutions/Solution')
        for solutionNode in solutionNodes:
            parent.append('<span>Solution %s has %s components located</span>'%(solutionNode.findall('.//ISOL')[0].text,solutionNode.findall('.//NCOMPONENTS')[0].text))
            table = parent.addTable(xmlnode=solutionNode, select='Component')
            for header, selector in [("Rot Func Z-score","RFZ"),("Trans Func Z-score","TFZ"),("Packing clashes","PAK"),("Log likelihood gain","LLG")]:
                table.addData(title = header, select = selector)

    def addSolutionsGraph(self, parent = None):
        if parent is None: parent = self
        solutionNodes = self.xmlnode.findall('.//Solutions/Solution')
        graph = parent.addFlotGraph(xmlnode=self.xmlnode, title = 'Solutions from phaser', select = './/Solutions/Solution',style="width:300px; height:200px;")
        graph.addData(title='Solution', select='ISOL')
        plot = graph.addPlotObject()
        plot.append('xintegral','true')
        plot.append('showlegend','false')
        for number, header, selector, colour, rightAxis  in [(2, "Components_found","NCOMPONENTS",'blue','false'),(3, "Overall_LLG","overallLLG",'green','true')]:
            graph.addData(title = header, select = selector)
            plot.append('yrange', rightaxis=rightAxis)
            plotLine = plot.append('plotline',xcol=1,ycol=number,rightaxis=rightAxis,colour=colour)

    def addSearchResults(self, parent = None):
        import string
        if parent is None: parent = self
        searchNodes = self.xmlnode.findall('.//Search')
        for iSearchNode in range(len(searchNodes)):
            searchNode = searchNodes[len(searchNodes)-(1+iSearchNode)]
            searchDiv = parent.addDiv(style='float:left; border: 1px solid black;width:1024px;')
            searchSummaryNodes = searchNode.findall('.//CCP4Summary')
            if len(searchSummaryNodes)> 0:
                #Unclear why I now have to put a div in to contain things....
                summaryDiv = searchDiv.addDiv(style="width:250px; border-width: 0px; border-color: black; border-style:solid; float:left;")
                summaryDiv.addPre(text=searchSummaryNodes[0].text)
            
            graphNodes = searchNode.findall('.//CCP4Table')
            if len(graphNodes) > 0:
                
                Labels = []
                Zs = []
                LLGs = []
                Rs = []
                
                #graphGroup = searchDiv.addFlotGraphGroup(style='width:300px; height:200px; float:left;')
                graphGallery = searchDiv.addObjectGallery(tableWidth='200px',contentWidth='300px',height='250px', style='float:left;')
                firstGraph = True
                for graphNode in graphNodes:
                    #seed a pimpleGraph frmo this...easier to interrogate column values.
                    asPimpleGraph = pimpleGraph(xmlnode=graphNode)
                    graph = graphGallery.addFlotGraph(xmlnode = graphNode, title=graphNode.get('title'),style='width:290px;height:240px;border:0px solid green;',initiallyDrawn=firstGraph, withLaunch = False)
                    firstGraph=False
                    graph.addPimpleData(xmlnode=graphNode)
                    
                    Labels.append(graphNode.get('title').split()[0])
                    
                    Z = '-'
                    zCol = asPimpleGraph.columnWithHeader('Z-Score')
                    if zCol is not None and len(zCol) > 0: Z=zCol[0]
                    Zs.append(Z)
                    
                    LLG = '-'
                    llgCol = asPimpleGraph.columnWithHeader('LLG')
                    if llgCol is not None and len(llgCol) > 0: LLG=llgCol[0]
                    else:
                        llgCol = asPimpleGraph.columnWithHeader('final-LLG')
                        if llgCol is not None and len(llgCol) > 0: LLG=llgCol[0]
                    LLGs.append(LLG)
                    
                    R = '-'
                    rCol = asPimpleGraph.columnWithHeader('final-R')
                    if rCol is not None and len(rCol) > 0:R=rCol[0]
                    Rs.append(R)
                
                summaryTable = summaryDiv.addTable(style='clear:both;')
                summaryTable.addData(data=Labels,title='Phase')
                summaryTable.addData(data=Zs,title='Z-Score')
                summaryTable.addData(data=LLGs,title='LLG')
                summaryTable.addData(data=Rs,title='R-factor')
    
        clearingDiv = parent.addDiv(style="clear:both;")
            
    def addTwinDetection(self, parent = None):
        if parent is None: parent = self
        summaryNodes = self.xmlnode.findall('.//CCP4Summary')
        for summaryNode in summaryNodes:
            if summaryNode.text is not None and 'tNCS/Twin Detection Table' in summaryNode.text:
                parent.addPre(text = summaryNode.text, style="clear:both;")
                intensityDistributionNodes = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table[@title="Intensity distribution for Data"]')
                if len(intensityDistributionNodes) > 0:
                    graphsDiv = parent.addDiv(style="clear:both;")
                    for intensityDistributionNode in intensityDistributionNodes:
                        intensityDistributionGraph = graphsDiv.addFlotGraph( xmlnode=intensityDistributionNode, title=intensityDistributionNode.get("title") , style = "width:300px; height:200px; float:left;")
                        intensityDistributionGraph.addPimpleData(xmlnode=intensityDistributionNode)
        clearingDiv = parent.addDiv(style="clear:both;")


    def addSummaries(self, parent = None):
        if parent is None: parent = self
        #Unclear why I now have to put a div in to contain things....
        summaryDiv = parent.addDiv()
        summaryNodes = self.xmlnode.findall('.//CCP4Summary')
        for summaryNode in summaryNodes:
            summaryDiv.addPre(text = summaryNode.text)

