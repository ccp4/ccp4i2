from report.CCP4ReportParser import *
import sys

class phaser_MR_RNP_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_MR_RNP'
    RUNNING = True
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.drawContent(jobStatus, self)
        
    def drawContent(self, jobStatus=None, parent=None):
        #print 'self.xmlnode',self.xmlnode
        if parent is None: parent=self
        if jobStatus is not None and jobStatus.lower() == 'running':
            self.addProgressBar(parent=parent)
        self.drawWarnings(parent=parent)
        
        bestSolFold = parent.addFold(label='Elements and scores of current solution',initiallyOpen=True)
        if len(self.xmlnode.findall('PhaserCurrentBestSolution')) > 0:
            self.addBestSolution(parent=bestSolFold)
        else:
            bestSolFold.addText(text='No solution found yet', style='color:red;')
        if len(self.xmlnode.findall('PhaserMrSolutions/Solutions')) > 0:
            compareSolutionsFold = parent.addFold(label='Comparison of solutions',initiallyOpen=True)
            self.addResults(parent=compareSolutionsFold)

        try: cellContentAnalysisTableNode = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table[@title="Cell Content Analysis"]')[0].text
        except: cellContentAnalysisTableNode = None
        try: dataStatisticsNode = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table[@title="Intensity distribution for Data"]')[0].text
        except: dataStatisticsNode = None
        if cellContentAnalysisTableNode is not None or dataStatisticsNode is not None:
            analysisFold = parent.addFold(label='Analysis of composition and data',initiallyOpen=False)
            if cellContentAnalysisTableNode is not None:
                cellContentAnalysisTableGraph = analysisFold.addGraph( xmlnode=cellContentAnalysisTableNode, title=cellContentAnalysisTableNode.get("title") , style = "width:300px; height:200px;")
                cellContentAnalysisTableData = cellContentAnalysisTableGraph.addPimpleData(xmlnode=cellContentAnalysisTableNode)
            if dataStatisticsNode is not None:
                dataStatisticsNodeGraph = analysisFold.addGraph( xmlnode=dataStatisticsNode, title=dataStatisticsNode.get("title") , style = "width:300px; height:200px;")
                dataStatisticsNodeData = dataStatisticsNodeGraph.addPimpleData(xmlnode=dataStatisticsNode)
    
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

    def addSearchProblem(self, parent=None):
        pass
    
    def addBestSolution(self, parent=None):
        if parent is None: parent = self
        solutionNode = self.xmlnode.findall('PhaserCurrentBestSolution/Solution')[0]
        if solutionNode is not None:
            if self.xmlnode.findall('PhaserCurrentBestSolution/spaceGroup')[0] is not None:
                parent.append('<span>Current best solution <br/> spacegroup %s </span>'%(self.xmlnode.findall('PhaserCurrentBestSolution/spaceGroup')[0].text))
            if len(solutionNode.findall('./overallLLG')) > 0:
                parent.append('<span>Refined overall LLG: %s </span>'%(solutionNode.findall('./overallLLG')[0].text))
            table = parent.addTable(xmlnode=solutionNode, select='Component')
            for header, selector in [("Ensemble name","Name"),("Rot Func Z-score","RFZ"),("Trans Func Z-score","TFZ"),("Packing clashes","PAK"),("Log likelihood gain","LLG")]:
                table.addData(title = header, select = selector)

    def addResults(self, parent=None):
        if parent is None: parent=self
        resultsNodes = self.xmlnode.findall('PhaserMrSolutions/Solutions')
        if len(resultsNodes) > 0:
            dataLabelTuples = [('HALL','Space group'), ('TFZ','Trans. Z-score'), ('TFZeq','Refined Trans. Z-score'), ('ORIG_LLG','Initial LLG'), ('LLG','Refined LLG'),  ('ORIG_R','Initial Rfactor'), ('R','Refined Rfactor'),('PAK','Clashes')]
            if len(resultsNodes[0].findall('Solution')) > 1:
                parent.addText(text='Multiple solutions found')
                parent.append('<br/>')
                graph = parent.addFlotGraph(title='Solutions',xmlnode=resultsNodes[0], select='Solution',style='height:200px; width:600px;')
                plotObject = graph.addPlotObject()
                plotObject.append('title','Solutions')
                plotObject.append('plottype','xy')
                plotObject.append('yrange', rightaxis='true')
                data=[i for i in range(len(resultsNodes[0].findall('Solution')))]
                graph.addData(title='Solution number'.replace(' ','_').replace('.','_'), data=data )
                iLine = 0
                for selector,label in dataLabelTuples:
                    if selector != 'HALL':
                        graph.addData(title = label.replace(' ','_').replace('.','_'), select = selector)
                        if selector == 'R' or selector == 'ORIG_R':
                            plotLine = plotObject.append('plotline',xcol=1,ycol=iLine+2,rightaxis='true')
                        else:
                            plotLine = plotObject.append('plotline',xcol=1,ycol=iLine+2)
                        iLine += 1
            elif len(resultsNodes[0].findall('Solution')) == 1:
                parent.addText(text='Unique solution found :-)')
            elif len(resultsNodes[0].findall('Solution')) == 0:
                parent.addText(text='No solutions found',style='color:red;')
                return
            table = parent.addTable(xmlnode=resultsNodes[0], select='Solution',style='width:600px;')
            for selector,label in dataLabelTuples:
                if selector == 'HALL': table.addData(title = label, select = selector)
                else: table.addData(title = label, select = selector, expr='"{0:.2f}".format(float(x))')

    def addProgressBar(self, parent=None):
        if parent is None: parent=self
        label, max, value  = ('No action in progress', 0, 0)
        currentActivityNodes = self.xmlnode.findall('CurrentActivity')
        if len(currentActivityNodes)>0:
            currentActivityNode = currentActivityNodes[0]
            label = currentActivityNode.findall('label')[0].text
            max = currentActivityNode.findall('max')[0].text
            value = currentActivityNode.findall('value')[0].text
        parent.addProgress(style="width:500px;border:1px solid green;",value=value,max=max,label=label)

    def addSmartieGraphs(self, parent=None):
        if parent is None: parent=self
        reportFold = parent.addFold(label='Plots from PHASER output',initiallyOpen=False)
        graphTableList = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table')
        graphgroup = reportFold.addFlotGraphGroup(style="width:500px;  height:300px;")
        for graphTableNode in graphTableList:
            graph = graphgroup.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title") )
            graph = graph.addPimpleData(xmlnode=graphTableNode)
        parent.addDiv(style='clear:both;')



