from __future__ import print_function
from report.CCP4ReportParser import *
import sys

class phaser_MR_FTF_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_MR_FTF'
    RUNNING = True
    SEPARATEDATA=True
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)
        self.outputXml = self.jobStatus is not None and self.jobStatus.lower() == 'running'
        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.drawContent(jobStatus, self)
        
    def drawContent(self, jobStatus=None, parent=None):
        self.addDiv(style='clear:both;')
        #print 'self.xmlnode',self.xmlnode
        if parent is None: parent=self
        #print jobStatus
        if jobStatus is not None and jobStatus.lower() == 'running':
            self.outputXml = True
            self.addProgressBar(parent=parent)
        
        self.drawAdvisories(parent=parent)
        self.drawWarnings(parent=parent)
        
        bestSolFold = parent.addFold(label='Elements and scores of current solution', brief='Current soln.', initiallyOpen=True)
        if len(self.xmlnode.findall('PhaserCurrentBestSolution')) > 0:
            self.addBestSolution(parent=bestSolFold)
        else:
            bestSolFold.addText(text='No solution found yet', style='color:red;',outputXML=self.outputXml, internalId='BestSolText')
        if len(self.xmlnode.findall('PhaserMrSolutions/Solutions')) > 0:
            compareSolutionsFold = parent.addFold(label='Comparison of solutions', brief='All solns', initiallyOpen=True)
            self.addResults(parent=compareSolutionsFold)
        analysisFold = parent.addFold(label='Analysis of composition and data', brief='Comp/data', initiallyOpen=False)
        searchPathFold = parent.addFold(label='Search strategy employed by PHASER', brief='Search tree', initiallyOpen=False)
        indentText = "...."
        addIndent = False
        removeIndent = False
        import re
        for iNode, summaryNode in enumerate(self.xmlnode.findall('Summary')):
            destination = searchPathFold
             #Make a result folder for any summmary blocks in the XML
            try:
                moduleLines = [line for line in summaryNode.text.split('\n') if "Phaser Module:" in line]
                summaryBlockTitle = moduleLines[0][19:-9].strip()
                
                if summaryBlockTitle == "MOLECULAR REPLACEMENT TRANSLATION FUNCTION":
                    try:
                        summaryBlockTitle += (': '+str(re.findall(r'^.*Top TFZ.*$',summaryNode.text,re.MULTILINE)[0]).strip())
                    except:
                        pass
                summaryBlockTitle = summaryBlockTitle[0:1]+summaryBlockTitle[1:].lower()
                if destination is not None:
                    phaserSummaryFold = destination.addFold(label=indentText + summaryBlockTitle, initiallyOpen=False)
                    phaserSummaryFold.addPre(text=summaryNode.text, outputXml=self.outputXml, internalId='OverallSummaryText'+str(iNode))
                try:
                    if addIndent: indentText += '....'
                    elif removeIndent: indentText = indentText[4:]
                except:
                    print('Overrun on manipulation of indent')
                addIndent = False
                removeIndent = False
            except:
                #Possible exceptions include there being no text in the summaryNode
                pass

        self.addSmartieGraphs(parent=parent)

    def drawAdvisories(self, parent=None):
        if parent is None: parent = self
        advisories = self.xmlnode.findall('.//PhaserAdvisories/Advisory')
        if len(advisories)>0:
            parent.append('<br/>')
            for advisory in advisories:
                parent.addText(text=advisory.text,style='color:orange;',outputXML=self.outputXml, internalId='AdvisoriesText')
                parent.append('<br/>')
        return

    def drawWarnings(self, parent=None):
        if parent is None: parent = self
        warnings = self.xmlnode.findall('.//PhaserWarnings/Warning')
        #warningsFolder = parent.addFold(label='PHASER warnings', initiallyOpen=True)
        if len(warnings)>0:
            parent.append('<br/>')
            for warning in warnings:
                parent.addText(text=warning.text,style='color:red;',outputXML=self.outputXml, internalId='WarningsText')
                parent.append('<br/>')
        return

    def addSearchProblem(self, parent=None):
        pass
    
    def addBestSolution(self, parent=None):
        if parent is None: parent = self
        solutionNode = self.xmlnode.findall('PhaserCurrentBestSolution/Solution')[0]
        
        summaryText = ""
        if len(solutionNode.findall('./spaceGroup')) > 0:
            spaceGroup = solutionNode.findall('./spaceGroup')[0].text
            summaryText += 'Current best solution has spacegroup %s'% spaceGroup
        if len(solutionNode.findall('./overallLLG')) > 0:
            summaryText +=' Refined overall LLG: %s'%(solutionNode.findall('./overallLLG')[0].text)
        if len(summaryText)>0: parent.addText(text=summaryText,outputXml=self.outputXml,internalId='BestSolutionSummary')
        
        table = parent.addTable(xmlnode=solutionNode, outputXml=self.outputXml, internalId='BestSolutionBreakdown')
        for header, selector in [("Ensemble name","Component/Name"),("Rot Func Z-score","Component/RFZ"),("Trans Func Z-score","Component/TFZ")]:
            if len (solutionNode.findall(selector)) > 0:
                table.addData(title = header, data=[component.text for component in solutionNode.findall(selector)],tip="Stored in xml column {0}".format(selector))
                #print header, [component.text for component in solutionNode.findall(selector)]

    def addResults(self, parent=None):
        if parent is None: parent=self
        resultsNodes = self.xmlnode.findall('PhaserMrSolutions/Solutions')
        if len(resultsNodes) > 0:
            dataLabelTuples = [('SPG','Space group'), ('TFZ','Trans. Z-score')]
            if len(resultsNodes[0].findall('Solution')) > 1:
                parent.addText(text='Multiple solutions found',outputXML=self.outputXml, internalId='MultipleSolText')
                parent.append('<br/>')
                graph = parent.addFlotGraph(title='Solutions',xmlnode=resultsNodes[0], select='Solution',style='height:200px; width:600px;',outputXml=self.outputXml, internalId='PhaserSolutions')
                plotObject = graph.addPlotObject()
                plotObject.append('title','Solutions')
                plotObject.append('plottype','xy')
                data=[i for i in range(len(resultsNodes[0].findall('Solution')))]
                graph.addData(title='Solution number'.replace(' ','_').replace('.','_'), data=data )
                iLine = 0
                for selector,label in dataLabelTuples:
                    if selector != 'SPG':
                        graph.addData(title = label.replace(' ','_').replace('.','_'), select = selector)
                        plotLine = plotObject.append('plotline',xcol=1,ycol=iLine+2)
                        iLine += 1
            elif len(resultsNodes[0].findall('Solution')) == 1:
                parent.addText(text='Unique solution found :-)',outputXML=self.outputXml, internalId='UniqueSolText')
            elif len(resultsNodes[0].findall('Solution')) == 0:
                parent.addText(text='No solutions found',style='color:red;',outputXML=self.outputXml, internalId='NoSolText')
                return
            table = parent.addTable(xmlnode=resultsNodes[0], select='Solution',style='width:600px;', outputXml=self.outputXml, internalId='ResultsTable')
            for selector,label in dataLabelTuples:
                if selector == 'SPG': table.addData(title = 'Space group', select = selector)
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
        parent.addProgress(style="width:500px;border:1px solid green;",value=value,max=max,internalId='PhaserProgress',outputXml=self.outputXml, label=label)
            
    def addSmartieGraphs(self, parent=None):
        if parent is None: parent=self
        reportFold = parent.addFold(label='Plots from PHASER output',initiallyOpen=True)
        graphTableList = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table')

        gallery = reportFold.addObjectGallery(height='270px',contentWidth='500px',tableWidth='300px',style='float:left;width:820px;')
        for iGraph, graphTableNode in enumerate(graphTableList):
            graph = gallery.addFlotGraph( xmlnode=graphTableNode, title=graphTableNode.get("title"), internalId='SmartiePlot'+str(iGraph), outputXml=self.outputXml, label=graphTableNode.get("title"),style="width:500px;" )
            graphDataNodes = graphTableNode.findall('data')
            wasTruncated = False
            for graphDataNode in graphDataNodes:
                if graphDataNode.text.count('\n') > 2000:
                    graphDataNode.text = '\n'.join(graphDataNode.text.split('\n')[0:1999])
                    wasTruncated = True
            if wasTruncated: graph.label = "TRUNCATED:"+ graph.label
            graph.addPimpleData(xmlnode=graphTableNode, usePlotly=False)
        parent.addDiv(style='clear:both;')
        return
        



