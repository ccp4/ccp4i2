from __future__ import print_function
from report.CCP4ReportParser import *
import sys

class phaser_MR_AUTO_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_MR_AUTO'
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
            self.outputXml = False
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
        '''
        cellContentAnalysisTableNode = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table[@title="Cell Content Analysis"]')[0].text
        dataStatisticsNode = self.xmlnode.findall('.//CCP4ApplicationOutput/CCP4Table[@title="Intensity distribution for Data"]')[0].text
        if cellContentAnalysisTableNode is not None or dataStatisticsNode is not None:
        
            if cellContentAnalysisTableNode is not None:
                cellContentAnalysisTableGraph = analysisFold.addFlotGraph( xmlnode=cellContentAnalysisTableNode, title=cellContentAnalysisTableNode.get("title") , style = "width:300px; height:200px;", outputXml=True, internalId='CellContentAnalysis')
                cellContentAnalysisTableData = cellContentAnalysisTableGraph.addPimpleData(xmlnode=cellContentAnalysisTableNode)
            if dataStatisticsNode is not None:
                dataStatisticsNodeGraph = analysisFold.addFlotGraph( xmlnode=dataStatisticsNode, title=dataStatisticsNode.get("title") , style = "width:300px; height:200px;", outputXml=True, internalId='DataAnalysis')
                dataStatisticsNodeData = dataStatisticsNodeGraph.addPimpleData(xmlnode=dataStatisticsNode)
        '''
    
        comFileFold = parent.addFold(label='COM file for this run', brief='COM file', initiallyOpen=False)
        
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
                if summaryBlockTitle == 'AUTOMATED MOLECULAR REPLACEMENT':
                    if 'Number of search ensembles' in summaryNode.text:
                        summaryBlockTitle = 'Search ensembles and spacegroups'
                        destination = analysisFold
                    elif 'Composition Table' in summaryNode.text:
                        summaryBlockTitle = 'Composition table'
                        destination = analysisFold
                    elif 'eLLG Values Computed for All Data' in summaryNode.text:
                        summaryBlockTitle = 'eLLG Values Computed for All Data'
                        destination = analysisFold
                    elif 'Search: Next component' in summaryNode.text:
                        summaryBlockTitle = 'Next component search'
                        summaryBlockTitle += [str(text) for text in re.findall(r'^.*Ensemble.*$',summaryNode.text,re.MULTILINE) if '*' in text][0]
                        addIndent = True
                    elif 'Search Order (next search *):' in summaryNode.text:
                        indentText = ''
                        phaserSummaryFold = destination.addFold(label=indentText + 'New search order', initiallyOpen=False)
                        indentText += '....'
                        summaryBlockTitle = 'Next component search:'
                        summaryBlockTitle += [str(text) for text in re.findall(r'^.*Ensemble.*$',summaryNode.text,re.MULTILINE) if '*' in text][0]
                        addIndent = True
                    elif 'Translation Function has peaks' in summaryNode.text:
                        summaryBlockTitle = 'Definite translation result'
                    elif 'Composition Table' in summaryNode.text:
                        destination = analysisFold
                        summaryBlockTitle = 'Composition table'
                    elif 'No definite solutions found' in summaryNode.text:
                        summaryBlockTitle = 'Failed search'
                        removeIndent = True
                    elif 'fixed components' in summaryNode.text:
                        summaryBlockTitle = 'Solutions to refine'
                        try:
                            summaryBlockTitle += (': '+str(re.findall(r'^.*Number of solutions.*$',summaryNode.text,re.MULTILINE)[0]).strip())
                        except:
                            pass
                    elif 'New Best LLG' in summaryNode.text:
                        summaryBlockTitle = ""
                        for bestLine in re.findall(r'^.*Best.*$',summaryNode.text,re.MULTILINE):
                            summaryBlockTitle += (bestLine.replace('Solution','Soln').replace('Search Component','Comp.').replace('Current is Best Soln','').replace('**','*') + ' ')
                    elif 'Definite solutions found' in summaryNode.text:
                        summaryBlockTitle = 'Definite solution'
                        removeIndent = True
                    elif 'End of possible ensembles' in summaryNode.text:
                        summaryBlockTitle = str(re.findall(r'^.*Search.*$',summaryNode.text,re.MULTILINE)[0]).strip()
                    elif 'Purge all but the amalgamated solutions' in summaryNode.text:
                        summaryBlockTitle = 'Purge smaller amalgamation groups'
                    elif 'ssful amalgamation' in summaryNode.text:
                        summaryBlockTitle = 'End of amalgamation attempt'
                        removeIndent = True
                else:
                    if "Translational non-crystallographic symmetry".upper() in summaryBlockTitle or \
                       "anisotropy correction".upper() in summaryBlockTitle or \
                       "cell content analysis".upper() in summaryBlockTitle or \
                       "Read data from mtz file" .upper() in summaryBlockTitle :
                        destination=analysisFold
                    elif 'PREPROCESSOR'.upper() in summaryBlockTitle:
                        comFileFold.addPre(text=summaryNode.text, outputXml=self.outputXml, internalId='ComFileText')
                        destination = None
                    if summaryBlockTitle == "MOLECULAR REPLACEMENT AMALGAMATION":
                        summaryBlockTitle = 'Begin amalgamation attempt'
                        addIndent = True
                    elif summaryBlockTitle == "MOLECULAR REPLACEMENT PACKING ANALYSIS":
                        try:
                            summaryBlockTitle += (': '+str(re.findall(r'^.*accepted of.*$',summaryNode.text,re.MULTILINE)[0]).strip())
                        except:
                            pass
                    elif summaryBlockTitle == "MOLECULAR REPLACEMENT TRANSLATION FUNCTION":
                        try:
                            summaryBlockTitle += (': '+str(re.findall(r'^.*Top TFZ.*$',summaryNode.text,re.MULTILINE)[0]).strip())
                        except:
                            pass
                    elif summaryBlockTitle == "MOLECULAR REPLACEMENT ROTATION FUNCTION":
                        try:
                            summaryBlockTitle += (': '+str(re.findall(r'^.*Top RF.*$',summaryNode.text,re.MULTILINE)[0]).strip())
                        except:
                            pass
                    elif summaryBlockTitle == "MOLECULAR REPLACEMENT REFINEMENT AND PHASING":
                        try:
                            summaryBlockTitle += (': '+str(re.findall(r'^.*Top LLG.*$',summaryNode.text,re.MULTILINE)[0]).strip())
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
                    
        # Finally a dump of all smartie graphs that were generated
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
        for header, selector in [("Ensemble name","Component/Name"),("Rot Func Z-score","Component/RFZ"),("Trans Func Z-score","Component/TFZ"),("Refined TFZ-equiv","Component/refTFZ"),("Packing clashes","Component/PAK"),("Log likelihood gain","Component/LLG"),("Overall LLG","Component/overallLLG")]:
            if len (solutionNode.findall(selector)) > 0:
                table.addData(title = header, data=[component.text for component in solutionNode.findall(selector)],tip="Stored in xml column {0}".format(selector))
                #print header, [component.text for component in solutionNode.findall(selector)]

    def addResults(self, parent=None):
        if parent is None: parent=self
        resultsNodes = self.xmlnode.findall('PhaserMrSolutions/Solutions')
        if len(resultsNodes) > 0:
            dataLabelTuples = [('SPG','Space group'), ('TFZ','Trans. Z-score'), ('TFZeq','Refined Trans. Z-score'), ('ORIG_LLG','Initial LLG'), ('LLG','Refined LLG'),  ('ORIG_R','Initial Rfactor'), ('R','Refined Rfactor'),('PAK','Clashes')]
            if len(resultsNodes[0].findall('Solution')) > 1:
                parent.addText(text='Multiple solutions found',outputXML=self.outputXml, internalId='MultipleSolText')
                parent.append('<br/>')
                graph = parent.addFlotGraph(title='Solutions',xmlnode=resultsNodes[0], select='Solution',style='height:200px; width:600px;',outputXml=self.outputXml, internalId='PhaserSolutions')
                plotObject = graph.addPlotObject()
                plotObject.append('title','Solutions')
                plotObject.append('plottype','xy')
                plotObject.append('yrange', rightaxis='true')
                data=[i for i in range(len(resultsNodes[0].findall('Solution')))]
                graph.addData(title='Solution number'.replace(' ','_').replace('.','_'), data=data )
                iLine = 0
                for selector,label in dataLabelTuples:
                    if selector != 'SPG':
                        graph.addData(title = label.replace(' ','_').replace('.','_'), select = selector)
                        if selector == 'R' or selector == 'ORIG_R':
                            plotLine = plotObject.append('plotline',xcol=1,ycol=iLine+2,rightaxis='true')
                        else:
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
        



