from ccp4i2.report import Report


class phaser_MR_PAK_report(Report):
    TASKNAME = 'phaser_MR_PAK'
    RUNNING = False
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
        
        self.drawAdvisories(parent=parent)
        self.drawWarnings(parent=parent)
        
        bestSolFold = parent.addFold(label='Elements and scores of current solution', brief='Current soln.', initiallyOpen=True)
        if len(self.xmlnode.xpath('PhaserCurrentBestSolution')) > 0:
            self.addBestSolution(parent=bestSolFold)
        else:
            bestSolFold.addText(text='No solution found yet', style='color:red;',outputXML=self.outputXml, internalId='BestSolText')
        if len(self.xmlnode.xpath('PhaserMrSolutions/Solutions')) > 0:
            compareSolutionsFold = parent.addFold(label='Comparison of solutions', brief='All solns', initiallyOpen=True)
            self.addResults(parent=compareSolutionsFold)
        analysisFold = parent.addFold(label='Analysis of composition and data', brief='Comp/data', initiallyOpen=False)
        searchPathFold = parent.addFold(label='Search strategy employed by PHASER', brief='Search tree', initiallyOpen=False)
        indentText = "...."
        addIndent = False
        removeIndent = False
        import re
        for iNode, summaryNode in enumerate(self.xmlnode.xpath('Summary')):
            destination = searchPathFold
             #Make a result folder for any summmary blocks in the XML
            try:
                moduleLines = [line for line in summaryNode.text.split('\n') if "Phaser Module:" in line]
                summaryBlockTitle = moduleLines[0][19:-9].strip()
                
                if summaryBlockTitle == "MOLECULAR REPLACEMENT PACKING ANALYSIS":
                    try:
                        summaryBlockTitle += (': '+str(re.findall(r'^.*accepted of.*$',summaryNode.text,re.MULTILINE)[0]).strip())
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


    def drawAdvisories(self, parent=None):
        if parent is None: parent = self
        advisories = self.xmlnode.xpath('//PhaserAdvisories/Advisory')
        if len(advisories)>0:
            parent.append('<br/>')
            for advisory in advisories:
                parent.addText(text=advisory.text,style='color:orange;',outputXML=self.outputXml, internalId='AdvisoriesText')
                parent.append('<br/>')
        return

    def drawWarnings(self, parent=None):
        if parent is None: parent = self
        warnings = self.xmlnode.xpath('//PhaserWarnings/Warning')
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
        solutionNode = self.xmlnode.xpath0('PhaserCurrentBestSolution/Solution')
        
        summaryText = ""
        if len(solutionNode.xpath('./spaceGroup')) > 0:
            spaceGroup = solutionNode.xpath('./spaceGroup')[0].text
            summaryText += 'Current best solution has spacegroup %s'% spaceGroup
        if len(solutionNode.xpath('./overallLLG')) > 0:
            summaryText +=' Refined overall LLG: %s'%(solutionNode.xpath('./overallLLG')[0].text)
        if len(summaryText)>0: parent.addText(text=summaryText,outputXml=self.outputXml,internalId='BestSolutionSummary')
        
        table = parent.addTable(xmlnode=solutionNode, outputXml=self.outputXml, internalId='BestSolutionBreakdown')
        for header, selector in [("Ensemble name","Component/Name"),("Rot Func Z-score","Component/RFZ"),("Trans Func Z-score","Component/TFZ")]:
            if len (solutionNode.xpath(selector)) > 0:
                table.addData(title = header, data=[component.text for component in solutionNode.xpath(selector)],tip="Stored in xml column {0}".format(selector))
                #print header, [component.text for component in solutionNode.xpath(selector)]

    def addResults(self, parent=None):
        if parent is None: parent=self
        resultsNodes = self.xmlnode.xpath('PhaserMrSolutions/Solutions')
        if len(resultsNodes) > 0:
            dataLabelTuples = [('SPG','Space group'), ('TFZ','Trans. Z-score')]
            if len(resultsNodes[0].xpath('Solution')) > 1:
                parent.addText(text='Multiple solutions found',outputXML=self.outputXml, internalId='MultipleSolText')
                parent.append('<br/>')
            elif len(resultsNodes[0].xpath('Solution')) == 1:
                parent.addText(text='Unique solution found :-)',outputXML=self.outputXml, internalId='UniqueSolText')
            elif len(resultsNodes[0].xpath('Solution')) == 0:
                parent.addText(text='No solutions found',style='color:red;',outputXML=self.outputXml, internalId='NoSolText')
                return
            table = parent.addTable(xmlnode=resultsNodes[0], select='Solution',style='width:600px;', outputXml=self.outputXml, internalId='ResultsTable')
            for selector,label in dataLabelTuples:
                if selector == 'SPG': table.addData(title = 'Space group', select = selector)
                else: table.addData(title = label, select = selector, expr='"{0:.2f}".format(float(x))')



