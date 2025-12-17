from ccp4i2.report.CCP4ReportParser import *


class phaser_MR_FRF_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'phaser_MR_FRF'
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
        
        if len(self.xmlnode.findall('PhaserMrSolutions/Solutions')) > 0:
            compareSolutionsFold = parent.addFold(label='Comparison of solutions',initiallyOpen=True)
            self.addResults(parent=compareSolutionsFold)
       
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
                if summaryBlockTitle == "MOLECULAR REPLACEMENT ROTATION FUNCTION":
                    try:
                        summaryBlockTitle += (': '+str(re.findall(r'^.*Top RF.*$',summaryNode.text,re.MULTILINE)[0]).strip())
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

    def addResults(self, parent=None):
        if parent is None: parent=self
        resultsNodes = self.xmlnode.findall('PhaserMrSolutions/Solutions')
        if len(resultsNodes) > 0:
            dataLabelTuples = [('EULER','Euler Angles'), ('RFZ','Rotation Function Z-score'), ('DEEP','Deep Search')]
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
                    if selector == 'RFZ':
                        graph.addData(title = label.replace(' ','_').replace('.','_'), select = selector)
                        plotLine = plotObject.append('plotline',xcol=1,ycol=iLine+2)
                        iLine += 1
            elif len(resultsNodes[0].findall('Solution')) == 1:
                parent.addText(text='Unique solution found :-)')
            elif len(resultsNodes[0].findall('Solution')) == 0:
                parent.addText(text='No solutions found',style='color:red;')
                return
            table = parent.addTable(xmlnode=resultsNodes[0], select='Solution',style='width:600px;')
            for selector,label in dataLabelTuples:
                if selector == 'DEEP': table.addData(title = label, select = selector)
                elif selector == 'EULER': table.addData(title = label, select = selector, function=self.format_euler)
                else: table.addData(title = label, select = selector, expr='"{0:.2f}".format(float(x))')
    
    def format_euler(self, euler):
        return ', '.join(["{0:.2f}".format(float(x)) for x in euler.replace('(','').replace(')','').split(',')])


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



