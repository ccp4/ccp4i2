import sys
from ccp4i2.report.CCP4ReportParser import *

class coot_stepped_refine_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'coot_stepped_refine'
    RUNNING = True
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return

        progressGraph = self.addGraph(title="Running stepped refinment",select=".//Coot_fit_residues/Table/row",style="height:250px; width:400px;float:left;")
        progressGraph.addData(title="Residue_number",         select="Col_0")
        progressGraph.addData(title="Initial_bond_deviation", select="Col_1")
        progressGraph.addData(title="Final_bond_deviation",   select="Col_2")
        plot = progressGraph.addPlotObject()
        plot.append('title','Running stepped refienment')
        plot.append('plottype','xy')
        plot.append('xintegral','true')
        for coordinate, colour in [(2,'blue'),(3,'green')]:
            plotLine = plot.append('plotline',xcol=1,ycol=coordinate,colour=colour)

    def addSmartieGraph(self, select, parent=None, title=None, style=None, plotTuples=None):
        if parent is None:
            parent = self
        cctableNode = self.xmlnode.findall(select)[0]
        if title is None:
            title = cctableNode.get('title')
        graph = parent.addGraph(title=title, select=select, style=style)
        
        headerElements = cctableNode.findall('Header')
        
        labelsToColumnsDict = {}
        for headerElement in headerElements:
            label = headerElement.attrib['label']
            identifier = headerElement.attrib['identifier']
            labelsToColumnsDict[label] = identifier
        iCol = 1
        labelsToGraphIcolDict = {}
        for label, identifier in list(labelsToColumnsDict.items()):
            graph.addData (title=label,  select="row/"+identifier )
            labelsToGraphIcolDict[label] = iCol
            iCol += 1
        
        if plotTuples is None:
            plotTuples = []
            graphs = cctableNode.findall('Graph')
            for smartieGraph in graphs:
                plotTitle = smartieGraph.get('title')
                columns = smartieGraph.findall('Column')
                columnsDict = {}
                for column in columns:
                    columnsDict[column.get('positionInList')] = column.get('label')
                plotLineTuples = []
                colorArray = ['red','green','blue','cyan','magenta','yellow','black','grey']
                for i in range (1, len(columns)):
                    plotLineTuples.append((columnsDict[str(i)], columnsDict[str(0)], colorArray[i-1]))
                plotTuples.append((plotTitle, plotLineTuples))
        
        for plotTitle, plotLineTuples in plotTuples:
            plotObject = graph.addPlotObject()
            plotObject.append('title',plotTitle)
            plotObject.append('plottype','xy')
            for plotLineTuple in plotLineTuples:
                yLabel, xLabel, colour = plotLineTuple
                xcol = labelsToGraphIcolDict[xLabel]
                ycol = labelsToGraphIcolDict[yLabel]
                plotLine = plotObject.append('plotline',xcol=xcol,ycol=ycol)
                plotLine.append('colour',colour)

