from ccp4i2.report import Report


class csymmatch_report(Report):
    TASKNAME = 'csymmatch'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        if jobStatus is not None and jobStatus.lower() == 'nooutput': return
        self.drawContent(jobStatus, self)

    def drawContent(self, jobStatus, parent=None):
        if parent is None: parent = self

        results = parent.addResults()

        if len(self.xmlnode.findall(".//Csymmatch/ChangeOfHand")) > 0:
          if "Y" in self.xmlnode.findall(".//Csymmatch/ChangeOfHand")[0].text:
            results.append('A change of hand was applied. If the model is more than a substructure, this probably means that no suitable match was found')

        if len(self.xmlnode.findall(".//Csymmatch/ChangeOfOrigin")) > 0:
          results.append('A change of origin was applied, with fractional coordinates '+self.xmlnode.findall(".//Csymmatch/ChangeOfOrigin")[0].text)

        segmentNodes = self.xmlnode.findall(".//Csymmatch/Segment")
        if len(segmentNodes) > 0:
            results.append('The structure was grouped into '+str(len(segmentNodes))+' segments for symmetry matching.')

            detailFold = parent.addFold(label='Transformations and scores')
            detailTable = detailFold.addTable(title='transformations and scores',select=".//Csymmatch/Segment")
            detailTable.addData(title='Range', select='Range')
            detailTable.addData(title='Operator', select='Operator')
            detailTable.addData(title='Shift', select='Shift')
            detailTable.addData(title='Score', select='Score')
            
            graphFold = parent.addFold(label='Normalized scores plot')
            progressGraph = graphFold.addFlotGraph(title="Per segment normalized score",select=".//Csymmatch/Segment",style="height:250px; width:600px;float:left;border:0px solid white;")
            segmentNumbers = [i for i in range(len(segmentNodes))]
            progressGraph.addData(title="Segment number", data=segmentNumbers)
            progressGraph.addData(title="Normalized_score", select="Score")
            plot = progressGraph.addPlotObject()
            plot.append('title','Normalized scores of segments')
            plot.append('plottype','xy')
            plot.append('xintegral','true')
            for coordinate, colour in [(2,'blue')]:
                plotLine = plot.append('plotline',xcol=1,ycol=coordinate,colour=colour)
            parent.addDiv(style='clear:both')
