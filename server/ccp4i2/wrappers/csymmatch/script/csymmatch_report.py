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

        # The program's XML is <Csymmatch> itself, and phaser_pipeline hands
        # over that element too: ".//Csymmatch/..." looks only below it, so
        # every one of these found nothing and the report said nothing.
        root = self.xmlnode.getroot() if hasattr(self.xmlnode, 'getroot') else self.xmlnode
        node = root if root.tag == 'Csymmatch' else root.find('.//Csymmatch')
        if node is None:
            return

        hand = node.find('ChangeOfHand')
        if hand is not None and "Y" in (hand.text or ''):
            results.append('A change of hand was applied. If the model is more than a substructure, this probably means that no suitable match was found')

        origin = node.find('ChangeOfOrigin')
        if origin is not None:
            results.append('A change of origin was applied, with fractional coordinates '+(origin.text or '').strip())

        segmentNodes = node.findall('Segment')
        if len(segmentNodes) > 0:
            results.append('The structure was grouped into '+str(len(segmentNodes))+' segments for symmetry matching.')

            detailFold = parent.addFold(label='Transformations and scores')
            detailTable = detailFold.addTable(title='transformations and scores',selectNodes=segmentNodes)
            detailTable.addData(title='Range', select='Range')
            detailTable.addData(title='Operator', select='Operator')
            detailTable.addData(title='Shift', select='Shift')
            detailTable.addData(title='Score', select='Score')
            
            graphFold = parent.addFold(label='Normalized scores plot')
            progressGraph = graphFold.addFlotGraph(title="Per segment normalized score",selectNodes=segmentNodes,style="height:250px; width:600px;float:left;border:0px solid white;")
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
