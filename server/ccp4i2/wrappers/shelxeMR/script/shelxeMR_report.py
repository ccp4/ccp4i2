from ccp4i2.report import Report

SHELMR_DYN = True

class shelxeMR_report(Report):

    TASKNAME = 'shelxeMR'
    RUNNING = SHELMR_DYN
    SEPARATEDATA=True
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return

        self.outputXml = self.jobStatus is not None and self.jobStatus.lower().count('running')
        if self.jobStatus is not None and not self.jobStatus.lower().count('running'): self.outputXml = False
        
        self.defaultReport()
    
    def defaultReport(self, parent=None):
        if parent is None: parent = self

        results = self.addResults()
        
        # The best trace, in numbers: the plot alone gave none.
        best = self.xmlnode.find('.//RunInfo/BestCycle') if self.xmlnode is not None else None
        if best is not None and best.findtext('BestCC'):
            try:
                residues = int(float(best.findtext('ChainLen')) * int(best.findtext('NumChains')))
                parent.append(
                    "<p>Best trace: cycle %s, CC %s%% for the traced structure against the "
                    "data, %d residues in %s chain(s).</p>"
                    % (best.findtext('BCycle'), best.findtext('BestCC'), residues,
                       best.findtext('NumChains')))
            except (TypeError, ValueError):
                pass

        graph_height = 300
        graph_width = 500
        
        graph = parent.addFlotGraph( title="Results by Shelxe Trace Cycle", select=".//RunInfo/Cycle",style="height:%dpx; width:%dpx; float:left; border:0px;" % (graph_height, graph_width),outputXml=self.outputXml,internalId="SummaryGraph" )
        graph.addData (title="Cycle",  select="NCycle" )
        graph.addData (title="Corr.Coef.", select="CorrelationCoef")
        graph.addData (title="Average chain length", select="AverageChainLen")
        
        p = graph.addPlotObject()
        p.append('title', 'Correlation Coefficient by Trace Cycle')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        p.append('ylabel','Corr.Coef.')
        
        l = p.append('plotline',xcol=1,ycol=2)
        l.append('label','By residue')
        l.append('colour','teal')
        
        p = graph.addPlotObject()
        p.append('title', 'Average Chain Length by Trace Cycle')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        p.append('ylabel','<Chain length>')
        
        l = p.append('plotline',xcol=1,ycol=3)
        l.append('label','By residue')
        l.append('colour','red')
