from ccp4i2.report import Report


class acorn_report(Report):
    TASKNAME= "acorn"
    RUNNING = False
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        #Report.__init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw) # This caused epic headaches, along with the commandScript
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)  # Note if you want the weird coot error report swap these lines out
        
        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None: parent = self

        results = self.addResults()
        
        parent.append("<p>Results for Acorn Run</p>")
        # The numbers, not only the plot: where the correlation started, its
        # best and where it ended (cycle 0 is the starting point, before any
        # dynamic density modification).
        cycles = [(c.findtext("NCycle"), c.findtext("CorrelationCoef"))
                  for c in self.xmlnode.findall(".//RunInfo/Cycle")] if self.xmlnode is not None else []
        cycles = [(int(n), float(cc)) for n, cc in cycles if n and cc and int(n) > 0]
        if cycles:
            best = max(cycles, key=lambda c: c[1])
            parent.addText(text="Correlation coefficient %.3f after cycle %d, best %.3f "
                           "(cycle %d), final %.3f after %d cycles." % (
                               cycles[0][1], cycles[0][0], best[1], best[0],
                               cycles[-1][1], cycles[-1][0]))
        graph_height = 300
        graph_width = 500
        
        graph = parent.addFlotGraph( title="Results by Cycle", select=".//RunInfo/Cycle",style="height:%dpx; width:%dpx; float:left; border:0px;" % (graph_height, graph_width) )
        graph.addData (title="Cycle",  select="NCycle" )
        graph.addData (title="Cycle",  select="CorrelationCoef" )
        
        p = graph.addPlotObject()
        p.append('title', 'Correlation Coefficient by Cycle')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        p.append('ylabel','Corr.Coef.')
        
        l = p.append('plotline',xcol=1,ycol=2)
        l.append('label','Correlation Coefficient')
        l.append('colour','teal')
