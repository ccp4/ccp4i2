from report.CCP4ReportParser import *
import sys

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
