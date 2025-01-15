import os

from ....report.CCP4ReportParser import Report


class buster_report(Report):

    TASKNAME = 'buster'
    USEPROGRAMXML = True
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.defaultReport()
        
    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.addResults()
        parent.append("Buster Refinement complete.")
        bestR = float(self.xmlnode.findall('.//RunInfo/Best/R')[0].text)
        bestRFree = float(self.xmlnode.findall('.//RunInfo/Best/RFree')[0].text)
        sbestR = f'{round(bestR, 3):.3g}'
        sbestRFree = f'{round(bestRFree, 3):.3g}'
        parent.append("<p>Best refinement in BUSTER reached with R/Rfree: %s / %s </p>"%(sbestR,sbestRFree))
        inst = "<h3>BUSTER Plots</h3>"
        inst += "<h3>R/RFree, Log-Likelihood, and RMS Values vs Cycle:</h3>"
        parent.append(inst)
        # Add graphs to report
        graph_height = 450
        graph_width = 600
        graph = parent.addFlotGraph( title="Results by Cycle", select=".//RunInfo/Cycle",style="height:%dpx; width:%dpx; float:left; border:0px;" % (graph_height, graph_width) )
        graph.addData (title="Total Cycles (summed)", select="NCycle" )
        graph.addData (title="R-Factor", select="RFact" )
        graph.addData (title="R-Free",  select="RFree" )
        graph.addData (title="Log-Li. Gain",  select="LLG" )
        graph.addData (title="Log-Li. Gain Free",  select="LLGF" )
        graph.addData (title="RMS Bonds",  select="RMSB" )
        graph.addData (title="RMS Angle",  select="RMSA" )
        # Add R-Factor plot
        p = graph.addPlotObject()
        p.append('title', 'R-Factors by Buster Cycle (inc. all large cycles)')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        p.append('ylabel','R-Factors')
        l = p.append('plotline',xcol=1,ycol=2)
        l.append('label','R-Factor')
        l.append('colour','gold')
        l = p.append('plotline',xcol=1,ycol=3)
        l.append('label','R-Free')
        l.append('colour','lightblue')
        # Add Log-likelihood plot
        p = graph.addPlotObject()
        p.append('title', 'Log-Likelihood Gain by Cycle')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        p.append('ylabel','LLG')
        l = p.append('plotline',xcol=1,ycol=4)
        l.append('label','LLG')
        l.append('colour','green')
        l = p.append('plotline',xcol=1,ycol=5)
        l.append('label','LLGF')
        l.append('colour','lightgreen')
        # Add RMS Values plot
        p = graph.addPlotObject()
        p.append('title', 'RMS Bonds & Angles by Cycle')
        p.append('plottype','xy')
        p.append('xintegral','true')
        p.append('xlabel','Cycle')
        p.append('ylabel','RMS')
        l = p.append('plotline',xcol=1,ycol=6)
        l.append('label','RMS Bonds (x100)')
        l.append('colour','red')
        l = p.append('plotline',xcol=1,ycol=7)
        l.append('label','RMS Angles')
        l.append('colour','lightgreen')

        self.addDiv(style='clear:both;')

        fileRoot = self.jobInfo['fileroot']
        if os.path.exists(os.path.join(fileRoot,"Plots","summary.png")):
            summary_Url = "Plots/summary.png"
            summaryFold = self.addFold(label="summary.png", initiallyOpen=False)
            summaryFold.append("<img src=\"%s\"></img>"%(summary_Url))
            
        if os.path.exists(os.path.join(fileRoot,"Plots","summary_LL.png")):
            summary_LL_Url = "Plots/summary_LL.png"
            summary_LL_Fold = self.addFold(label="summary_LL.png", initiallyOpen=False)
            summary_LL_Fold.append("<img src=\"%s\"></img>"%(summary_LL_Url))
