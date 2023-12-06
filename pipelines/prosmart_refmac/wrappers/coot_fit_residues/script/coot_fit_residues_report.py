import sys
from report.CCP4ReportParser import *

class coot_fit_residues_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'coot_fit_residues'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        if jobStatus is not None and jobStatus.lower() == 'nooutput':
            return
        
        progressGraph = self.addGraph(title="Running fit residues",select=".//Coot_fit_residues/Table/row",style="height:250px; width:400px;float:left;")
        progressGraph.addData(title="Residue_number",         select="Col_0")
        progressGraph.addData(title="Initial_bond_deviation", select="Col_1")
        progressGraph.addData(title="Final_bond_deviation",   select="Col_2")
        plot = progressGraph.addPlotObject()
        plot.append('title','Running fit residues')
        plot.append('plottype','xy')
        plot.append('xintegral','true')
        for coordinate, colour in [(2,'blue'),(3,'green')]:
            plotLine = plot.append('plotline',xcol=1,ycol=coordinate,colour=colour)

