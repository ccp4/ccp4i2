from ccp4i2.report.CCP4ReportParser import *
import sys

class coot_script_lines_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'coot_script_lines'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        self.addText(text='Happily finished')
        try:
            pdbsWrittenPath = './/number_output_pdbs'
            if len(xmlnode.findall(pdbsWrittenPath)) > 0:
                pdbsWrittenString = xmlnode.findall(pdbsWrittenPath)[0].text
                self.addText(text='Number of PDBs written: ' + pdbsWrittenString)
            cifsWrittenPath = './/number_output_cifs'
            if len(xmlnode.findall(cifsWrittenPath)) > 0:
                cifsWrittenString = xmlnode.findall(cifsWrittenPath)[0].text
                self.addText(text='Number of CIFs written: ' + cifsWrittenString)
        except:
            pass

        tableRows = self.xmlnode.findall('.//Table/row')
        if len(tableRows) > 0:
            self.append('<br/>')
            progressGraph = self.addGraph(title="Running fit residues",select=".//Table/row",style="height:250px; width:400px;")
            progressGraph.addData(title="Residue_number",         select="Col_0")
            progressGraph.addData(title="Initial_bond_deviation", select="Col_1")
            progressGraph.addData(title="Final_bond_deviation",   select="Col_2")
            plot = progressGraph.addPlotObject()
            plot.append('title','Running fit residues')
            plot.append('plottype','xy')
            plot.append('xintegral','true')
            for coordinate, colour in [(2,'blue'),(3,'green')]:
                plotLine = plot.append('plotline',xcol=1,ycol=coordinate,colour=colour)
        newDiv = self.addDiv(style='clear:both;')
