"""
    i2Dimple_report.py: CCP4 GUI Project
    
    This library is free software: you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    version 3, modified in accordance with the provisions of the
    license to address the requirements of UK law.
    
    You should have received a copy of the modified GNU Lesser General
    Public License along with this library.  If not, copies may be
    downloaded from http://www.ccp4.ac.uk/ccp4license.php
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.
    """

from ccp4i2.report.CCP4ReportParser import Report
import sys

class i2Dimple_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'i2Dimple'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        clearingDiv = self.addDiv(style="clear:both;")
        self.addDefaultReport(self)
        
    def addDefaultReport(self, parent=None):
        if parent is None: parent=self
        
        # Report here if Dimple's pointless run identified need for a reindexing
        reindexNodes = self.xmlnode.findall(".//REINDEX")
        if len(reindexNodes) > 0:
            newFold = parent.addFold(label="POINTLESS result", initiallyOpen=True)
            newFold.addPre(style="font-size:125%; font-color:red;", text="DIMPLE identified a need to reindex.")
            reindexText = "New reflection and FreeR (if given) have been output with operator {}".format(reindexNodes[0].text)
            newFold.addPre(style="font-size:125%; font-color:red;", text=reindexText)
        
        #Add a summary to flag where PHASER may have been used
        if len(self.xmlnode.findall("PHASER")) > 0:
            newFold = parent.addFold(label="PHASER used to solve structure", initiallyOpen=True)
            newFold.addPre(text=self.xmlnode.findall("PHASER")[0].text)
        
        #Add a tabular and graphical analysis of the refmac progress
        if len(self.xmlnode.findall("REFMAC")) > 0:
            newFold = parent.addFold(label="Refmac progress", initiallyOpen=True)
            
            #Add a summary table
            tableDiv = newFold.addDiv(style="float:left;width:300px;");
            tableDiv.addText(text="Initial and final values",style="font-size:125%");
            table = tableDiv.addTable(transpose=True,style="float:right;")
            table.addData(title="Parameter",data=["Initial","Final"])
            table.addData(title="R-value",data=[self.xmlnode.findall(".//new_cycle/r_factor")[0].text,
                                                self.xmlnode.findall(".//new_cycle/r_factor")[-1].text])
            table.addData(title="R-free",data=[self.xmlnode.findall(".//new_cycle/r_free")[0].text,
                                                self.xmlnode.findall(".//new_cycle/r_free")[-1].text])
            table.addData(title="rmsBOND",data=[self.xmlnode.findall(".//new_cycle/rmsBOND")[0].text,
                                                self.xmlnode.findall(".//new_cycle/rmsBOND")[-1].text])
            table.addData(title="rmsANGLE",data=[self.xmlnode.findall(".//new_cycle/rmsANGLE")[0].text,
                                                self.xmlnode.findall(".//new_cycle/rmsANGLE")[-1].text])
            table.addData(title="rmsCHIRAL",data=[self.xmlnode.findall(".//new_cycle/rmsCHIRAL")[0].text,
                                                self.xmlnode.findall(".//new_cycle/rmsCHIRAL")[-1].text])
                        
            #add a summary graph
            progressGraph = newFold.addFlotGraph(title="Refmac progress", select=".//stats_vs_cycle", style="width:400px;border:0px;float:left;")
            newFold.addDiv(style="clear:both;")
            
            #Create a data column for the cycle number
            nCycles = len(self.xmlnode.findall(".//REFMAC/Overall_stats/stats_vs_cycle/new_cycle"))
            #progressGraph.addData(title="Cycle", data=[iCycle for iCycle in range(nCycles)])
            progressGraph.addData(title="Cyle", select="new_cycle/cycle")
            progressGraph.addData(title="R_Factor", select="new_cycle/r_factor")
            progressGraph.addData(title="R_Free",  select="new_cycle/r_free")
            progressGraph.addData(title="rmsBOND",  select="new_cycle/rmsBOND")
            progressGraph.addData(title="rmsANGLE",  select="new_cycle/rmsANGLE")
            progressGraph.addData(title="rmsCHIRAL",  select="new_cycle/rmsCHIRAL")
            
            p = progressGraph.addPlotObject()

            p.append('title', 'R-factors and geometry VS cycle')
            p.append('plottype','xy')
            p.append('xintegral','true')
            p.append('xlabel','Cycle')
            p.append('ylabel','R-factors')
            
            l = p.append('plotline',xcol=1,ycol=2)
            l.append('label','R-factor')
            l.append('colour','blue')
            
            l = p.append('plotline',xcol=1,ycol=3)
            l.append('label','Free-R')
            l.append('colour','green')

            l = p.append('plotline',xcol=1,ycol=4,rightaxis='true')
            l.append('label','rmsBOND')
            l.append('colour','red')
            
            p.append('yrange', rightaxis='false', min='0.1',max='0.6')
            p.append('yrange', rightaxis='true', min='0.0',max='0.05')
        
        #Add a tabular and graphical analysis of any blobs
        if len(self.xmlnode.findall("find-blobs")) > 0:
            newFold = parent.addFold(label="Blobs found", initiallyOpen=True)
            
            #Add a summary table
            tableDiv = newFold.addDiv(style="float:left;width:300px;");
            tableDiv.addText(text="Blobs found:",style="font-size:125%");
            table = tableDiv.addTable(style="float:right;")
            score = [e.text for e in self.xmlnode.findall(".//find-blobs/Blob/score")]
            x = [e.text for e in self.xmlnode.findall(".//find-blobs/Blob/x")]
            y = [e.text for e in self.xmlnode.findall(".//find-blobs/Blob/y")]
            z = [e.text for e in self.xmlnode.findall(".//find-blobs/Blob/z")]
            table.addData(title="Score",data=score)
            table.addData(title="X",data=x)
            table.addData(title="Y",data=y)
            table.addData(title="Z",data=z)
            
            #add a summary graph
            blobGraph = newFold.addFlotGraph(title="Blob score distribution", select=".//find-blobs", style="width:400px;border:0px;float:left;")
            
            #Create a data column for the cycle number
            nBlobs = len(self.xmlnode.findall(".//Blob"))
            
            blobGraph.addData(title="Rank", data=[iBlob for iBlob in range(nBlobs)])
            blobGraph.addData(title="Score", select="Blob/score")
            
            p = blobGraph.addPlotObject()

            p.append('title', 'Score versus rank')
            p.append('plottype','xy')
            p.append('xintegral','true')
            p.append('xlabel','Rank')
            p.append('ylabel','score')
            
            l = p.append('plotline',xcol=1,ycol=2)
            l.append('label','Score')
            l.append('colour','blue')
            
            newFold.addDiv(style="clear:both;")

        if len(self.xmlnode.findall("LogText")) > 0:
            newFold = parent.addFold(label="Summary", initiallyOpen=True)
            newFold.addPre(text = self.xmlnode.findall("LogText")[0].text)

