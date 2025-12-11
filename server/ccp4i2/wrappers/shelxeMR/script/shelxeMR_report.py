#=======================================================================================
#
#    shelxeMR_report.py : shelxeMR_report(Report)
#    
#    Author  : Kyle Stevenson,STFC
#    Created : 14th April 2016, KJS
#
#    Class to create reports for MR solutions using Shelxe
#
#=======================================================================================

from ccp4i2.report.CCP4ReportParser import *
import sys

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
        
        parent.append("<p>Results for Shelxe for MR Run</p>")
        
        if not SHELMR_DYN:
       	    bestCC = float(self.xmlnode.findall('.//RunInfo/BestCycle/BestCC')[0].text)
       	    parent.append("<p>The optimal Correlation Coefficient during the run was found to be :- %f </p>"%(bestCC))

        graph_height = 300
        graph_width = 500
        
        graph = parent.addFlotGraph( title="Results by Shelxe Trace Cycle", select=".//RunInfo/Cycle",style="height:%dpx; width:%dpx; float:left; border:0px;" % (graph_height, graph_width),outputXml=self.outputXml,internalId="SummaryGraph" )
        graph.addData (title="Cycle",  select="NCycle" )
        graph.addData (title="Corr.Coef.", select="CorrelationCoef")
        graph.addData (title="Corr.Coef.", select="AverageChainLen")
        
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
        
"""
     CTaskShelxeMR.py: CCP4 GUI Project
     Copyright (C) 2015 STFC

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
    
