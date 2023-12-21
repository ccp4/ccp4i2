"""
    sheetbend_report.py: CCP4 GUI Project
    
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

from __future__ import print_function
from report.CCP4ReportParser import Report
import sys

class sheetbend_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'sheetbend'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.defaultReport()
        
    def defaultReport(self, parent=None):
        if parent is None: parent = self

        clearingDiv = parent.addDiv(style="clear:both;")
        parent.append( "<p>Note: R factors and free R factors are only comparable for cycles where the resolution is the same.</p>" )

        tableDiv = parent.addDiv(style="float:left;border:0px;")
        table = tableDiv.addTable(select=".", transpose=False, id='cycles') 
        try:
          for title,select,expr in [[ "Cycle" , "Cycles/Cycle/Number", "x" ],
                                    [ "Resolution" , "Cycles/Cycle/Resolution", "x" ],
                                    [ "R<sub>Work</sub>"  ,   "Cycles/Cycle/Rwork", "round(x,3)"],
                                    [ "R<sub>Free</sub>" ,    "Cycles/Cycle/Rfree","round(x,3) if float(x)>0.0 else '-' " ] ]:
            table.addData(title=title,select=select,expr=expr)
        except Exception as e:
            print( "ERROR sheetbend_report ", e, file=sys.stderr )

        clearingDiv = parent.addDiv(style="clear:both;")
