import sys

from ccp4i2.report.CCP4ReportParser import Report


class sheetbend_report(Report):
    TASKNAME = 'sheetbend'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)

        if jobStatus is None or jobStatus.lower() == 'nooutput': return
        self.defaultReport()
        
    def defaultReport(self, parent=None, select="."):
        if parent is None: parent = self

        clearingDiv = parent.addDiv(style="clear:both;")
        parent.append( "<p>Note: R factors and free R factors are only comparable for cycles where the resolution is the same.</p>" )

        tableDiv = parent.addDiv(style="float:left;border:0px;")
        print("##################################################")
        print("##################################################")
        print(select)
        print("##################################################")
        print("##################################################")
        table = tableDiv.addTable(select=select, transpose=False, id='cycles') 
        try:
          for title,select,expr in [[ "Cycle" , "Cycles/Cycle/Number", "x" ],
                                    [ "Resolution" , "Cycles/Cycle/Resolution", "x" ],
                                    [ "R<sub>Work</sub>"  ,   "Cycles/Cycle/Rwork", "round(x,3)"],
                                    [ "R<sub>Free</sub>" ,    "Cycles/Cycle/Rfree","round(x,3) if float(x)>0.0 else '-' " ] ]:
            table.addData(title=title,select=select,expr=expr)
        except Exception as e:
            print( "ERROR sheetbend_report ", e, file=sys.stderr )

        clearingDiv = parent.addDiv(style="clear:both;")
