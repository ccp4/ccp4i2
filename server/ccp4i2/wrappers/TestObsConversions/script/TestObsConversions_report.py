from ccp4i2.report.CCP4ReportParser import *
import sys

class TestObsConversions_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'TestObsConversions'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style="clear:both;")
        self.defaultReport()
        
    def defaultReport(self, parent=None):
        if parent is None: parent = self

        table = parent.addTable(title='Analysis of files', select ='.//TestObsConversions/File')
        for property in ('Role','ContentFlag','Columns','Path'):
            table.addData(title=property,select=property)
        
