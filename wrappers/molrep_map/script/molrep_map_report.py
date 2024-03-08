from report.CCP4ReportParser import *
import sys

class molrep_map_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'molrep_map'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        self.addText(text='Happily finished')
       


