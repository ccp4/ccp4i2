from report.CCP4ReportParser import *
import sys
#from lxml import etree

class splitMtz_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'splitMtz'
    RUNNING = False
    USEPROGRAMXML=False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        
        
