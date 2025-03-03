from report.CCP4ReportParser import *
import sys
from lxml import etree

class prosmart_martin_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'prosmart_martin'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        frameset = self.append('<iframe src="ProSMART_Results.html" style="width:760px;"/>')
