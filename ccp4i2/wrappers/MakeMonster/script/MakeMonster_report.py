from ccp4i2.report.CCP4ReportParser import *
import sys
from lxml import etree

class MakeMonster_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'MakeMonster'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)
