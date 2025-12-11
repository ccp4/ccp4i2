import sys
from lxml import etree
from ccp4i2.report.CCP4ReportParser import *

class reindex_processed_data_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'reindex_processed_data'
    RUNNING = False
    USEPROGRAMXML = False
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

