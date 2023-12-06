from report.CCP4ReportParser import *

class chltofom_report(Report):
    TASKNAME = 'chltofom'
    USEPROGRAMXML = False
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
      Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
