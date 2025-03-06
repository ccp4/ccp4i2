from report.CCP4ReportParser import *

class cif2mtz_report(Report):
    TASKNAME = 'cif2mtz'
    USEPROGRAMXML = False
    
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
      Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
