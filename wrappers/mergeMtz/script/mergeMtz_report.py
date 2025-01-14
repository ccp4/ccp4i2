from ....report import CCP4ReportParser


class mergeMtz_report(CCP4ReportParser.Report):
  TASKNAME = 'mergeMtz'
  USEPROGRAMXML = False
  def __init__(self,xmlnode=None,jobInfo={},**kw):
    CCP4ReportParser.Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
