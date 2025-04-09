from ....report.CCP4ReportParser import Report
from ...import_mosflm.script import import_mosflm_report


class imosflm_report(import_mosflm_report.import_mosflm_report):
    TASKNAME='imosflm'

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
