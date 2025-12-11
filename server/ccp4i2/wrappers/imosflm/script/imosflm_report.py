from ccp4i2.report.CCP4ReportParser import *
import sys
from lxml import etree
from ccp4i2.wrappers.import_mosflm.script import import_mosflm_report

class imosflm_report(import_mosflm_report.import_mosflm_report):
    TASKNAME='imosflm'

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__( self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
