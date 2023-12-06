from report.CCP4ReportParser import *
import sys
from lxml import etree
from wrappers.import_mosflm.script import import_mosflm_report

class imosflm_report(import_mosflm_report.import_mosflm_report):
    TASKNAME='imosflm'
