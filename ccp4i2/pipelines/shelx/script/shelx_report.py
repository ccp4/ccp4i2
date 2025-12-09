from ccp4i2.report.CCP4ReportParser import *
from ccp4i2.pipelines.crank2.script import crank2_report

class shelx_report(crank2_report.crank2_report):
  TASKNAME="shelx"
  RUNNING = True 
