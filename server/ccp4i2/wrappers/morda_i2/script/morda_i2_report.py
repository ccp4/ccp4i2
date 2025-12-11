
from ccp4i2.report.CCP4RvapiParser import RvapiReport

class morda_i2_report(RvapiReport):
  TASKNAME = 'morda_i2'
  WATCHED_FILE = 'program.xml'
  RVAPI_XML = 'program.xml'
  RVAPI_DIR = 'report'

  RVAPI_MERGE_TABS = False
  RVAPI_STD_HEADER = False

