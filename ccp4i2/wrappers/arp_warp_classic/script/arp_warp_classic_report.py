
from ....report.CCP4RvapiParser import RvapiReport

class arp_warp_classic_report(RvapiReport):
  TASKNAME = 'arp_warp_classic'

# WATCHED_FILE  = 'program.xml'
  WATCHED_FILE  = 'report/task.tsk'
  RVAPI_XML     = 'program.xml'
  RVAPI_DIR     = 'report'

  RVAPI_MERGE_TABS = True
  RVAPI_STD_HEADER = False
