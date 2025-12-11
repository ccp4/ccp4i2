#Just use the xia2_run_report code - this is to get round any
#places (such as CProjectViewer) that does not know that the taskname
#has changed
from ccp4i2.pipelines.import_xia2.wrappers.xia2_run.script import xia2_run_report

class import_xia2_report(xia2_run_report.xia2_run_report):
    TASKNAME = 'import_xia2'
    pass
