from ....report.CCP4ReportParser import Report


class ZZPluginNameZZ_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ZZPluginNameZZ'
    RUNNING = False

    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        self.addDiv(style="clear:both;")
