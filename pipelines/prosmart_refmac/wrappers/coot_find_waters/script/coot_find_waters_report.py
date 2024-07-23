from report.CCP4ReportParser import *

class coot_find_waters_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'coot_find_waters'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        watersFoundPath = './/WatersFound'
        watersFoundString = xmlnode.findall(watersFoundPath)[0].text
        self.addText(text='Number of waters found: ' + watersFoundString)
