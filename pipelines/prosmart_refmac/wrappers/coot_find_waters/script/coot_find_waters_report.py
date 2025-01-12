from ......report.CCP4ReportParser import Report


class coot_find_waters_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'coot_find_waters'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        watersFoundPath = './/WatersFound'
        if len(xmlnode.findall(watersFoundPath)) > 0:
            watersFoundString = xmlnode.findall(watersFoundPath)[0].text
            self.addText(text='Number of waters found: ' + watersFoundString)
        else:
#We can get here on Windows because the Coot log is not written to log.txt
            self.addText(text='Find waters complete. Output file contains any added waters')
