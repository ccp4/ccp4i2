<<<<<<<< HEAD:server/ccp4i2/wrappers/coot_find_waters/script/coot_find_waters_report.py
from ccp4i2.report import Report

========
from report.CCP4ReportParser import Report
>>>>>>>> ab6a32677 (Replacing Coot scripts with CHAPI (#106)):wrappers/coot_find_waters/script/coot_find_waters_report.py

class coot_find_waters_report(Report):
    TASKNAME = 'coot_find_waters'
    RUNNING = False

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        super(). __init__(xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        watersFoundPath = './/WatersFound'
        if len(xmlnode.findall(watersFoundPath)) > 0:
            watersFoundString = xmlnode.findall(watersFoundPath)[0].text
            self.addText(text='Number of waters found: ' + watersFoundString)
        else:
            self.addText(text='Find waters complete. Output file contains any added waters')
