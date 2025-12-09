from ccp4i2.report.CCP4ReportParser import *
import sys

class coot_rsr_morph_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'coot_rsr_morph'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        #watersFoundPath = './/coot_find_waters/WatersFound'
        #watersFoundString = xmlnode.findall(watersFoundPath)[0].text
        #self.addText(text='Number of waters found: ' + watersFoundString)
        self.addText(text='Coot real space morphing finished. Full reporting is not yet available in this task.')
