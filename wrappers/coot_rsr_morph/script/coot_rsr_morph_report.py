from report.CCP4ReportParser import *
import sys
import base64

class coot_rsr_morph_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'coot_rsr_morph'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        #watersFoundPath = './/coot_find_waters/WatersFound'
        #watersFoundString = xmlnode.findall(watersFoundPath)[0].text
        #self.addText(text='Number of waters found: ' + watersFoundString)
        clearingDiv = self.addDiv(style="clear:both;")
        summaryFold = self.addFold(label="Log File", initiallyOpen=False)
        logNodes = self.xmlnode.findall('Log')
        if len(logNodes)>0:
             logText = base64.b64decode(logNodes[0].text).decode()
             summaryDiv = summaryFold.addDiv()
             summaryDiv.append(logText.replace("\n","<br/>").replace("\033[1m","<b>").replace("\033[0m","</b>"))
