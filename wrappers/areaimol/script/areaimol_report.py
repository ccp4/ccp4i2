from report.CCP4ReportParser import *
import sys
import base64

class areaimol_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'areaimol'
    RUNNING = False
    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")

        if jobStatus not in ["Running", "Running remotely"]:
            try:
                summaryText = ""
                if len(self.xmlnode.findall(".//SummaryText"))>0:
                    xmlPath = './/SummaryText'
                    xmlNodes = self.xmlnode.findall(xmlPath)
                    for node in xmlNodes:
                        summaryText += base64.b64decode(node.text).decode()
                if summaryText:
                    fold = self.addFold(label="Summary", initiallyOpen=True)
                    fold.addPre(text=summaryText)
            except:
                pass

        fold = self.addFold(label="Areaimol log file")
        if len(self.xmlnode.findall(".//LogText"))>0:
            xmlPath = './/LogText'
            xmlNodes = self.xmlnode.findall(xmlPath)
            for node in xmlNodes:
                fold.addPre(text=base64.b64decode(node.text).decode())
