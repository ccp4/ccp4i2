import base64

from ccp4i2.report import Report


class pdbset_ui_report(Report):
    TASKNAME = 'pdbset_ui'
    RUNNING = False
    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(
            self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw
        )
        self.addDiv(style="clear:both;")
        if jobStatus in ["Running", "Running remotely"]:
            self.append("<p><b>The job is currently running.</b></p>")

        fold = self.addFold(label="Log file")
        if len(self.xmlnode.findall(".//LogText"))>0:
            xmlPath = './/LogText'
            xmlNodes = self.xmlnode.findall(xmlPath)
            for node in xmlNodes:
                fold.addPre(text=base64.b64decode(node.text).decode())
