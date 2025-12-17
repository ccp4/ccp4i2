
from ccp4i2.report import Report


class dials_image_report(Report):

    TASKNAME = 'dials_image'

    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        #Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,jobStatus=jobStatus,**kw)
        #if jobStatus is None or jobStatus.lower() is 'nooutput':
        return
        #self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.append("<p>Running the DIALS Image Viewer</p>")
