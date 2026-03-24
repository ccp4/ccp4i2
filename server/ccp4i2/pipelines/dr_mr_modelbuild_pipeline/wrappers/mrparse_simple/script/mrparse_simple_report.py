import os

from ccp4i2.report import Report


class mrparse_simple_report(Report):

    TASKNAME = 'mrparse_simple'
    USEPROGRAMXML = False
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, jobStatus=jobStatus, **kw)
        if jobStatus is None or jobStatus.lower() == 'nooutput':
            return
        self.outputXml = self.jobStatus is not None and self.jobStatus.lower().count('running')
        if self.jobStatus is not None and not self.jobStatus.lower().count('running'):
            self.outputXml = False
        self.defaultReport()
        return

    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        parent.append("<p>Finished running MrParse</p>")
        basepath = self.jobInfo['fileroot']
        mrparse_rep = os.path.join(basepath, "mrparse_0", 'mrparse.html')

        ResultsI2Folder = parent.addFold(label='MrParse Reports', initiallyOpen=True)
        if os.path.exists(mrparse_rep):
            ResultsI2Folder.addFileLink(
                label='Open MrParse Results',
                relativePath='mrparse_0/mrparse.html',
                fileType='html',
            )
        else:
            ResultsI2Folder.append('<p>MrParse report not found</p>')
        return
