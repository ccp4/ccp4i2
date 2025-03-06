
from report.CCP4ReportParser import *

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
        if os.path.exists(mrparse_rep):
            projectid = self.jobInfo.get("projectid", None)
            jobNumber = self.jobInfo.get("jobnumber", None)

            mrparseurl = (
                "/database/?getProjectJobFile?projectId="
                + projectid
                + "?fileName=mrparse_0/mrparse.html?jobNumber="
                + jobNumber
            )

        ResultsI2Folder = parent.addFold(label='MrParse Reports', initiallyOpen=True)
        ResultsI2Folder.append('<br></br>')
        ResultsI2Folder.append('<span style="font-size:110%">Click on the '
                                    'following link to display the'
                                    'browser report for the MrParse job '
                                    '</span>')
        ResultsI2Folder.append('<br></br>')
        #ResultsI2Folder.append('<a href="{0}">Open Results</a>'.format(mrparseurl))
        ResultsI2Folder.append('<a href="{0}">Open Results</a>'.format(mrparse_rep))
        return
