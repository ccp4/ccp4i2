
from ccp4i2.report.CCP4ReportParser import *

class findmyseq_report(Report):

    TASKNAME = 'findmyseq'
    USEPROGRAMXML = True
    SEPARATEDATA = True

    def __init__(self, xmlnode=None, jobInfo={}, jobStatus=None, **kw):
        Report.__init__(self, xmlnode=xmlnode, jobInfo=jobInfo, **kw)
        self.defaultReport()
        
    def defaultReport(self, parent=None):
        if parent is None:
            parent = self
        jobDirectory = self.jobInfo['fileroot']
        parent.addResults()
        parent.append("FindMySequence run completed.")
    
        tableFoldmr = parent.addFold(label='FindMySequence Log', initiallyOpen=False)
        tableFoldmr.append('Results from FindMySequence<br/>')


        # This is the log file
        if not os.path.isfile(os.path.join(jobDirectory, "log.txt")):
            parent.append("FindMySequence results will appear here soon...")
        else: 
            alog=open(os.path.join(jobDirectory, "log.txt"), "r")
            lines="".join(alog.readlines())
            #lines=alog.readlines()
            alog.close()
            tableFoldmr.append("<pre>%s</pre>" % lines)
        return
