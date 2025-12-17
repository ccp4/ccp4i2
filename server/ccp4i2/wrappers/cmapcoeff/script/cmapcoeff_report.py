import sys

from ccp4i2.report import Report


class cmapcoeff_report(Report):
    TASKNAME = 'cmapcoeff'
    
    def __init__(self,xmlnode=None,jobInfo={},**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)

        results = self.addResults()
        results.append ( 'Please find below the output files. If you want to do a peak search, you can select <i>Manual model rebuilding</i>. ' ) 
