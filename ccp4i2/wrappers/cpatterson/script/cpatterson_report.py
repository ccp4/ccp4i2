from report.CCP4ReportParser import *
import sys

class cpatterson_report(Report):
    TASKNAME = 'cpatterson'
    
    def __init__(self,xmlnode=None,jobInfo={},**kw):
        Report.__init__(self,xmlnode=xmlnode,jobInfo=jobInfo,**kw)
        self.addText(text='The cpatterson job has finished.')
        

