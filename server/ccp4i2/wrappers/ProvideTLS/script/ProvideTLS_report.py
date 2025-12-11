from ccp4i2.report.CCP4ReportParser import *
import sys

class ProvideTLS_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ProvideTLS'
    RUNNING = False
    def __init__(self,xmlnode=None,jobInfo={},jobStatus=None,**kw):
        Report. __init__(self,xmlnode=xmlnode,jobInfo=jobInfo, **kw)
        
        sequenceFold = self.addFold(label="TLS provided:", initiallyOpen=True)
        for sequenceNode in self.xmlnode.findall('.//TLSProvided'):
            preDiv = sequenceFold.addPre(text = sequenceNode.text)
