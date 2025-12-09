from ccp4i2.report.CCP4ReportParser import *
import sys
import math
from ccp4i2.pipelines.pisapipe.wrappers.pisa_xml.script.pisa_xml_report import pisa_xml_report

class pisapipe_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'pisapipe'
    RUNNING = True
    
    def __init__(self,*args,**kws):
        Report.__init__(self, *args, **kws)
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.defaultReport()
    
    def defaultReport(self,parent=None):
        if parent is None: parent=self
        pisa_xmlNode = self.xmlnode.findall('pisa_xml')[0]
        if pisa_xmlNode is not None:
            my_pisa_xml_report = pisa_xml_report(xmlnode=pisa_xmlNode)
            my_pisa_xml_report.defaultReport(parent=self)

