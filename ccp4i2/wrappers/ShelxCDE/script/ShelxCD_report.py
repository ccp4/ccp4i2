from ccp4i2.report.CCP4ReportParser import *
import sys
import math
from ccp4i2.wrappers.ShelxCDE.script.ShelxCDEBaseReport import ShelxCDEBaseReport

class ShelxCD_report(ShelxCDEBaseReport):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ShelxCD'
    RUNNING = True
    SEPARATEDATA=True
    
    def __init__(self,*args,**kws):
        super(ShelxCD_report,self).__init__(*args,**kws)
    
    def defaultReport(self, parent=None):
        if parent is None: parent = self
        
        datasetNodes = self.xmlnode.findall('.//Dataset')
        if len(self.xmlnode.findall('.//Shelxd'))==0:
            self.shelXCReport(parent, initiallyOpen=True )
        else:
            self.shelXCReport(parent, initiallyOpen=False )
            self.shelXDReport(parent, initiallyOpen=True)
