from report.CCP4ReportParser import *
import sys
import math
from wrappers.ShelxCDE.script import ShelxCDEBaseReport

class ShelxCE_report(ShelxCDEBaseReport.ShelxCDEBaseReport):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'ShelxCE'
    RUNNING = True
    SEPARATEDATA=True
    
    def defaultReport(self, parent=None):
        if parent is None: parent = self
        
        if len(self.xmlnode.findall('Shelxe'))==0:
            self.shelXCReport(parent, initiallyOpen=True )
        else:
            self.shelXCReport(parent, initiallyOpen=False )
            self.shelXEReport(parent, initiallyOpen=True)
        
