from ccp4i2.wrappers.ShelxCDE.script import ShelxCDEBaseReport


class ShelxCE_report(ShelxCDEBaseReport.ShelxCDEBaseReport):
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
