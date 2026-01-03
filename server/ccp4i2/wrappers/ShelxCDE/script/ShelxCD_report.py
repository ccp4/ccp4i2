from ccp4i2.wrappers.ShelxCDE.script.ShelxCDEBaseReport import ShelxCDEBaseReport


class ShelxCD_report(ShelxCDEBaseReport):
    TASKNAME = 'ShelxCD'
    RUNNING = True
    SEPARATEDATA=True

    def defaultReport(self, parent=None):
        if parent is None: parent = self
        
        datasetNodes = self.xmlnode.findall('.//Dataset')
        if len(self.xmlnode.findall('.//Shelxd'))==0:
            self.shelXCReport(parent, initiallyOpen=True )
        else:
            self.shelXCReport(parent, initiallyOpen=False )
            self.shelXDReport(parent, initiallyOpen=True)
