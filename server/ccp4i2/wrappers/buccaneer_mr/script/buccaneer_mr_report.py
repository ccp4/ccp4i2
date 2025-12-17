
from ccp4i2.report import Report


class buccaneer_mr_report(Report):
    TASKNAME = 'buccaneer_mr'
    RUNNING = False

    def __init__(self, *args, **kws):
        Report.__init__(self, *args, **kws)
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None: parent = self
        self.addDiv(style='clear:both;')
        self.summaryTable(self)

    def summaryTable(self, parent=None):
        try:
          cres = float(self.xmlnode.findall('.//BuccaneerResult/Cycles/Cycle/CompletenessByResiduesBuilt')[-1].text)
          cchn = float(self.xmlnode.findall('.//BuccaneerResult/Cycles/Cycle/CompletenessByChainsBuilt')[-1].text)
          frgb = float(self.xmlnode.findall('.//BuccaneerResult/Cycles/Cycle/FragmentsBuilt')[-1].text)
          chnb = float(self.xmlnode.findall('.//BuccaneerResult/Cycles/Cycle/ChainsBuilt')[-1].text)
          resb = float(self.xmlnode.findall('.//BuccaneerResult/Cycles/Cycle/ResiduesBuilt')[-1].text)
          ress = float(self.xmlnode.findall('.//BuccaneerResult/Cycles/Cycle/ResiduesSequenced')[-1].text)
          parent.append( "<p>%d residues were built in %d fragments. Of these, %d residues were assigned to the sequence.</p><p>The number of chains is estimated to be %d. Of these chains, %5.1f%% of the residues have been built.<br/>Of the residues that were built, %5.1f%% were assigned to a chain.</p>"%(resb,frgb,ress,chnb,100*cchn,100*cres) )
        except:
          parent.append( "<p>Encountered problems trying to interpret your latest buccaneer run.</p>" )

        progressTable = parent.addTable(title='Progress through buccaneeer cycle',select='//BuccaneerResult/Cycles/Cycle')
        nameMap = [('Number','Cycle'),('CompletenessByResiduesBuilt','Completeness<br/>by res.'),('CompletenessByChainsBuilt','Completeness<br/>by chain'),('ChainsBuilt','Chains<br/>built'),('FragmentsBuilt','Fragments<br/>built'),('ResiduesUnique','Unique<br/>residues'),('ResiduesBuilt','Residues<br/>built'),('ResiduesSequenced','Residues<br/>sequenced'),('ResiduesLongestFragment','Longest<br/>fragment')]
        for key,value in nameMap:
            progressTable.addData(title = value, select=key)
