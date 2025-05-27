import os

from lxml import etree

from ....core.CCP4PluginScript import CPluginScript
from ....pimple import MGQTmatplotlib
from ....smartie import smartie


class cphasematch(CPluginScript):

    TASKMODULE = 'expt_data_utility'
    TASKTITLE = 'Match and analyse phases to reference set'
    TASKNAME = 'cphasematch'
    TASKCOMMAND = 'cphasematch'
    DESCRIPTION = '''Compare phases from different sources (with option change of origin/hand)'''
    TASKVERSION= 0.0

    def processInputFiles ( self ):
        from ....core import CCP4XtalData
        inp = self.container.inputData
        colgrps = [ ['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN],
                    'ABCD1', 'ABCD2' ] 
        self.hklin, columns, error = self.makeHklin0 ( colgrps )


    def makeCommandAndScript(self):
        inp = self.container.inputData
        out = self.container.outputData

        self.hklout = os.path.join(self.workDirectory,"hklout.mtz")

        self.appendCommandLine([ '-stdin' ])

        self.appendCommandScript( '-mtzin ' + self.hklin )
        self.appendCommandScript( '-mtzout ' + self.hklout )
        self.appendCommandScript ([ '-colin-fo F_SIGF_F,F_SIGF_SIGF' ])
        if self.container.inputData.ABCD1.contentFlag == self.container.inputData.ABCD1.CONTENT_FLAG_HL:
          self.appendCommandScript ([ '-colin-hl-1 ABCD1_HLA,ABCD1_HLB,ABCD1_HLC,ABCD1_HLD' ])
        else:
          self.appendCommandScript ([ '-colin-phifom-1 ABCD1_PHI,ABCD1_FOM' ])
        if self.container.inputData.ABCD2.contentFlag == self.container.inputData.ABCD2.CONTENT_FLAG_HL:
          self.appendCommandScript ([ '-colin-hl-2 ABCD2_HLA,ABCD2_HLB,ABCD2_HLC,ABCD2_HLD' ])
        else:
          self.appendCommandScript ([ '-colin-phifom-2 ABCD2_PHI,ABCD2_FOM' ])
        self.appendCommandScript( '-colout i2' )

        return CPluginScript.SUCCEEDED


    def processOutputFiles(self):
        error = self.splitHklout( [ 'ABCDOUT' ], [ 'i2.ABCD.A,i2.ABCD.B,i2.ABCD.C,i2.ABCD.D' ] )
        self.container.outputData.ABCDOUT.annotation = 'shifted ABCD'

        lines = open(self.makeFileName('LOG')).readlines()
        for i in range(len(lines)):
          if "Overall" in lines[i]:
            dphi = lines[i+2].split()[3]
            wdphi1 = lines[i+2].split()[4]
            wdphi2 = lines[i+2].split()[5]
            fcorr = lines[i+2].split()[6]
            ecorr = lines[i+2].split()[7]
            self.container.outputData.PERFORMANCE.phaseError = float(dphi)
            self.container.outputData.PERFORMANCE.weightedPhaseError = float(wdphi1)
            self.container.outputData.PERFORMANCE.reflectionCorrelation = float(ecorr)

            self.xmlnode = ET.Element('cphasematch')
            ET.SubElement(self.xmlnode,'phaseError').text = dphi
            ET.SubElement(self.xmlnode,'weightedPhaseError1').text = wdphi1
            ET.SubElement(self.xmlnode,'weightedPhaseError2').text = wdphi2
            ET.SubElement(self.xmlnode,'reflectionFCorrelation').text = fcorr
            ET.SubElement(self.xmlnode,'reflectionECorrelation').text = ecorr

            smartieNode = ET.SubElement(self.xmlnode,'SmartieGraphs')
            self.scrapeSmartieGraphs(smartieNode)

            ET.ElementTree(self.xmlnode).write(self.makeFileName('PROGRAMXML'))

            return CPluginScript.SUCCEEDED

        return CPluginScript.FAILED


    def scrapeSmartieGraphs(self, smartieNode):
        logfile = smartie.parselog(self.makeFileName('LOG'))
        for smartieTable in logfile.tables():
            if smartieTable.ngraphs() > 0:
                tableetree = MGQTmatplotlib.CCP4LogToEtree(smartieTable.rawtable())
                smartieNode.append(tableetree)
        return
