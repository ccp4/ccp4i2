from report.CCP4ReportParser import *
import sys
#from lxml import etree
import xml.etree.ElementTree as etree
import math
from ccp4i2.wrappers.refmac_i2.script.refmac_report import refmac_report
class PrepareDeposit_report(Report):
    # Specify which gui task and/or pluginscript this applies to
    TASKNAME = 'PrepareDeposit'
    RUNNING = True

    def __init__(self, *args, **kws):
        '''
        kws['additionalJsFiles'] = ['jsme/jsme.nocache.js']
        kws['additionalScript'] = 'function jsmeOnLoad() { jsmeApplet = new JSApplet.JSME("jsme_container", "380px", "340px");'
        '''
        Report.__init__(self, *args, **kws)
        if self.jobStatus is None or self.jobStatus.lower() == 'nooutput': return
        self.addDiv(style='clear:both;')
        self.defaultReport()

    def defaultReport(self, parent=None):
        if parent is None: parent=self
        alignChainNodes = self.xmlnode.findall('.//AlignChain')
        for alignChainNode in alignChainNodes:
            chainIdNode = alignChainNode.findall('ChainId')[0]
            alignChainFold = parent.addFold(label="Alignments for chain "+chainIdNode.text)
            alignChainFold.addPre(text = alignChainNode.findall('Alignment')[0].text)
        if len(self.xmlnode.findall('REFMAC'))>0:
            rmReport = refmac_report(xmlnode=self.xmlnode.findall('REFMAC')[0],jobStatus='nooutput')
            rmReport.addSummary(parent=self)

