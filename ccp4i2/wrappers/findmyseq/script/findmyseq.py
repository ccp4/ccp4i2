import os
import json
from lxml import etree
from core.CCP4PluginScript import CPluginScript
from core import CCP4XtalData
from core import CCP4ErrorHandling
from core import CCP4Utils
from core import CCP4Modules

class findmyseq(CPluginScript):
    TASKMODULE = 'bioinformatics'         # Gui location
    TASKTITLE = 'Find My Sequence'    # Short title for Gui
    TASKNAME = 'findmyseq'            # Task name - same as class name
    TASKCOMMAND = 'findmysequence'    # executable
    TASKVERSION = 1.0                 # plugin version
    COMTEMPLATE = None                # The program com file template
    COMTEMPLATEFILE = None            # Name of file containing com file template
    WHATNEXT = ['modelcraft']         # after ?
    PERFORMANCECLASS = 'CRefinementPerformance'
    ASYNCHRONOUS = False
    MAINTAINER = 'kyle.stevenson@stfc.ac.uk'
    
    ERROR_CODES = { 101 : {'description' : 'Blank ' \
                           'In case needed (Prob not needed here)', 
                           'severity':CCP4ErrorHandling.SEVERITY_ERROR } }

    def __init__(self, *args, **kwargs):
        self.outjfile = None
        self.mapin = None
        self.pdbin = None # Purely for the main chain (side chain not used)
        self.lseqdb = None
        self.seqout = None
        CPluginScript.__init__(self, *args, **kwargs)

    def process(self):
        CPluginScript.process(self)

    def processInputFiles(self):
        self.pdbin = self.container.inputData.XYZIN.fullPath.__str__()
        if self.container.inputData.LSEQDB.isSet():
            self.lseqdb = self.container.inputData.LSEQDB.fullPath.__str__()
        if self.container.inputData.XYZIN.isSelectionSet():
            self.pdbin = os.path.join(self.workDirectory,'XYZIN_selected_atoms.pdb')
            rv = self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.pdbin)
        else:
            self.pdbin = str( self.container.inputData.XYZIN)
        self.mapin = self.container.inputData.FPHI.fullPath.__str__()
        cols1 = []
        self.container.inputData.F_SIGF.loadFile()
        self.container.inputData.F_SIGF.setContentFlag()
        self.bFData = self.container.inputData.F_SIGF.contentFlag == 2 or self.container.inputData.F_SIGF.contentFlag == 4
        self.bIData = self.container.inputData.F_SIGF.contentFlag == 1 or self.container.inputData.F_SIGF.contentFlag == 3
        if self.bIData:
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN])
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN])
        if self.bFData:
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN])
        cols1.append(['FPHI', None])
        self.hklin, __, errorb = self.makeHklInput(cols1, extendOutputColnames=True, useInputColnames=True)
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        # LoadJsonOut Later (I'm still missing the json from trunk version of ccp4)
        self.outjfile = os.path.join(self.getWorkDirectory(), "fmsresults.json")
        self.loadJsonOut()
        if self.seqout:
            self.container.outputData.SEQOUT = self.seqout
            self.container.outputData.SEQOUT.annotation = "Best sequence file from FindMySequence"
        # xml (should not be necessary for graphs but images will not load)
        rootNode = etree.Element("FindMySeq")
        #xmlRI = etree.SubElement(rootNode, "RunInfo")
        # Save xml
        xmlfile = open(self.xmlout, 'wb')
        xmlString= etree.tostring(rootNode, pretty_print=True)
        xmlfile.write(xmlString)
        xmlfile.close()
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, container=None):
        self.appendCommandLine("--mtzin")
        self.appendCommandLine(str(self.hklin))
        self.appendCommandLine("--labin")
        self.appendCommandLine("F,PHI")
        self.appendCommandLine("--modelin")
        self.appendCommandLine(str(self.pdbin))
        if self.lseqdb:
            self.appendCommandLine("--db")
            self.appendCommandLine(str(self.lseqdb))
        self.appendCommandLine("--jsonout")
        self.appendCommandLine("fmsresults.json")
        # The xml isn't needed as such, but unfortunately the original code fails if not there...
        self.xmlout = self.makeFileName('PROGRAMXML')
        return CPluginScript.SUCCEEDED

    def loadJsonOut(self):
        # Load json file & write output sequence
        jfi = open(self.outjfile)
        data = json.load(jfi)
        jfi.close()
        # loop over FindMySeq results.
        bestseq = None
        bestseqe = 101.0
        for key, item in data.items():
            if float(item.get("evalue")) < bestseqe:
                bestseq = str(item.get("sequence_id"))
                bestseqe = float(item.get("evalue"))
                bestkey = key
        if bestseq:
            seqid = "".join(("> ", str(data.get(bestkey).get("sequence_id")), "\n"))
            seqtrs = str(data.get(bestkey).get("sequence"))
            # Save the sequence data to output file.
            self.seqout = os.path.join(self.getWorkDirectory(), 'seqout.fasta')
            sfo = open(self.seqout, 'w')
            sfo.write(seqid)
            sfo.write(seqtrs)
            sfo.close()

