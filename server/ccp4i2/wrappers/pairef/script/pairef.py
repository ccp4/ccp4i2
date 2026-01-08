import os

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript


class pairef(CPluginScript):
    TASKMODULE = 'refinement'         # Gui menu location
    TASKTITLE = 'Pairef'        # Short title for Gui
    TASKNAME = 'pairef'               # Task name - same as class name
    TASKCOMMAND = 'pairef'            # The command to run the executable
    TASKVERSION = 1.1                 # plugin version
    COMTEMPLATE = None                # The program com file template
    COMTEMPLATEFILE = None            # Name of file containing com file template
    PERFORMANCECLASS = 'CPairefPerformance'
    MAINTAINER = 'stuart.mcnicholas@york.ac.uk'
    
    ERROR_CODES = { 101 : {'description' : 'Blank for now, may need this ',
                           'severity':CCP4ErrorHandling.SEVERITY_ERROR } }

    def __init__(self, *args, **kwargs):
        self.seqin = None
        self.hklin = None
        self.pdbin = None
        self.tlsin = None
        self.unmer = None
        self.comfile = None
        CPluginScript.__init__(self, *args, **kwargs)

    def processInputFiles(self): 
        self.pdbin = self.container.inputData.XYZIN.fullPath.__str__()
        if self.container.inputData.TLSIN.isSet():
            self.tlsin = self.container.inputData.TLSIN.fullPath.__str__()
        if self.container.inputData.UNMERGED.isSet():
            self.unmer = self.container.inputData.UNMERGED.fullPath.__str__()
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
        if self.container.inputData.FREERFLAG.isSet():
            cols1.append(['FREERFLAG', None])
        self.hklin, __, errorb = self.makeHklInput(cols1, extendOutputColnames=True, useInputColnames=True)
        if errorb.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, container=None):
        self.appendCommandLine("--XYZIN")
        self.appendCommandLine(str(self.pdbin))
        self.appendCommandLine("--HKLIN")
        self.appendCommandLine(str(self.hklin))
        if self.container.inputParameters.FIXED_TLS:
            self.appendCommandLine("--TLSIN") # See also TLSCYC below
            self.appendCommandLine(str(self.tlsin))
        if self.unmer:
            self.appendCommandLine("--unmerged")
            self.appendCommandLine(str(self.unmer))
        if self.container.inputParameters.INIRES > 0.0:
            self.appendCommandLine("-i")
            self.appendCommandLine(str(self.container.inputParameters.INIRES))
        # Pre-refinement
        if self.container.inputParameters.USE_PREREF:
            self.appendCommandLine("--prerefinement-ncyc")
            self.appendCommandLine(str(self.container.inputParameters.NPRECYCLES))
            if self.container.inputParameters.USE_SHAKE:
                self.appendCommandLine("--prerefinement-shake-sites")
                self.appendCommandLine(str(self.container.inputParameters.SHAKE))
            if self.container.inputParameters.RESETBFAC:
                self.appendCommandLine("--prerefinement-reset-bfactor")
        # Define shells
        if str(self.container.inputParameters.SH_TYPE) == 'manual':
            self.appendCommandLine("-r")
            self.appendCommandLine(str(self.container.inputParameters.MANSHELL))
        elif str(self.container.inputParameters.SH_TYPE) == 'semi':
            self.appendCommandLine("-n")
            self.appendCommandLine(str(self.container.inputParameters.NSHELL))
            self.appendCommandLine("-s")
            self.appendCommandLine(str(self.container.inputParameters.WSHELL))
        if not self.container.inputParameters.AUTO_WGT:
            self.appendCommandLine("-w")
            self.appendCommandLine(str(self.container.inputParameters.WGT_TRM))
        if self.container.inputParameters.FIXED_TLS:
            self.appendCommandLine("--tls-ncyc")
            self.appendCommandLine(str(self.container.inputParameters.TLSCYC))
        if self.container.inputParameters.COMPLETE:
            self.appendCommandLine("--complete")
        self.appendCommandLine("--ncyc")
        self.appendCommandLine(str(self.container.inputParameters.NCYCLES))

        if self.container.inputData.DICT.isSet():
            self.dictin = self.container.inputData.DICT.fullPath.__str__()
            self.appendCommandLine("--libin")
            self.appendCommandLine(str(self.dictin))
        if self.container.inputParameters.REFMAC_KEYWORD_FILE.isSet():
            self.comfile = self.container.inputParameters.REFMAC_KEYWORD_FILE.fullPath.__str__()
            self.appendCommandLine("--comfile")
            self.appendCommandLine(str(self.comfile))

        self.xmlout = self.makeFileName('PROGRAMXML')
        rootNode = etree.Element("Pairef")
        # Save xml
        xmlfile = open(self.xmlout, 'wb')
        xmlString= etree.tostring(rootNode, pretty_print=True)
        xmlfile.write(xmlString)
        xmlfile.close()
        
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        cutoff = '0.0'
        cfn = os.path.join(self.getWorkDirectory(), "pairef_project" ,"PAIREF_cutoff.txt")
        with open(cfn) as cf:
            cutoff = cf.read().strip()
        self.container.outputData.PERFORMANCEINDICATOR.cutoff.set(cutoff)

        rootNode = etree.Element("Pairef")
        etree.SubElement(rootNode, "Cutoff").text = cutoff

        # Save xml
        xmlfile = open(self.xmlout, 'wb')
        xmlString= etree.tostring(rootNode, pretty_print=True)
        xmlfile.write(xmlString)
        xmlfile.close()
        return CPluginScript.SUCCEEDED
