import os
import re
import shutil
import subprocess
from lxml import etree

from core.CCP4PluginScript import CPluginScript
from core import CCP4XtalData
from core import CCP4ErrorHandling
from core import CCP4Utils
from core import CCP4Modules

class buster(CPluginScript):
    TASKMODULE = 'refinement'         # Gui menu location
    TASKTITLE = 'Refinement with Buster (Global Phasing Limited)'        # Short title for Gui
    TASKNAME = 'buster'               # Task name - same as class name
    TASKCOMMAND = 'refine'            # The command to run the executable
    TASKVERSION = 1.0                 # plugin version
    COMTEMPLATE = None                # The program com file template
    COMTEMPLATEFILE = None            # Name of file containing com file template
    WHATNEXT = ['buster', 'prosmart_refmac', 'buccaneer_build_refine_mr', 'modelcraft']
    PERFORMANCECLASS = 'CRefinementPerformance'
    ASYNCHRONOUS = False
    MAINTAINER = 'kyle.stevenson@stfc.ac.uk'
    
    ERROR_CODES = { 101 : {'description' : 'Failed to initialise BUSTER, do you have BUSTER installed & the i2 preferences setup ' \
                           'to point to the correct BUSTER installation folder (or have run the setup script for BUSTER) ?', 
                           'severity':CCP4ErrorHandling.SEVERITY_ERROR } }

    def __init__(self, *args, **kwargs):
        self.seqin = None
        self.hklin = None
        CPluginScript.__init__(self, *args, **kwargs)

    def process(self):
        goodtogo = False
        # Need to first check that Buster is switch on / present. Rem to put Buster in no-Windows list before release.
        bpres_act = shutil.which('refine')
        if bpres_act:
            goodtogo = True
        elif CCP4Modules.PREFERENCES().BUSTERDIR.exists():
            scriplo = os.path.join(CCP4Modules.PREFERENCES().BUSTERDIR.__str__(), 'setup.sh')
            print(scriplo)
            self.source_script(scriplo)
            goodtogo = True
        if goodtogo:
            CPluginScript.process(self)
        else:
            # Failed to find BUSTER. Flag problem & also write advice into stdout.
            self.appendErrorReport(101)
            print("ERROR REPORTED : Failed to find BUSTER installation.\n"
                  "                 If you have installed BUSTER, please ensure you either run the setup script\n" 
                  "                 for BUSTER, before running ccp4i2, or point i2 to the BUSTER installation folder\n"
                  "                 in user preferences.")
            CPluginScript.process(self)

    def source_script(self, script):
        pipe = subprocess.Popen(". %s; env" % script, stdout=subprocess.PIPE, shell=True)
        output = pipe.communicate()[0]
        env = dict((line.decode('utf-8').split("=", 1) for line in output.splitlines()))
        os.environ.update(env)
    
    def processInputFiles(self): 
        self.pdbin = self.container.inputData.XYZIN.fullPath.__str__()
        self.dictin = None
        if self.container.inputData.DICT.isSet():
            self.dictin = self.container.inputData.DICT.fullPath.__str__()
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
        # Run CAD to change FREER flag name to something refine likes
        self.hklin, __, errorb = self.makeHklInput(cols1, extendOutputColnames=True, useInputColnames=True)
        binc = os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir().__str__(), 'bin', 'cad' ))
        self.outfilec = os.path.join(os.path.split(self.hklin)[0], 'rcadout.mtz')
        self.logfc = os.path.join(os.path.split(self.hklin)[0], 'rcadout.log')
        arglist = ['hklin1', self.hklin, 'hklout', self.outfilec]
        # Extra input needed to be read-in by cad
        if self.bIData:
            inputText = "LABIN FILE 1 E1=F_SIGF_F E2=F_SIGF_SIGF E3=F_SIGF_I E4=F_SIGF_SIGI E5=FREERFLAG_FREER\n"
            inputText += ("LABOUT FILE 1 E1=F_SIGF_F E2=F_SIGF_SIGF E3=F_SIGF_I E4=F_SIGF_SIGI E5=FREER")
        elif self.bFData:
            inputText = "LABIN FILE 1 E1=F_SIGF_F E2=F_SIGF_SIGF E3=FREERFLAG_FREER\n"
            inputText += ("LABOUT FILE 1 E1=F_SIGF_F E2=F_SIGF_SIGF E3=FREER")
        # Fire the hkl conversion up.
        pid = CCP4Modules.PROCESSMANAGER().startProcess(binc, arglist, inputText=inputText, logFile=self.logfc)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        extCde = CCP4Modules.PROCESSMANAGER().getJobData(pid, 'exitCode')
        if not(status == 0 and os.path.exists(self.outfilec)):
            return CPluginScript.FAILED
        if errorb.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        # Keep this similar to Refmac5 output 
        pdbfile_bus = os.path.join(self.getWorkDirectory(), "refine.pdb")
        cif_coord_bus = os.path.join(self.getWorkDirectory(), "BUSTER_model.cif")
        cif_refln_bus = os.path.join(self.getWorkDirectory(), "BUSTER_refln.cif")
        if os.path.exists(cif_refln_bus):
            self.container.outputData.REFLFN_CIFFILE.set(cif_refln_bus)
            self.container.outputData.REFLFN_CIFFILE.annotation.set('MMCIF reflections from refinement ("BUSTER_refln.cif")')
        if os.path.exists(cif_coord_bus):
            self.container.outputData.CIFFILE.set(cif_coord_bus)
            self.container.outputData.CIFFILE.annotation.set('MMCIF model from refinement ("BUSTER_model.cif")')
        if os.path.exists(pdbfile_bus):
            self.container.outputData.XYZOUT = pdbfile_bus
        self.container.outputData.XYZOUT.annotation = 'Model from refinement'
        # MTZ i2 int conv.
        mtzfile_bus = os.path.join(self.getWorkDirectory(), "refine.mtz")
        self.container.outputData.ABCDOUT.annotation = 'Calculated phases from refinement'
        #self.container.outputData.ABCDOUT.contentFlag = CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL
        self.container.outputData.FPHIOUT.annotation = 'Weighted map from refinement'
        self.container.outputData.DIFFPHIOUT.annotation = 'Weighted difference map from refinement'
        outFiles = ['FPHIOUT', 'DIFFPHIOUT', 'ABCDOUT']
        outCols =  ['2FOFCWT,PH2FOFCWT', 'FOFCWT,PHFOFCWT', 'HLA,HLB,HLC,HLD']
        rep = self.splitHklout(miniMtzsOut=outFiles, programColumnNames=outCols, infile=mtzfile_bus)
        if rep.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        # Extract what is needed for the i2 Buster report.
        blfilep = self.makeFileName('LOG')
        blfile = open(blfilep, 'r')
        bltxt = blfile.read()
        gop = re.findall(r"(best refinement in BUSTER reached\D*)(\d*\.?\d+)/(\d*\.?\d+)", bltxt)
        linee = bltxt.splitlines()
        blfile.close()
        # Final R/Rfree from log
        rrfr = [gop[0][1], gop[0][2]]
        # Set Perfm. Output
        self.container.outputData.PERFORMANCEINDICATOR.RFactor.set(str(rrfr[0]))
        self.container.outputData.PERFORMANCEINDICATOR.RFree.set(str(rrfr[1]))
        # Extract graphs from logfile
        graphf = False
        allcyc = []
        for linea in linee:
            if graphf:
                numin = re.findall(r"([+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+))|([+\-]?\d*\.?\d+)|(-)", linea) # want 0, 3 , 4 (cyc, r, rf)
                if numin:
                    try:
                        llg = float("".join(numin[7]))
                        llgf = float("".join(numin[8]))
                    except:
                        llg = 0
                        llgf = 0
                    cyc = [int("".join(numin[0])), float("".join(numin[3])), float("".join(numin[4])),
                           llg, llgf, float("".join(numin[10]))*100.0, float("".join(numin[11]))]
                    allcyc.append(cyc)
                else:
                    graphf = False
            fig = re.search(r"Ncyc\s+Total\s+Grms\s+Rfact\s+Rfree", linea)
            if fig:
                graphf = True
        # xml (should not be necessary for graphs but images will not load)
        rootNode = etree.Element("Buster")
        xmlRI = etree.SubElement(rootNode, "RunInfo")
        xmlbcyc = etree.SubElement(xmlRI, "Best")
        etree.SubElement(xmlbcyc, "R").text = str(rrfr[0])
        etree.SubElement(xmlbcyc, "RFree").text = str(rrfr[1])
        # Graphs into xml format recogn. by i2
        for ij, cycle in enumerate(allcyc):
            xmlcyc = etree.SubElement(xmlRI, "Cycle")
            etree.SubElement(xmlcyc, "NCycle").text = str(ij)
            etree.SubElement(xmlcyc, "RFact").text = str(cycle[1])
            etree.SubElement(xmlcyc, "RFree").text = str(cycle[2])
            etree.SubElement(xmlcyc, "LLG").text = str(cycle[3])
            etree.SubElement(xmlcyc, "LLGF").text = str(cycle[4])
            etree.SubElement(xmlcyc, "RMSB").text = str(cycle[5])
            etree.SubElement(xmlcyc, "RMSA").text = str(cycle[6])
        # Save xml
        xmlfile = open(self.xmlout, 'wb')
        xmlString= etree.tostring(rootNode, pretty_print=True)
        xmlfile.write(xmlString)
        xmlfile.close()
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, container=None):
        self.appendCommandLine("-p")
        self.appendCommandLine(str(self.pdbin))
        self.appendCommandLine("-m")
        if self.dictin:
            self.appendCommandLine("-l")
            self.appendCommandLine(str(self.dictin))
        self.appendCommandLine(str(self.outfilec))
        self.appendCommandLine("-nbig")
        self.appendCommandLine(str(self.container.inputParameters.NBCYCLES))
        self.appendCommandLine("-nsmall")
        self.appendCommandLine(str(self.container.inputParameters.NSCYCLES))
        if self.container.inputParameters.AUTO_NCS:
            self.appendCommandLine("-autoncs")
        if self.container.inputParameters.RBR:
            self.appendCommandLine("-RB")
        if self.container.inputParameters.AUTO_NCS:
            self.appendCommandLine("-TLS")
        # Water treatment
        wtrt = str(self.container.inputParameters.WAT)
        if wtrt == "ON":
            self.appendCommandLine("-WAT")
        elif wtrt == "MAN":
            self.appendCommandLine("-WAT")
            self.appendCommandLine(str(self.container.inputParameters.WATCYC))
        elif wtrt == "OFF":
            self.appendCommandLine("-noWAT")
        self.xmlout = self.makeFileName('PROGRAMXML')
        return CPluginScript.SUCCEEDED
