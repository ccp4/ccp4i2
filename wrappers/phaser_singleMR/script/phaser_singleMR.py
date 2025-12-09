from __future__ import print_function
#=======================================================================================
#
#    shelxeMR.py : phaserSingleMR(CPluginScript)
#    
#    Author  : Kyle Stevenson,STFC
#    Created : 10th July 2018, KJS
#
#    Complementary Gui Class for running phaser in Single Atom mode.
#    Handles the command and script processing.
#
#=======================================================================================

import os
import sys
from lxml import etree
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4XtalData
from ccp4i2.core import CCP4ErrorHandling


class phaser_singleMR(CPluginScript):

    TASKMODULE = 'molecular_replacement'      # Gui menu location
    TASKTITLE = 'Single Atom MR'              # Short title for Gui
    TASKNAME = 'phaser_singleMR'              # Task name - same as class name
    TASKCOMMAND = 'phaser'                    # The command to run the executable
    TASKVERSION = 1.0                         # plugin version
    COMTEMPLATE = None
    COMTEMPLATEFILE = None
    PERFORMANCECLASS = 'CExpPhasPerformance'
    ASYNCHRONOUS = True
    MAINTAINER = 'Kyle.Stevenson@stfc.ac.uk'   #  I think this is mostly ok now

    def __init__(self, *args, **kwargs):
        self.hklin1 = None
        self.bFData = None
        self.bIData = None
        self.xmlout = None
        CPluginScript.__init__(self, *args, **kwargs)

    def process(self):
        CPluginScript.process(self)

    def processInputFiles(self):
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
        self.hklin1, __, error1 = self.makeHklInput(cols1, extendOutputColnames=True, useInputColnames=True)
        # Need to join up the previous i2-only mini-mtz's to get back the prior mtz file ...
        self.seqin = os.path.join(self.workDirectory, 'input_seq.fasta')
        self.container.inputData.ASUFILE.writeFasta(self.seqin)
        if error1.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        num_sol = 1
        for i in range(1, num_sol + 1):
            xyzout = os.path.join(self.getWorkDirectory(), "SingleMR." + str(i) + ".pdb")
            if os.path.exists(xyzout):
                self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                self.container.outputData.XYZOUT[-1].setFullPath(xyzout)
                self.container.outputData.XYZOUT[-1].annotation = 'Positioned coordinates for solution ' + str(i)
            else:
                self.appendErrorReport(201, xyzout)
                return CPluginScript.FAILED
            hklout = os.path.join(self.getWorkDirectory(), "SingleMR." + str(i) + ".mtz")
            if os.path.exists(hklout):
                self.container.outputData.HKLOUT.append(self.container.outputData.HKLOUT.makeItem())
                self.container.outputData.HKLOUT[-1].setFullPath(hklout)
            else:
                self.appendErrorReport(201, hklout)
                return CPluginScript.FAILED
        self.splitHkloutList(miniMtzsOut=['MAPOUT', 'ABCDOUT'],
                             programColumnNames=['FWT,PHWT', 'HLA,HLB,HLC,HLD'],
                             outputBaseName=['MAPOUT', 'ABCDOUT'],
                             outputContentFlags=[1, CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL],
                             infileList=self.container.outputData.HKLOUT)
        for indx in range(len(self.container.outputData.MAPOUT)):
            self.container.outputData.MAPOUT[indx].annotation = 'Map for solution ' + str(indx + 1)
            self.container.outputData.MAPOUT[indx].contentFlag = 1
            self.container.outputData.MAPOUT[indx].subType = 1
            self.container.outputData.ABCDOUT[indx].annotation = 'H-L Co-efficients' + str(indx + 1)
        self.parseLogfile()
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, container=None):
        # nb. the -xml flag doesn't work in phaser & the html flag does not work in smartie.
        self.appendCommandScript("TITLE SingleMR")
        self.appendCommandScript("MODE MR_ATOM")
        self.appendCommandScript("ROOT SingleMR")
        self.appendCommandScript("HKLIN \"%s\" "%str(self.hklin1))
        if self.bIData:
            self.appendCommandScript("LABIN I=F_SIGF_I SIGI=F_SIGF_SIGI")
        if self.bFData:
            self.appendCommandScript("LABIN F=F_SIGF_F SIGF=F_SIGF_SIGF")  # Should be I or F..
        if not self.bIData and not self.bFData:
            print("phaser_singleMR : Something has gone horribly wrong with the input") # Should really trigger an exception here...
        if self.container.inputData.RESOL_ON:
            self.appendCommandScript("RESOLUTION HIGH %s LOW %s"%(str(self.container.inputData.RESOL_HI),
                                                                  str(self.container.inputData.RESOL_LO)))
        # self.appendCommandScript("RESOLUTION")  # The resolution.
        self.appendCommandScript("SGALTERNATIVE SELECT NONE")  # Need to fix the space group so user can run it manually.
        # Define the search
        self.appendCommandScript("ENSEMBLE SINGLEATOM ATOM %s"%str(self.container.inputData.SINGLE_ATOM_TYPE)) # Needs the atom & number
        self.appendCommandScript("SEARCH ENSEMBLE SINGLEATOM NUMBER %s"%str(self.container.inputData.SINGLE_ATOM_NUM))
        self.appendCommandScript("LLGCOMPLETE SCATTERER %s"%str(self.container.inputData.LLG_COMPL_ATOMTYP)) # Needs the LLG atom
        if self.container.inputData.LLG_COMPL_ATMSEP_ON:
            self.appendCommandScript("LLGCOMPLETE CLASH %s"%str(self.container.inputData.LLG_COMPL_ATMSEP))
        if self.container.inputData.LLG_COMPL_SIGCO_ON:
            self.appendCommandScript("LLGCOMPLETE SIGMA %s"%str(self.container.inputData.LLG_COMPL_SIGCO))
        if self.container.inputData.LLG_COMPL_MAXCYC_ON:
            self.appendCommandScript("LLGCOMPLETE NCYC %s"%str(self.container.inputData.LLG_COMPL_MAXCYC))
        # Need the ASU Unit to be defined
        self.appendCommandScript("COMPOSITION BY ASU")   # Obviously can be something else.
        self.appendCommandScript("COMPOSITION PROTEIN SEQ \"%s\" NUMBER 1"%str(self.seqin))
        # Expert switches
        if self.container.inputDataExp.EXP_PACKCT_ON:
            if str(self.container.inputDataExp.EXP_PACKCT_TYPE) != "ALL":
                self.appendCommandScript("PACK SELECT %s CUTOFF %s"%(str(self.container.inputDataExp.EXP_PACKCT_TYPE),
                                                                     str(self.container.inputDataExp.EXP_PACKCT_AMT)))
                if self.container.inputDataExp.EXP_QKPACK_ON:
                    self.appendCommandScript("PACK QUICK ON")
                else:
                    self.appendCommandScript("PACK QUICK OFF")
            else:
                self.appendCommandScript("PACK SELECT ALL")
        if self.container.inputDataExp.EXP_TRAN_SRCHPK_ON:
            if str(self.container.inputDataExp.EXP_TRAN_SRCHPK_TYPE) != "ALL":
                self.appendCommandScript("PEAKS TRA SELECT %s CUTOFF %s"%(str(self.container.inputDataExp.EXP_TRAN_SRCHPK_TYPE),
                                                                          str(self.container.inputDataExp.EXP_TRAN_SRCHPK_AMT)))
            else:
                self.appendCommandScript("PEAKS TRA SELECT ALL")
        # PURGE TRA ENABLE ON  & then PURGE TRA PERCENT 75  & also   PURGE TRA NUMBER 52
        if self.container.inputDataExp.EXP_PURGE_TRANPK_ON:
            self.appendCommandScript("PURGE TRA ENABLE ON")
            self.appendCommandScript("PURGE TRA PERCENT %s"%str(self.container.inputDataExp.EXP_PURGE_TRANPK_PER))
            self.appendCommandScript("PURGE TRA NUMBER %s"%str(self.container.inputDataExp.EXP_PURGE_TRANPK_NUM))
        # RESOLUTION AUTO HIGH 0.95 LOW 11.0
        if self.container.inputDataExp.EXP_RESRAN_HRREFINE_ON:
            self.appendCommandScript("RESOLUTION AUTO HIGH %s LOW %s"%(str(self.container.inputDataExp.EXP_RESRAN_HRREFINE_HI),
                                                                       str(self.container.inputDataExp.EXP_RESRAN_HRREFINE_LO)))
        self.xmlout = self.makeFileName('PROGRAMXML')
        return CPluginScript.SUCCEEDED

    def parseLogfile(self):
        logfile = self.makeFileName('LOG')
        # Load ccp4 smartie and qtr code
        from ccp4i2.core import CCP4Utils
        smpth = os.path.join(CCP4Utils.getCCP4Dir(), 'share', 'smartie')
        sys.path.append(smpth)
        import smartie
        import qtrgeneric
        from qtrgeneric import CLReader, LogConverter
        # Use ccp4 automated output. Note this is replacing I2's internal xml file handling.
        infi = CLReader({'REP_XML': 'program.xml', 'REP_XRT': 'program.xrt', 'LOGFILE': logfile})
        smin = smartie.parselog(logfile)
        lconv = LogConverter()
        lconv.convert(smin, infi)
        lconv.xmltree.write(os.path.join(self.getWorkDirectory(), infi.xml))
        xrto = open(os.path.join(self.getWorkDirectory(), infi.xrt), "w")
        from xml.etree import ElementTree as ET
        xrto.write((ET.tostring(lconv.xrttree.getroot()).decode()).replace("ns0", "xrt"))
        xrto.close()
        return
