import os
import re
import shutil

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling, CCP4Modules, CCP4Utils, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript


class shelxeMR(CPluginScript):
    TASKMODULE = 'model_building'      # Gui menu location
    TASKTITLE = 'ShelxeMR'             # Short title for Gui
    TASKNAME = 'shelxeMR'              # Task name - same as class name
    TASKCOMMAND = 'shelxe'             # The command to run the executable
    TASKVERSION = 1.0                  # plugin version
    COMTEMPLATE = None                 # The program com file template
    COMTEMPLATEFILE = None             # Name of file containing com file template
    PERFORMANCECLASS = 'CExpPhasPerformance'
    filecaught = False

    def __init__(self, *args, **kwargs):
        CPluginScript.__init__(self, *args, **kwargs)
        self.MTZ = None
        self.shelx_fname = "shelxrun"

    # Override functions that are needed for run
    def process(self):
        if not self.filecaught:
            jobDirectory = CCP4Modules.PROJECTSMANAGER().db().jobDirectory(jobId=self.jobId)
            self.watchDirectory(jobDirectory, handler=self.handleDirectoryChanged)
        super().process()

    def handleOutputChanged(self, logfile):
        self.parseLogfile()

    def handleDirectoryChanged(self, dirname):
        if not self.filecaught:
            logf = self.makeFileName('LOG')
            if os.path.exists(logf):
                self.watchFile(logf, handler=self.handleOutputChanged, minDeltaSize=34, unwatchWhileHandling=True)
                self.filecaught = True

    def processInputFiles(self):
        print("Processing the input files(KJS-24/09-shelxeMR)")
        # Only two things need to be done here. The .pdb file needs renaming & the genHKL function needs firing up.
        cols = [['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]]
        if self.container.inputData.FREERFLAG.isSet():
            cols.append('FREERFLAG')
        self.hklin, columns, error = self.makeHklInput(cols, extendOutputColnames=True, useInputColnames=True)
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        # Rename the .pdb file to .pda
        pdb_fullpath = self.container.inputData.XYZIN.fullPath.__str__()
        print("Using mtz file ", self.hklin)
        print("Copying ", pdb_fullpath)
        print("to ", os.path.join(self.getWorkDirectory(), self.shelx_fname + '.pda'))
        shutil.copy(pdb_fullpath, os.path.join(self.getWorkDirectory(), self.shelx_fname + '.pda'))
        shutil.copy(self.hklin, os.path.join(self.getWorkDirectory(), self.shelx_fname + '.mtz'))
        # Process the mtz file in order to generate a shelxe friendly hkl file.
        self.genHKL()
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        print("Processing the output files(KJS-24/09-shelxeMR)")
        pdbfile_out = os.path.join(self.getWorkDirectory(), "shelxrun.pdb")
        if os.path.exists(pdbfile_out):
            self.container.outputData.XYZOUT = pdbfile_out
        self.container.outputData.XYZOUT.annotation = "Co-ordinate file for model built by Shelxe"
        self.convToMTZ()
        # Parse the shelxe logfile & then split the output mtz into mini-mtz's
        self.parseLogfile()
        outputFiles = ['FPHIOUT']
        outputColumns = ['FWT,PHWT']
        error = self.splitHklout(outputFiles, outputColumns, os.path.join(self.getWorkDirectory(), 'sftools_out.mtz'))
        if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED

    def parseLogfile(self):
        from mrbump.parsers.parse_shelxe import ShelxeLogParser

        # Parse the shelxe logfile
        rootNode = etree.Element("shelxeMR")
        sxlog = ShelxeLogParser(self.makeFileName('LOG'))
        xmlRI = etree.SubElement(rootNode, "RunInfo")
        xmlbcyc = etree.SubElement(xmlRI, "BestCycle")
        etree.SubElement(xmlbcyc, "BCycle").text = str(sxlog.cycle)
        etree.SubElement(xmlbcyc, "BestCC").text = str(sxlog.CC)
        etree.SubElement(xmlbcyc, "ChainLen").text = str(sxlog.avgChainLength)
        etree.SubElement(xmlbcyc, "NumChains").text = str(sxlog.numChains)
        for iCyc, cycDat in enumerate(sxlog.cycleData):
            xmlcyc = etree.SubElement(xmlRI, "Cycle")
            etree.SubElement(xmlcyc, "NCycle").text = str(iCyc + 1)
            etree.SubElement(xmlcyc, "CorrelationCoef").text = str(cycDat[0])
            etree.SubElement(xmlcyc, "AverageChainLen").text = str(cycDat[1])
            etree.SubElement(xmlcyc, "NumChains").text = str(cycDat[3])
        xmlfile = open(self.xmlout, 'wb')
        xmlString= etree.tostring(rootNode, pretty_print=True)
        xmlfile.write(xmlString)

    def makeCommandAndScript(self, container=None):
        print("Constructing command script (KJS-24/09-shelxeMR)")
        pdbfile_name = self.shelx_fname + '.pda'
        self.appendCommandLine(pdbfile_name)
        self.appendCommandLine("-a%s"%(str(self.container.inputData.NTCYCLES)))
        self.appendCommandLine("-m%s"%(str(self.container.inputData.NMCYCLES)))
        if self.container.inputData.SALPHELICE:
            self.appendCommandLine("-q")
        if self.container.inputData.SBETA and self.container.inputData.SANTIBETA:
            self.appendCommandLine("-B3")
        elif self.container.inputData.SBETA and not self.container.inputData.SANTIBETA:
            self.appendCommandLine("-B2")
        elif not self.container.inputData.SBETA and self.container.inputData.SANTIBETA:
            self.appendCommandLine("-B1")
        self.appendCommandLine("-s%s"%(str(self.container.inputData.FSOLVENT)))
        if self.container.inputData.PRUNRES:
            self.appendCommandLine("-o")
        if self.container.inputData.USENFOLD:
            self.appendCommandLine("-n")
        self.appendCommandLine("-t%s"%(str(self.container.inputData.TIMEFAC)))
        self.xmlout = self.makeFileName('PROGRAMXML')
        return CPluginScript.SUCCEEDED

##================================ Local functions used in data run =========================================
    # Generate .hkl file from .mtz reflection file.
    def genHKL(self):
        from mrbump.file_info import MTZ_parse

        # Define the binary file (with path), as well as the mtz2various logfile
        binf = os.path.normpath(os.path.join( CCP4Utils.getCCP4Dir().__str__(), 'bin', 'mtz2various' ))
        # input & output file (same name as model+.hkl). Also logfile name.
        infile = os.path.join(self.getWorkDirectory(),self.shelx_fname + '.mtz')
        outfile = os.path.join(self.getWorkDirectory(),self.shelx_fname + '.hkl')
        arglist = ['HKLIN', infile, 'HKLOUT', outfile]
        logFile = os.path.normpath(os.path.join(self.getWorkDirectory(), 'mtz2various.log'))
        self.mtz = MTZ_parse.MTZ_parse()
        self.mtz.run_mtzdmp(infile)
        # Extra input needed to be read-in by mtz2various
        inputText = "LABIN FP=%s SIGFP=%s FREE=%s\n" % (self.mtz.F, self.mtz.SIGF, self.mtz.FreeR_flag)
        inputText += ('output SHELX\n')
        inputText += ('FSQUARED')
        # Fire the hkl conversion up.
        pid = CCP4Modules.PROCESSMANAGER().startProcess(binf, arglist, inputText=inputText, logFile=logFile)
        status = CCP4Modules.PROCESSMANAGER().getJobData(pid)
        exitCode = CCP4Modules.PROCESSMANAGER().getJobData(pid, 'exitCode')
        if status == 0 and os.path.exists(outfile):
            return CPluginScript.SUCCEEDED
        self.appendErrorReport(exitCode) # was 201 ?!?
        return CPluginScript.FAILED

    def convToMTZ(self):
        binf2m = os.path.normpath(os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'bin', 'f2mtz'))
        bincad = os.path.normpath(os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'bin', 'cad'))
        binsft = os.path.normpath(os.path.join(CCP4Utils.getCCP4Dir().__str__(), 'bin', 'sftools'))
        pdb_fullpath = self.container.inputData.XYZIN.fullPath.__str__()
        infile1 = os.path.join(self.getWorkDirectory(), self.shelx_fname + '.phs')
        outfile1 = os.path.join(self.getWorkDirectory(), 'f2mtz_phs.mtz')
        # This is an example of the format; need to get actual cell dims. from one of the ascii data files proper.
        # TITLE   f2mtz
        # cell    73.530  39.060  23.150  90.00  90.00  90.00
        # symm    P212121
        # labout  H K L F FOM PHI SIGF
        # CTYPOUT H H H F W P Q
        # pname   f2mtz
        # dname   f2mtz
        # END
        infile = open(pdb_fullpath, "r")
        cell = "1.0  1.0  1.0  90.0  90.0  90.0"
        symm = "P1"
        for cline in infile:
            if re.match(r"CRYST1", cline):
                cell = "  ".join(cline.split()[1:7])
                lenn = min(len(cline),66)
                symm = "".join(cline[55:lenn].split())
                break
        infile.close()
        arglist1 = ['hklin ', infile1, 'hklout', outfile1]
        inputText1 = 'TITLE   f2mtz\n'
        inputText1 += ('cell    %s\n') % cell
        inputText1 += ('symm    %s\n') % symm
        inputText1 += ('labout  H K L F FOM PHI SIGF\n')
        inputText1 += ('CTYPOUT H H H F W P Q\n')
        inputText1 += ('pname   f2mtz\n')
        inputText1 += ('dname   f2mtz\n')
        inputText1 += ('END')
        # LABIN FILE 1 E1=F E2=SIGF E3=PHI E4=FOM
        # LABIN FILE 2 E1=FreeR_flag
        # LABOUT FILE 1 E1=F E2=SIGF E3=PHI_SHELXE E4=FOM_SHELXE
        # LABOUT FILE 2 E1=FreeR_flag
        # RESOLUTION OVERALL 200.0 2.30
        # END
        outfile2 = os.path.join(self.getWorkDirectory(),'cad_out.mtz')
        arglist2 = ['hklin1', outfile1, 'hklin2', self.hklin, 'hklout', outfile2]
        inputText2 = 'LABIN FILE 1 E1=F E2=SIGF E3=PHI E4=FOM\n'
        inputText2 += ('LABIN FILE 2 E1=%s\n') % self.mtz.FreeR_flag
        inputText2 += ('LABOUT FILE 1 E1=F E2=SIGF E3=PHI_SHELXE E4=FOM_SHELXE\n')
        inputText2 += ('LABOUT FILE 2 E1=%s\n') % self.mtz.FreeR_flag
        inputText2 += ('RESOLUTION OVERALL 200.0 %.2lf\n') % float(self.mtz.resolution)
        inputText2 += ('END')
        # mode batch
        # read "cad_out.mtz" mtz
        # CALC F COL FWT = COL F COL FOM_SHELXE *
        # CALC P COL PHWT = COL PHI_SHELXE 0 +
        # write "..\toxd_shelxe4mr_soln1.mtz" mtz
        # EXIT
        # YES
        arglist3 = []
        inputText3 = 'mode batch\n'
        txtFrg1 = os.path.normpath(os.path.join(self.getWorkDirectory(), 'cad_out.mtz'))
        txtFrg2 = os.path.normpath(os.path.join(self.getWorkDirectory(), 'sftools_out.mtz'))
        inputText3 += ('read "%s" mtz\n') % txtFrg1
        inputText3 += ('CALC F COL FWT = COL F COL FOM_SHELXE *\n')
        inputText3 += ('CALC P COL PHWT = COL PHI_SHELXE 0 +\n')
        inputText3 += ('write "%s" mtz\n') % txtFrg2
        inputText3 += ('EXIT\n')
        inputText3 += ('YES')
        logFile1 =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'f2mtzCh.log'))
        logFile2 =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'cad.log'))
        logFile3 =  os.path.normpath(os.path.join(self.getWorkDirectory(), 'sftools.log'))
        pid1 = CCP4Modules.PROCESSMANAGER().startProcess(binf2m, arglist1, inputText=inputText1, logFile=logFile1)
        stat1 = CCP4Modules.PROCESSMANAGER().getJobData(pid1)
        exCd1 = CCP4Modules.PROCESSMANAGER().getJobData(pid1, 'exitCode')
        pid2 = CCP4Modules.PROCESSMANAGER().startProcess( bincad, arglist2, inputText=inputText2, logFile=logFile2)
        stat2 = CCP4Modules.PROCESSMANAGER().getJobData(pid2)
        exCd2 = CCP4Modules.PROCESSMANAGER().getJobData(pid2, 'exitCode')
        pid3 = CCP4Modules.PROCESSMANAGER().startProcess(binsft, arglist3, inputText=inputText3, logFile=logFile3)
        stat3 = CCP4Modules.PROCESSMANAGER().getJobData(pid3)
        exCd3 = CCP4Modules.PROCESSMANAGER().getJobData(pid3, 'exitCode')
        return CPluginScript.SUCCEEDED
