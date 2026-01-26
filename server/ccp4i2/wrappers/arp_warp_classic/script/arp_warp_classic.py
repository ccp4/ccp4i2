
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core import CCP4Utils
from ccp4i2.core import CCP4ErrorHandling, CCP4XtalData
from ccp4i2.core import CCP4Modules
import os, sys, shutil

class arp_warp_classic(CPluginScript):

    TASKTITLE = 'ARP/WARP classic'
    TASKNAME = 'arp_warp_classic'
    MAINTAINER = 'andrey.lebedev@stfc.ac.uk'

    TASKHOME = os.path.join(CCP4Utils.getCCP4I2Dir(), 'wrappers', 'arp_warp_classic')
    TASKCOMMAND = sys.executable
    PERFORMANCECLASS = 'CModelBuildPerformance'

    @staticmethod
    def getAsuParams(asuObj):
        maxResidues = asuResidues = asuCopies = 0
        for seqObj in asuObj.fileContent.seqList:
            selectionMode = asuObj.qualifiers('selectionMode')
            if selectionMode == 0 or (not asuObj.selection.isSet()) or asuObj.selection[name]:
                seqResidues = int(seqObj.numberOfResidues(countMulti=True))
                asuResidues += seqResidues
                if seqResidues > maxResidues:
                    asuCopies = int(seqObj.nCopies)

        return asuResidues, asuCopies

    def processInputFiles(self):
      self.report_dir = 'report'
      self.res_file = os.path.join(self.report_dir, 'awa.res')

      cinp =  self.container.inputData
      cpar =  self.container.controlParameters
      cols = list()
      self.script = list()

      arp_mode = cpar.AWA_ARP_MODE.get()
      free_is_set = cinp.AWA_FREE.isSet()
      seqin_is_set = cinp.AWA_SEQIN.isSet()
      ref_mode = cpar.AWA_REF_MODE.get()

      cols.append(['AWA_FOBS', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN])
      self.script.append(('fp', 'AWA_FOBS_F'))
      self.script.append(('sigfp', 'AWA_FOBS_SIGF'))

      if free_is_set:
        cols.append('AWA_FREE')
        self.script.append(('freelabin', 'AWA_FREE_FREER'))

      if arp_mode == 'WARPNTRACEPHASES':
        cols.append(['AWA_PHINI', CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM])
        self.script.append(('phibest', 'AWA_PHINI_PHI'))
        self.script.append(('fom', 'AWA_PHINI_FOM'))

      else:
        self.script.append(('modelin', str(cinp.AWA_MODELIN)))

      if seqin_is_set:
        tmpSeqFile = os.path.join(self.workDirectory,'ARP_SEQIN.txt')
        cinp.AWA_SEQIN.writeArpPir(tmpSeqFile,writeMulti=True)
        self.script.append(('seqin',tmpSeqFile ))
        nres, nmol = self.getAsuParams(cinp.AWA_SEQIN)
        # nres_x = int(cinp.AWA_SEQIN.numberOfResidues(countMulti=True))
        # assert nres_x == nres
        self.script.append(('residues', nres))
        self.script.append(('cgr', nmol))

      if ref_mode == 'AWA_SAD':
        cols.append(['AWA_FOBS', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FPAIR])
        self.script.append(('phaselabin', 'F+=AWA_FOBS_F(+) SIGF+=AWA_FOBS_SIGF(+) F-=AWA_FOBS_F(-) SIGF-=AWA_FOBS_SIGF(-)'))
        self.script.append(('heavyin', cinp.AWA_HEAVYIN.get()))
        ano_option = cpar.AWA_ANO_OPTION.get()
        if ano_option == 'LAMBDA':
          self.script.append(('sadcard', 'ANOM WAVE ' + str(cpar.AWA_SCAT_LAMBDA_AWA_SAD.get())))

        elif ano_option == 'MEASURED':
          scat_atom = cpar.AWA_SCAT_ATOM.get()
          scat_fp_sad = cpar.AWA_SCAT_FP_AWA_SAD.get()
          scat_fdp_sad = cpar.AWA_SCAT_FDP_AWA_SAD.get()
          self.script.append(('sadcard', 'ANOM FORM %s %s %s' %(scat_atom, scat_fp_sad, scat_fdp_sad)))

      elif ref_mode == 'AWA_HL':
        cols.append(['AWA_PHREF', CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL])
        self.script.append(('phaselabin', 'HLA=AWA_PHREF_HLA HLB=AWA_PHREF_HLB HLC=AWA_PHREF_HLC HLD=AWA_PHREF_HLD'))

      elif ref_mode == 'AWA_PHASED':
        cols.append(['AWA_PHREF', CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM])
        self.script.append(('phaselabin', 'PHIB=AWA_PHREF_PHI FOM=AWA_PHREF_FOM'))

      self.hklin, columns, error = self.makeHklin0(cols)
      if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
        return CPluginScript.FAILED




      self.script.append(('buildingcycles', cpar.AWA_BIG_CYCLES.get()))
      small_cycles = cpar.AWA_SMALL_CYCLES.get()
      self.script.append(('restrref', small_cycles))
      use_cond = cpar.AWA_USE_COND.get()
      force_cond = cpar.AWA_FORCE_COND.get()
      self.script.append(('restraints', 2 if force_cond and use_cond else 1 if use_cond else 0))
      if cpar.AWA_FAKE_DATA.get():
        self.script.append(('fakedata', '0.33 0.75 1'))

      self.script.append(('ncsrestraints', int(cpar.AWA_NCS_RESTRAINTS.get())))
      self.script.append(('ncsextension', int(cpar.AWA_NCS_EXTENSION.get())))
      self.script.append(('loops', int(cpar.AWA_LOOPS.get())))
      self.script.append(('side', cpar.AWA_SIDE_AFTER.get() if seqin_is_set and cpar.AWA_BUILD_SIDE.get() else -1))
      if ref_mode == 'AWA_SAD':
        self.script.append(('is_semet', int(cpar.AWA_IS_SEMET.get())))

      self.script.append(('albe', int(cpar.AWA_ALBE.get())))
      skip_build = cpar.AWA_SKIP_BUILD.get()
      skip_cycles = cpar.AWA_SKIP_CYCLES.get()
      self.script.append(('cycskip', skip_cycles* small_cycles if skip_build else 0))
      self.script.append(('multit', cpar.AWA_MULTITRACE.get()))





      if arp_mode == 'WARPNTRACEMODEL':
        self.script.append(('freebuild', int(cpar.AWA_FREEBUILD.get())))
        self.script.append(('flatten', int(cpar.AWA_FLATTEN.get())))

      self.script.append(('fsig', cpar.AWA_ADDATOM_SIGMA.get()))
      self.script.append(('rsig', cpar.AWA_REMATOM_SIGMA.get()))
      self.script.append(('upmore', cpar.AWA_UP_UPDATE.get()))
      if ref_mode != 'AWA_SAD':
        self.script.append(('twin', int(cpar.AWA_TWIN.get())))

      self.script.append(('rrcyc', cpar.AWA_NCYCLES.get()))
      if ref_mode in ('AWA_HL', 'AWA_PHASED'):
        self.script.append(('phaseref', 'PHAS SCBL ' + str(cpar.AWA_PHASE_BLUR.get())))

      weight_mode = cpar.AWA_WEIGHT_MODE.get()
      self.script.append(('wmat', 'AUTO' if weight_mode == 'AUTO' else 'MATRIX'))
      if weight_mode != 'AUTO':
        self.script.append(('weightv', cpar.AWA_WMAT.get()))

      self.script.append(('ridgerestraints', int(cpar.AWA_RIDGE_RESTRAINTS.get())))
      scale = cpar.AWA_SCALE.get()
      self.script.append(('scaleopt', scale + ' LSSC ANIS' if cpar.AWA_SCANIS.get() else scale))
      self.script.append(('scalml', 'SCAL MLSC FREE' if free_is_set and cpar.AWA_REFMAC_REF_SET.get() else 'SCAL MLSC'))
      self.script.append(('solvent', int(cpar.AWA_SOLVENT.get())))
      self.script.append(('reportdir', self.report_dir))
      self.script.append(('resfile', self.res_file))

      if cpar.AWA_MOCKYES.get():
        self.script.append(('mockdir', os.path.join(self.TASKHOME, 'mock_data', 'PSP')))
        self.script.append(('pause', cpar.AWA_MOCKPAUSE.get()))
        if cpar.AWA_JSRVIEW.get():
          report_path = os.path.join(self.workDirectory, self.report_dir)
          if not os.path.exists(report_path):
            os.mkdir(report_path)

          jsrview = os.path.join(os.environ['CCP4'], 'libexec', 'jsrview')
          if sys.platform.startswith('win'):
            jsrview += '.exe'

          cmd = [jsrview, os.path.join(report_path, 'index.html')]
          import subprocess
          subprocess.Popen(cmd)

      self.script.append(('workdir', self.workDirectory))
      self.script.append(('jobId', 'PSP'))
      self.script.append(('datafile', self.hklin))
      self.script.append(('xmlout', self.makeFileName('PROGRAMXML')))

      return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
      for key, value in self.script:
        if not value is None:
          tvalue = str(value).strip()
          if tvalue:
            self.appendCommandScript(key + ' ' + tvalue)

      self.appendCommandLine('-m')
      self.appendCommandLine('pyrvapi_ext.parsers.arpwarp')
      return CPluginScript.SUCCEEDED

    def  processOutputFiles(self):
      awa_prefix, mtz_ext, ending = os.path.basename(self.hklin).rpartition('.mtz')
      assert mtz_ext == '.mtz' and not ending
      awa_oroot = os.path.join(self.workDirectory, 'PSP', awa_prefix + '_warpNtrace')
      del awa_prefix, mtz_ext, ending

      annotation = self.jobNumberString() + ' 2mFo-DFc map coefficients from ARP/WARP'
      self.container.outputData.FPHIOUT.annotation = annotation
      self.container.outputData.FPHIOUT.subType = 1

      annotation = self.jobNumberString() + ' mFo-DFc map coefficients from ARP/WARP'
      self.container.outputData.DIFFPHIOUT.annotation = annotation
      self.container.outputData.DIFFPHIOUT.subType = 2

      xyzwrk = awa_oroot + '.pdb'
      if os.path.isfile(xyzwrk):
        #shutil.copy(xyzwrk, xyzout)
        # LP - make sure there are no 'DUM' atoms in the file - save any to 'dummy_atoms.pdb'
        #from ccp4i2.core import CCP4ModelData
        #pdbobj = CCP4ModelData.CPdbDataFile(xyzwrk)
        #pdbobj.removeDummyAtoms(str(self.container.outputData.XYZOUT),os.path.join(self.workDirectory,'XYZOUT_dummy_atoms.pdb'))
        # AL - the above code is broken and replaced with a hack below;
        # the arp-warp output with dummy atoms added to report.
        with open(xyzwrk) as istream:
          pdb_records = istream.read()

        import re
        match_dummy = re.search('^ATOM .+ DUM +DUM ', pdb_records, flags=re.M)
        if match_dummy:
          xyzdum = str(self.container.outputData.XYZDUM.fullPath)
          annotation = self.jobNumberString() + ' ARP/WARP model with dummy atoms'
          self.container.outputData.XYZDUM.annotation = annotation
          self.container.outputData.XYZDUM.subType = 1
          with open(xyzdum, 'w') as ostream:
            ostream.write(pdb_records)

          pdb_records = pdb_records[:match_dummy.start()]

        xyzout = str(self.container.outputData.XYZOUT.fullPath)
        annotation = self.jobNumberString() + ' ARP/WARP model'
        self.container.outputData.XYZOUT.annotation = annotation
        self.container.outputData.XYZOUT.subType = 1
        with open(xyzout, 'w') as ostream:
          ostream.write(pdb_records)

      else:
        return CPluginScript.FAILED

      hklawa = awa_oroot + '.mtz'
      outputFiles = ['FPHIOUT','DIFFPHIOUT']
      outputColumns = ['FWT,PHWT','DELFWT,PHDELWT']
      error = self.splitHklout(outputFiles, outputColumns, infile=hklawa)

      if error.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
        return CPluginScript.FAILED

      try:
        with open(os.path.join(self.workDirectory, self.res_file)) as istream:
            frac_built, r_free, r_cryst = [float(x) for x in istream.read().split()]

        self.container.outputData.PERFORMANCE.RFactor.set(r_cryst)
        if frac_built > 0:
            self.container.outputData.PERFORMANCE.completeness.set(frac_built)

      except:
        pass

      try:
        if not CCP4Modules.PREFERENCES().RETAIN_DIAGNOSTIC_FILES:
          shutil.rmtree(os.path.join(self.workDirectory, 'PSP'))

      except:
        pass

      return CPluginScript.SUCCEEDED

