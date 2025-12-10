from __future__ import print_function
from future.utils import raise_
"""
     crank2.py: CCP4 GUI Project
     Copyright (C) 2010 University of York, Leiden University

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

from core.CCP4PluginScript import CPluginScript
from core import CCP4ErrorHandling, CCP4Utils, CCP4XtalData
from core import CCP4Modules
from pipelines.crank2.script import crank2_basepipe

import sys,os,shutil

crank2_path=os.path.join(CCP4Utils.getCCP4I2Dir(),'pipelines','crank2','crank2')
sys.path.append( crank2_path )

class crank2(CPluginScript):

  #TASKMODULE       = 'expt_phasing'
  TASKTITLE        = 'Crank2'
  SHORTTASKTITLE   = 'CRANK2'
  TASKNAME         = 'crank2'
  TASKCOMMAND      = 'crank2.py'
  TASKVERSION      = 0.02
  PERFORMANCECLASS = 'CExpPhasPerformance'
  ERROR_CODES = { 0: {'description': ' '}, }
  MAINTAINER = 'skubakp@gmail.com'


  def has_cont_attr(self,cont,strng):
    try:
      return hasattr(cont,strng)
    except:
      return False

  def CheckUse(self, cont, param):
    # use this function for all parameters that can be set by default or by user
    param_obj = getattr(cont, param)
    is_set = param_obj.isSet()
    has_user = self.has_cont_attr(cont, 'USER_'+param)
    user_val = getattr(cont, 'USER_'+param) if has_user else None
    print(f"[DEBUG CheckUse] param={param}, isSet={is_set}, has_USER={has_user}, USER_val={user_val}, value={param_obj}")
    if is_set and ( not has_user or user_val ):
      return True
    else:
      return False

  def GetProgramKeysOrArgs(self, key_str, args=False, onechar=False):
    keys = key_str.split(',')
    keys2 = []
    for k in keys:
      key,sep,val = k.strip().partition(' ')
      key,val=key.strip(),val.strip()
      if args:
        key,val=key.strip('-'),val.strip('-')
      if len(key)>1 and val=='' and onechar:
        val=key[1:]
        key=key[0]
      if val=='':
        val = 'True'
      if (not val.startswith("\"") and not key.startswith("\"")) or not val.endswith("\""):
        val = val.split()
      else:
        val = [val,]
      if key and args:
        for v in val:
          keys2.append('{};{}'.format(key,v))
      if key and not args:
        for v in val:
          keys2.append('{}:{}'.format(key,v))
    return keys2

  def process(self, container=None):

    defaults = False
    if container:
      self.container = container
      defaults = True
    # fixbrokenpluginname:
    try:
      if self.parent() and self.parent().pluginName.startswith('crank2_'):
        self.container.saveDataToXml(self.parent().comFilePath)
    except RuntimeError:
      pass
      #sys.exc_clear()
    basepipe = crank2_basepipe.crank2_basepipe()
    basepipe.SetBaseSteps(self.container)
    # for convenience
    inp  = self.container.inputData
    ctrl = self.container.controlParameters

    # Debug: show input file info
    print(f"[DEBUG crank2 process()] Starting to build crank_lines")
    print(f"[DEBUG crank2] inp.NON_MTZ = {inp.NON_MTZ}")
    print(f"[DEBUG crank2] F_SIGFanom type: {type(inp.F_SIGFanom)}")
    print(f"[DEBUG crank2] F_SIGFanom.baseName: {inp.F_SIGFanom.baseName}")
    print(f"[DEBUG crank2] F_SIGFanom.fullPath: {inp.F_SIGFanom.fullPath}")
    print(f"[DEBUG crank2] F_SIGFanom.fullPath.isSet(): {inp.F_SIGFanom.fullPath.isSet()}")
    print(f"[DEBUG crank2] F_SIGFanom.contentFlag: {inp.F_SIGFanom.contentFlag}")
    print(f"[DEBUG crank2] F_SIGFanom.contentFlag.isSet(): {inp.F_SIGFanom.contentFlag.isSet()}")

    crank_lines = []
    fpfpp = {}
    cell, spgr = '', ''
    if inp.NON_MTZ and (inp.USER_CELL_A or inp.USER_CELL_B or inp.USER_CELL_C or inp.USER_CELL_D or inp.USER_CELL_E or inp.USER_CELL_F):
      cell = 'cell={},{},{},{},{},{}'.format(inp.CELL_A,inp.CELL_B,inp.CELL_C,inp.CELL_D,inp.CELL_E,inp.CELL_F)
    if inp.USER_SPACEGROUP and inp.NON_MTZ:
      spgr = '"spgr={}"'.format(inp.SPACEGROUP)
    for i in range(4):
      si = str(i+1)  if i  else ''
      sif = si+'_nonmtz'  if inp.NON_MTZ  else si
      anom = getattr(inp,'F_SIGFanom'+sif)
      dn = getattr(inp,'DNAME'+si)
      dnstr = "dname="+str(dn)  if dn  else ""
      print(f"[DEBUG crank2] i={i}, F_SIGFanom{sif}: fullPath.isSet()={anom.fullPath.isSet()}, fullPath={anom.fullPath}, contentFlag.isSet()={anom.contentFlag.isSet()}, contentFlag={anom.contentFlag}")
      if anom.fullPath.isSet() and (anom.contentFlag.isSet() or inp.NON_MTZ) and (not i or getattr(inp,'MAD'+si)):
        saved_fpm = getattr(inp,'SAVED_FPMFILE'+si)
        if saved_fpm and defaults and os.path.isfile(str(saved_fpm)):
          crank_lines.append("fsigf plus {1} f={2} sigf={3} \"file={0}\" {4} {5}".format( saved_fpm, dnstr,
                             getattr(inp,'SAVED_FPLUS'+si),getattr(inp,'SAVED_SIGFPLUS'+si),cell,spgr ))
          crank_lines.append("fsigf minus {0} f={1} sigf={2} {3} {4}".format( dnstr,
                             getattr(inp,'SAVED_FMIN'+si),getattr(inp,'SAVED_SIGFMIN'+si),cell,spgr ))
          crank_lines.append("fsigf average {0} f={2} sigf={3} \"file={1}\" {4} {5}".format( dnstr, getattr(inp,'SAVED_FAVFILE'+si),
                             getattr(inp,'SAVED_FAVER'+si),getattr(inp,'SAVED_SIGFAVER'+si),cell,spgr ))
        elif anom.contentFlag.isSet() and anom.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_FPAIR:
          crank_lines.append("fsigf plus {} f=Fplus sigf=SIGFplus \"file={}\" {} {}".format(dnstr,anom.fullPath,cell,spgr))
          crank_lines.append("fsigf minus {} f=Fminus sigf=SIGFminus {} {}".format(dnstr,cell,spgr))
        elif inp.NON_MTZ or anom.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR:
          crank_lines.append("fsigf plus {1} i=Iplus sigi=SIGIplus \"file={0}\" {2} {3}".format(anom.fullPath,dnstr,cell,spgr))
          crank_lines.append("fsigf minus {} i=Iminus sigi=SIGIminus {} {}".format(dnstr,cell,spgr))
        if self.CheckUse(inp,'WAVELENGTH'+si) and crank_lines:
          crank_lines[-1] += " wavel={}".format(getattr(inp,'WAVELENGTH'+si))
        if self.CheckUse(inp,'FPRIME'+si) or self.CheckUse(inp,'FDPRIME'+si):
          fpfpp[dn] = [getattr(inp,'FPRIME'+si), getattr(inp,'FDPRIME'+si)]

    native = inp.F_SIGFnative_nonmtz  if inp.NON_MTZ  else inp.F_SIGFnative
    if inp.NATIVE and native.fullPath.isSet() and native.contentFlag.isSet():
      saved_aver = getattr(inp,'SAVED_FAVFILE_NATIVE')
      if saved_aver and os.path.isfile(str(saved_aver)) and (defaults or not '.mtz' in native.fullPath):
        crank_lines.append("fsigf average xname=native dname=nat f={} sigf={} \"file={}\" {}".format( \
          inp.SAVED_FAVER_NATIVE,inp.SAVED_SIGFAVER_NATIVE,inp.SAVED_FAVFILE_NATIVE,cell,spgr))
      elif native.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN:
        crank_lines.append("fsigf average xname=native dname=nat f=F sigf=SIGF \"file={}\" {} {}".format(native.fullPath,cell,spgr))
      if native.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN:
        crank_lines.append("fsigf average xname=native dname=nat i=I sigi=SIGI \"file={}\" {} {}".format(native.fullPath,cell,spgr))
      elif native.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR:
        crank_lines.append("fsigf plus xname=native dname=nat i=Iplus sigi=SIGIplus \"file={0}\" {1} {2}".format(native.fullPath,cell,spgr))
        crank_lines.append("fsigf minus xname=native dname=nat i=Iminus sigi=SIGIminus {} {}".format(cell,spgr))
      elif native.contentFlag == CCP4XtalData.CObsDataFile.CONTENT_FLAG_FPAIR:
        crank_lines.append("fsigf plus xname=native dname=nat f=Fplus sigf=SIGFplus \"file={0}\" {1} {2}".format(native.fullPath,cell,spgr))
        crank_lines.append("fsigf minus xname=native dname=nat f=Fminus sigf=SIGFminus {} {}".format(cell,spgr))
      #else:
      #  print 'WARNING: Native inputted as anomalous pairs, mean is needed!  Data ignored.'


    if inp.INPUT_PARTIAL:
      model = "model unknown"
      if inp.XYZIN.isSet():
        model += " \"file={0}\"".format(inp.XYZIN.fullPath)
        if self.CheckUse(ctrl,'COMB_PHDMMB_NCS_DET_MR') and ctrl.COMB_PHDMMB_NCS_DET_MR:
          model += " custom=ncs"
      if inp.ATOM_TYPE.isSet():
        model += " atomtype={0}".format(inp.ATOM_TYPE)
        if self.CheckUse(inp,'DNAME'):
          model += " d_name={0}".format(inp.DNAME)
        if self.CheckUse(inp,'FPRIME') and self.CheckUse(inp,'FDPRIME'):
          model += " fp={0} fpp={1}".format(inp.FPRIME, inp.FDPRIME)
      crank_lines.append("{0}".format(model))
    # why was this code here?  it inputs the partial model if partial model was unclicked...
    #elif inp.XYZIN.isSet():
    #  crank_lines.append("model partial \"file={0}\"".format(inp.XYZIN.fullPath))

    if inp.ATOM_TYPE.isSet():
      model                                 = "model substr atomtype={0}".format(inp.ATOM_TYPE)
      if inp.XYZIN_SUB.isSet():
        #substr_file = str(inp.XYZIN_SUB)
        #if os.path.basename(str(inp.XYZIN_SUB.fullPath)).startswith('SHELX_fa'):
        #  substr_file = os.path.join(self.workDirectory,'SHELX_fa.res')
        #  shutil.copy(str(inp.XYZIN_SUB),substr_file)
        model += " \"file={}\"".format(inp.XYZIN_SUB)
      if inp.XYZIN_SUB_RES.isSet():
        model += " \"file={}\"".format(inp.XYZIN_SUB_RES)
      for dn,f in fpfpp.items():
        if dn:
          model += " d_name={}".format(dn)
        model += " fp={} fpp={}".format(f[0], f[1])
      if not fpfpp and self.CheckUse(inp,'DNAME'):
          model += " d_name={0}".format(inp.DNAME)
      #if self.CheckUse(inp,'NUMBER_SUBSTRUCTURE') and inp.MONOMERS_ASYM.isSet():
      #  model                              += " exp_num_atoms={0}".format(inp.NUMBER_SUBSTRUCTURE*inp.MONOMERS_ASYM)
      if self.CheckUse(inp,'NUMBER_SUBSTRUCTURE'):
        model                              += " exp_num_atoms={0}".format(inp.NUMBER_SUBSTRUCTURE)
      if self.CheckUse(inp,'RESIDUES_MON'):
        model                             += " residues_mon={0}".format(inp.RESIDUES_MON)
      if self.CheckUse(inp,'MONOMERS_ASYM'):
        model                             += " monomers_asym={0}".format(inp.MONOMERS_ASYM)
      if self.CheckUse(inp,'SOLVENT_CONTENT'):
        model                             += " solvent_content={0}".format(inp.SOLVENT_CONTENT)
      crank_lines.append("{0}".format(model))

    # *** NSP - what about DNA???
    if inp.SEQIN.isSet() and inp.INPUT_SEQUENCE:
      
      seqFileName = os.path.join(self.workDirectory,'SEQIN.fasta')
      inp.SEQIN.writeFasta(seqFileName)
      sequence = "sequence \"file={0}\"".format(seqFileName)
      crank_lines.append("{0}".format(sequence))

    if inp.FPHIN_HL.isSet() and inp.INPUT_PHASES:
      if inp.FPHIN_HL.contentFlag==CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL:
        phases = "mapcoef \"file={0}\" hla=HLA hlb=HLB hlc=HLC hld=HLD".format(inp.FPHIN_HL.fullPath)
      else:
        phases = "mapcoef \"file={0}\" ph=PHI fom=FOM".format(inp.FPHIN_HL.fullPath)
      crank_lines.append(phases)

    if self.CheckUse(inp,'EXPTYPE'):
      crank_lines.append("target::{0}".format(inp.EXPTYPE))

    #if inp.REPLACE_MET_MSE.isSet() and inp.REPLACE_MET_MSE:
    #  crank_lines.append("replace_met_mse::1")

    if str(inp.FREE)=='new':
      no_out_next = 'no_output_to_next_step::True'  if not inp.INPUT_PARTIAL  else ''
      crank_lines.append("createfree {} fraction::{}".format(no_out_next, inp.FREE_RATIO*0.01))

    self.i2_shelxdir = None
    if hasattr(CCP4Modules.PREFERENCES(),'SHELXDIR') and CCP4Modules.PREFERENCES().SHELXDIR:
      self.i2_shelxdir = str(CCP4Modules.PREFERENCES().SHELXDIR)

    if basepipe.ToggleDetection() or (basepipe.ToggleShelxCDE() and (inp.XYZIN_SUB.isSet() or inp.XYZIN_SUB_RES.isSet())):
      faest = "faest"
      # this should be removed once we have a "dummy" 'program' keyowrd in crank2!
      if not ctrl.FAEST_PROGRAM.isSet() and self.i2_shelxdir and os.path.isfile(os.path.join(self.i2_shelxdir,'shelxc')):
        ctrl.FAEST_PROGRAM.set('shelxc')
      if ctrl.FAEST_PROGRAM.isSet():
        faest += " {}".format(ctrl.FAEST_PROGRAM)
        if ctrl.FAEST_PROGRAM.startswith('shelx') and self.i2_shelxdir:
          faest += ' "binary::{}"'.format(os.path.join(self.i2_shelxdir,str(ctrl.FAEST_PROGRAM)))
      crank_lines.append(faest)

      if basepipe.ToggleDetection():
        substrdet    = "substrdet"
        if self.CheckUse(ctrl,'SUBSTRDET_HIGH_RES_CUTOFF'):
          substrdet += " high_res_cutoff::{0}".format(ctrl.SUBSTRDET_HIGH_RES_CUTOFF) 
        if self.CheckUse(ctrl,'SUBSTRDET_HIGH_RES_CUTOFF_CCHALF'):
          substrdet += " high_res_cutoff_cchalf::{0}".format(bool(ctrl.SUBSTRDET_HIGH_RES_CUTOFF_CCHALF))
        if self.CheckUse(ctrl,'SUBSTRDET_HIGH_RES_CUTOFF_RADIUS') and str(ctrl.SUBSTRDET_PROGRAM)=='prasa':
          substrdet += " high_res_cutoff_radius::{0}".format(ctrl.SUBSTRDET_HIGH_RES_CUTOFF_RADIUS) 
        if self.CheckUse(ctrl,'SUBSTRDET_HIGH_RES_CUTOFF_STEP'):
          substrdet += " high_res_cutoff_step::{0}".format(ctrl.SUBSTRDET_HIGH_RES_CUTOFF_STEP) 
        if self.CheckUse(ctrl,'SUBSTRDET_THRESHOLD_STOP'):
          substrdet += ' threshold_stop::{0}'.format(ctrl.SUBSTRDET_THRESHOLD_STOP) 
        if self.CheckUse(ctrl,'SUBSTRDET_THRESHOLD_WEAK'):
          substrdet += ' threshold_weak::{0}'.format(ctrl.SUBSTRDET_THRESHOLD_WEAK) 
        if self.CheckUse(ctrl,'SUBSTRDET_NUM_TRIALS'):
          substrdet += " num_trials::{0}".format(ctrl.SUBSTRDET_NUM_TRIALS)
        if self.CheckUse(ctrl,'SUBSTRDET_MIN_DIST_ATOMS'):
          minus='-' if not str(ctrl.SUBSTRDET_MIN_DIST_ATOMS).startswith('-') else ''
          substrdet += " min_dist_atoms::{}{}".format(minus,ctrl.SUBSTRDET_MIN_DIST_ATOMS)
        if self.CheckUse(ctrl,'SUBSTRDET_MIN_DIST_SYMM_ATOMS'):
          substrdet += " min_dist_symm_atoms::{}".format(-0.1 if ctrl.SUBSTRDET_MIN_DIST_SYMM_ATOMS else 3)
        if self.CheckUse(ctrl,'SUBSTRDET_OPTIMIZE_SOL'):
          substrdet += " optimize_sol::{}".format(2 if ctrl.SUBSTRDET_OPTIMIZE_SOL else 0)
        if self.CheckUse(ctrl,'SUBSTRDET_NUM_THREADS'):
          substrdet += " num_threads::{0}".format(ctrl.SUBSTRDET_NUM_THREADS)
        if ctrl.SUBSTRDET_NUM_ATOMS: #and inp.NUMBER_SUBSTRUCTURE:
          substrdet += " num_atoms::{}".format(ctrl.SUBSTRDET_NUM_ATOMS) #inp.NUMBER_SUBSTRUCTURE)
        # this should be removed once we have a "dummy" 'program' keyowrd in crank2!
        if not ctrl.SUBSTRDET_PROGRAM.isSet() and self.i2_shelxdir and os.path.isfile(os.path.join(self.i2_shelxdir,'shelxd')):
          ctrl.SUBSTRDET_PROGRAM.set('shelxd')
        if ctrl.SUBSTRDET_PROGRAM.isSet():
          if ctrl.SUBSTRDET_PROGRAM.startswith('shelx') and self.CheckUse(inp,'SUBSTRDET_NUM_DSUL'):
            substrdet += " num_dsul::{}".format(inp.SUBSTRDET_NUM_DSUL)
          substrdet   += " {0}".format(ctrl.SUBSTRDET_PROGRAM)
          if ctrl.SUBSTRDET_PROGRAM.startswith('prasa') and ctrl.PRASA_MINPEAKS:
            substrdet += " minpeaks:{}".format(ctrl.PRASA_MINPEAKS)
          if ctrl.SUBSTRDET_PROGRAM.startswith('shelx') and self.i2_shelxdir:
            substrdet += ' "binary::{}"'.format(os.path.join(self.i2_shelxdir,str(ctrl.SUBSTRDET_PROGRAM)))
        if ctrl.KEYWORDS_SUBSTRDET.isSet():
          substrdet += ' '+' '.join( self.GetProgramKeysOrArgs(ctrl.KEYWORDS_SUBSTRDET) )
          #substrdet                          += '  - \n \"keywords={0}\"'.format(ctrl.KEYWORDS_SUBSTRDET)
        crank_lines.append(substrdet)

    if basepipe.TogglePeakSearch():
      refatompick                           = "refatompick"
      if self.CheckUse(ctrl,'REFATOMPICK_NUM_ITER'):
        refatompick                        += " num_iter::{0}".format(ctrl.REFATOMPICK_NUM_ITER)
      if self.CheckUse(ctrl,'REFATOMPICK_REFCYC'):
        refatompick                        += " refcyc::{0}".format(ctrl.REFATOMPICK_REFCYC)
      if self.CheckUse(ctrl,'REFATOMPICK_RMS_THRESHOLD'):
        refatompick                        += " rms_threshold::{0}".format(ctrl.REFATOMPICK_RMS_THRESHOLD)
      if self.CheckUse(ctrl,'REFATOMPICK_OCC_CUT'):
        refatompick                        += " occ_cut::{0}".format(ctrl.REFATOMPICK_OCC_CUT)
      #if self.CheckUse(ctrl,'REFATOMPICK_REF_PROGRAM'):
      #  refatompick                        += " ref {0}".format(ctrl.REFATOMPICK_REF_PROGRAM)
      #else:
      #  refatompick                        += " ref"   # Do we need this?? ***NSP
      crank_lines.append(refatompick)
      if inp.INPUT_PARTIAL and inp.PARTIAL_AS_SUBSTR:
        crank_lines.append("sepsubstrprot")

    if basepipe.ToggleShelxCDE():
      phdmmb    = "phdmmb"
      if self.CheckUse(ctrl,'PHDMMB_DMCYC'):
        phdmmb += " dmcyc::{0}".format(ctrl.PHDMMB_DMCYC)
      if self.CheckUse(ctrl,'PHDMMB_BIGCYC'):
        phdmmb += " bigcyc::{0}".format(ctrl.PHDMMB_BIGCYC)
      if self.CheckUse(ctrl,'PHDMMB_THRESHOLD_STOP'):
        phdmmb += " threshold_stop::{0}".format(ctrl.PHDMMB_THRESHOLD_STOP)
      if self.CheckUse(ctrl,'PHDMMB_THRESHOLD_HAND_STOP'):
        phdmmb += " threshold_hand_stop::{0}".format(ctrl.PHDMMB_THRESHOLD_HAND_STOP)
      if self.CheckUse(ctrl,'PHDMMB_THOROUGH_BUILD'):
        phdmmb += " thorough_build::{0}".format(bool(ctrl.PHDMMB_THOROUGH_BUILD))
      if self.CheckUse(inp,'SUBSTR_ATOMS_NATIVE'):
        phdmmb += " substr_in_native::{}".format(str(inp.SUBSTR_ATOMS_NATIVE))
      phdmmb                               += " shelxe"
      if self.i2_shelxdir:
        phdmmb += ' "binary::{}"'.format(os.path.join(self.i2_shelxdir,'shelxe'))
      if ctrl.ARGUMENTS_SHELXE.isSet():
        phdmmb += ' '+' '.join( self.GetProgramKeysOrArgs(ctrl.ARGUMENTS_SHELXE,args=True,onechar=True) )
      crank_lines.append(phdmmb)

    if basepipe.TogglePhasing():
      phas                                  = "phas"
      if self.CheckUse(ctrl,'PHAS_CYCLES'):
        phas                               += " cycles::{0}".format(ctrl.PHAS_CYCLES) 
      if ctrl.PHAS_PROGRAM.isSet():
        phas                               += " {0}".format(ctrl.PHAS_PROGRAM)
      if ctrl.KEYWORDS_PHAS.isSet():
        phas                               += '  - \n \" keywords={0}\"'.format(ctrl.KEYWORDS_PHAS)
      crank_lines.append(phas)

    if basepipe.ToggleHandDetermination() and ctrl.DO_HANDDET:
      handdet                               = "handdet"
      if ctrl.HANDDET_THRESHOLD_DISCRIM.isSet():
        handdet                            += " threshold_discrim::{0}".format(ctrl.HANDDET_THRESHOLD_DISCRIM)
      handdet                              += " dmfull"
      #if ctrl.HANDDET_DMCYC.isSet():
      #  handdet                            += " dmcyc::{0}".format(ctrl.HANDDET_DMCYC)
      if ctrl.HANDDET_DMFULL_DM_PROGRAM.isSet():
        handdet                            += " dm {0}".format(ctrl.HANDDET_DMFULL_DM_PROGRAM)
      else:
        handdet                            += " dm"   # Do we need this?? ***NSP
      if ctrl.KEYWORDS_HANDDET_DM.isSet():
        handdet                            += '  - \n \" keywords={0}\"\n'.format(ctrl.KEYWORDS_HANDDET_DM)
      if ctrl.HANDDET_DMFULL_PHCOMB_PROGRAM.isSet():
        handdet                            += " phcomb {0}".format(ctrl.HANDDET_DMFULL_PHCOMB_PROGRAM)
      else:
        handdet                            += " phcomb" # Do we need this?? ***NSP
      if ctrl.KEYWORDS_HANDDET_PHCOMB.isSet():
        handdet                            += '  - \n \" keywords={0}\"'.format(ctrl.KEYWORDS_HANDDET_PHCOMB)
      crank_lines.append(handdet)

    if basepipe.ToggleDensityModification():
      dmfull                                = "dmfull"
      if self.CheckUse(ctrl,'DMFULL_DMCYC'):
        dmfull                             += " dmcyc::{0}".format(ctrl.DMFULL_DMCYC)
      if self.CheckUse(ctrl,'DMFULL_THRESHOLD_STOP'):
        dmfull                             += " threshold_stop::{0}".format(ctrl.DMFULL_THRESHOLD_STOP)
      if ctrl.DMFULL_DM_PROGRAM.isSet():
        dmfull                             += " dm {0}".format(ctrl.DMFULL_DM_PROGRAM)
      else:
        dmfull                             += " dm"   # Do we need this?? ***NSP  # We do. ***PS
      if ctrl.KEYWORDS_DMFULL_DM.isSet():
        dmfull += ' '+' '.join( self.GetProgramKeysOrArgs(ctrl.KEYWORDS_DMFULL_DM,args=str(ctrl.DMFULL_DM_PROGRAM)!='solomon') )
      #  dmfull                             += '  - \n \" keywords={0}\"\n'.format(ctrl.KEYWORDS_DMFULL_DM)
      if ctrl.DMFULL_PHCOMB_PROGRAM.isSet():
        dmfull                             += " phcomb {0}".format(ctrl.DMFULL_PHCOMB_PROGRAM)
      else:
        dmfull                             += " phcomb" # Do we need this?? ***NSP
      if ctrl.KEYWORDS_DMFULL_PHCOMB.isSet():
        dmfull                             += '  - \n \" keywords={0}\"'.format(ctrl.KEYWORDS_DMFULL_PHCOMB)
      crank_lines.append(dmfull)

    if basepipe.ToggleModelBuilding():
      if basepipe.ToggleUseComb():
        model_build                       = "comb_phdmmb target::SAD"
        if self.CheckUse(ctrl,'COMB_PHDMMB_MINBIGCYC'):
          model_build                    += " minbigcyc::{0}".format(ctrl.COMB_PHDMMB_MINBIGCYC)
        if self.CheckUse(ctrl,'COMB_PHDMMB_MAXBIGCYC'):
          model_build                    += " maxbigcyc::{0}".format(ctrl.COMB_PHDMMB_MAXBIGCYC)
        if self.CheckUse(ctrl,'COMB_PHDMMB_NCS_DET') and basepipe.ToggleNCS():
          model_build                    += " ncs_det::{0}".format(ctrl.COMB_PHDMMB_NCS_DET._value)
        if self.CheckUse(ctrl,'COMB_PHDMMB_NCS_DET_MR') and basepipe.ToggleNCS():
          model_build                    += " ncs_det_mr::{0}".format(ctrl.COMB_PHDMMB_NCS_DET_MR._value)
        if self.CheckUse(ctrl,'COMB_PHDMMB_NUM_PARALLEL'):
          model_build                    += " num_parallel::{0}".format(ctrl.COMB_PHDMMB_NUM_PARALLEL)
        if self.CheckUse(ctrl,'COMB_PHDMMB_SKIP_INITIAL_BUILD'):
          model_build                    += " skip_initial_build::{0}".format(ctrl.COMB_PHDMMB_SKIP_INITIAL_BUILD._value)
        if self.CheckUse(ctrl,'COMB_PHDMMB_REBUILD_ONLY'):
          model_build                    += " rebuild_only::{0}".format(ctrl.COMB_PHDMMB_REBUILD_ONLY._value)
        if self.CheckUse(ctrl,'COMB_PHDMMB_START_SHELXE'):
          model_build                    += " start_shelxe::{0}".format(bool(ctrl.COMB_PHDMMB_START_SHELXE))
        if ctrl.COMB_PHDMMB_EXCLUDE_FREE!='never':
          if ctrl.COMB_PHDMMB_EXCLUDE_FREE=='last':
            model_build += " always_exclude_free::False"
          model_build += self.SetFree()
        if self.CheckUse(ctrl,'COMB_PHDMMB_DMFULL_DM_PROGRAM'):
          model_build                    += " dmfull dm {0}".format(ctrl.COMB_PHDMMB_DMFULL_DM_PROGRAM)
        else:
          model_build                    += " dmfull dm"   # Do we need this?? ***NSP
        if ctrl.KEYWORDS_COMB_DM.isSet():
          model_build += ' '+' '.join( self.GetProgramKeysOrArgs(ctrl.KEYWORDS_COMB_DM,args=str(ctrl.COMB_PHDMMB_DMFULL_DM_PROGRAM)!='solomon') )
          #model_build                  += '  - \n \" keywords={0}\"\n'.format(ctrl.KEYWORDS_COMB_DM)
        if ctrl.KEYWORDS_MB.isSet():
          model_build += ' mb buccaneer '+' '.join( self.GetProgramKeysOrArgs(ctrl.KEYWORDS_MB,args=True) )
        if ctrl.KEYWORDS_COMB_SHELXE.isSet() and self.CheckUse(ctrl,'COMB_PHDMMB_START_SHELXE'):
          model_build += ' shelxe \" keywords={0}\"\n'.format(self.GetProgramKeysOrArgs(ctrl.KEYWORDS_COMB_SHELXE,args=True,onechar=True))
      else:
        model_build                       = "mbref"
        if ctrl.MBREF_EXCLUDE_FREE:
          model_build += self.SetFree()
        if self.CheckUse(ctrl,'MBREF_BIGCYC'):
          model_build                    += " bigcyc::{0}".format(ctrl.MBREF_BIGCYC)
        if ctrl.MB_PROGRAM.isSet() and str(ctrl.MB_PROGRAM)=='arpwarp':
          model_build += " arpwarp"
        else:
          if ctrl.MBREF_REF_PROGRAM.isSet():
            model_build                    += " ref {0}".format(ctrl.MBREF_REF_PROGRAM)
          else:
            model_build                    += " ref"   # Do we need this?? ***NSP
          if ctrl.KEYWORDS_MB.isSet():
            model_build                    += '  - \n \" keywords={0}\"\n'.format(ctrl.KEYWORDS_MB_REF)
  
          if ctrl.MB_PROGRAM.isSet():
            model_build                    += " mb {0}".format(ctrl.MB_PROGRAM)
          else:
            model_build                    += " mb"   # Do we need this?? ***NSP
          if ctrl.KEYWORDS_MB.isSet():
            model_build                    += '  - \n \" keywords={0}\"\n'.format(ctrl.KEYWORDS_MB)
          if ctrl.KEYWORDS_COMB_SHELXE.isSet() and self.CheckUse(ctrl,'COMB_PHDMMB_START_SHELXE'):
            model_build                    += ' shelxe \" keywords={0}\"\n'.format(self.GetProgramKeysOrArgs(ctrl.KEYWORDS_COMB_SHELXE,args=True,onechar=True))

      crank_lines.append(model_build)

    if basepipe.ToggleRefine():
      ref                               = "ref target::mlhl"
      if self.CheckUse(ctrl,'REF_CYCLES'):
        ref                            += " cycles::{0}".format(ctrl.REF_CYCLES)
        # ***NSP when rfree is allowed from input, we need to allow it to be input here...
      if ctrl.REF_EXCLUDE_FREE:
        ref += self.SetFree()
      if ctrl.REF_PROGRAM.isSet():
        ref                            += " {0}".format(ctrl.REF_PROGRAM)
      if ctrl.KEYWORDS_REF.isSet():
        ref                            += '  - \n \" keywords={0}\"\n'.format(ctrl.KEYWORDS_REF)
      crank_lines.append(ref)

    inpfile = os.path.join(self.workDirectory, 'crank.inp')
    with open(inpfile,'w') as g:
      g.write( '\n'.join(crank_lines) )


    if (inp.F_SIGFanom.fullPath.isSet() or inp.F_SIGFanom_nonmtz.fullPath.isSet()) or (defaults and inp.SEQIN.isSet()):
#        (inp.SEQIN.isSet() or inp.RESIDUES_MON.isSet()):
      rvapi_style = str(ctrl.PRESENT_STYLE)  if (ctrl,'PRESENT_STYLE')  else 'rvapi_tree'
      i2natfile = os.path.join(self.workDirectory, 'i2_native_present')
      self.rvapi_converter=False
      if rvapi_style == 'i2_bigpage':
        with open(i2natfile,'w') as g:
          g.write('using i2 native presentation')
        self.rvapi_converter=True
      elif os.path.isfile(i2natfile):
        os.remove(i2natfile)
      import ccp4i2crank,common,traceback
      try:
        crank2 = ccp4i2crank.CallCrankFromCCP4i2(self, inpfile=inpfile, defaults=defaults, rvapi_style=rvapi_style)
        if not defaults and self.has_cont_attr(ctrl,"CLEANUP") and bool(ctrl.CLEANUP):
          self.CleanUp(crank2)
          self.reportStatus(CPluginScript.SUCCEEDED)
        return crank2
#      except common.CrankError as e:
#        self.appendErrorReport(1000, str(e))
#        try:
#          self.reportStatus(CPluginScript.FAILED)
#        except RuntimeError:
#          raise e
      except common.Unsuccessful as e:
        self.appendErrorReport(1000, str(e))
        try:
          self.reportStatus(CPluginScript.UNSATISFACTORY)
        except RuntimeError:
          raise_(e,None,sys.exc_info()[2])
      except Exception as e:
        err = str(e) +'\n\n'+str(traceback.format_exc())
        self.appendErrorReport(0, err)
        try:
          self.reportStatus(CPluginScript.FAILED)
        except RuntimeError:
          raise_(e,None,sys.exc_info()[2])
    else:
      return 0

  def CleanUp(self,crank):
    # pokus - delete all mtz that are not registered as i2 i/o
    all_cont = [crank.ccp4i2job.container,] + [ proc.ccp4i2job.container for proc in crank.processes if hasattr(proc,'ccp4i2job') ]
    protected = [ os.path.realpath(str(getattr(c.outputData,f))) for c in all_cont for f in c.outputData._dataOrder if str(getattr(c.outputData,f)) ]
    protected += [ os.path.realpath(str(getattr(c.inputData,f))) for c in all_cont for f in c.inputData._dataOrder if str(getattr(c.inputData,f)) ]
    count_removed=0
    for root,dirs,files in os.walk( CCP4Modules.PROJECTSMANAGER().jobDirectory(self.jobId) ):
      for f in files:
        if f.lower().endswith('.mtz') and os.path.realpath(os.path.join(root,f)) not in protected:
          os.remove( os.path.join(root,f) )
          count_removed+=1
    print('Removed',count_removed,'intermediary mtz files.')

  def SetFree(self):
    inp = self.container.inputData
    if str(inp.FREE)=='new':
      return " exclude obj_from=0,typ=freeR"
    elif str(inp.FREE)=='existing' and inp.FREERFLAG.isSet():
      return " exclude typ=freeR free=FREER \"file={}\"".format(inp.FREERFLAG.fullPath)
    return ""

  def setProgramVersion(self):
    with open(os.path.join(crank2_path,'VERSION')) as f:
      return f.read()
