#!/usr/bin/python
import os,sys,re,shutil,copy
from ..process import process, crvapi
from ..program import program
from .. import common, inout
par=common.parameter

class comb_phdmmb(process):
  name="combined iterative model building with density modification with phased refinement"
  short_name="combined model building"
  supported_procs = ["dmfull", "mb"]
  supported_params={}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True, share=True )
  supported_params['minbigcyc'] = par( desc='Minimal number of model building cycles', typ=int )
  supported_params['maxbigcyc'] = par( desc='Maximal number of model building cycles', typ=int )
  supported_params['low_res_par'] = par( desc='Use low resolution parameters in refinement and model building.', typ=bool, share=True )
  supported_params['refcyc_finish'] = par( desc='Number of refinement cycles in "finishing mode"', typ=int )
  supported_params['ncs_det'] = par( desc='Determine NCS from heavy atoms or partial model (Parrot only)', typ=bool, share=True )
  supported_params['ncs_det_ha'] = par( desc='Determine NCS from heavy atoms (Parrot only)', typ=bool, share=True )
  supported_params['ncs_det_mr'] = par( desc='Determine NCS from partial model (Parrot only)', typ=bool, share=True )
  supported_params['solvent_content'] = par( desc='(Expected) solvent fraction of crystal', typ=float, share=True )
  supported_params['solventmask_radius'] = par( desc='Use the specified solvent mask radius (Parrot only)', typ=float, share=True )
  supported_params['no_bias_est'] = par( desc='Do not estimate bias in density modification', typ=bool, share=True )
  supported_params['skip_initial_build'] = par( desc='Skip model bulding in the first "big" cycle', typ=bool )
  supported_params['start_dm'] = par( desc='Include initial phase improvement without model building', typ=bool )
  supported_params['rebuild_only'] = par( desc='Use correlation mode rebuilding cycles only"', typ=bool )
  supported_params['always_exclude_free'] = par( desc='Exclude inputted free refl. during the whole process (default True, if False then in last mb cycle only)', typ=bool )
  supported_params['optimize_solvent'] = par( desc='Disable solvent content optimization', typ=bool )
  supported_params['R_threshold'] = par( desc='Rcomb-factor threshold for switching to finishing mode', typ=float )
  supported_params['fom_threshold'] = par( desc='FOM threshold for switching to finishing mode', typ=float )
  supported_params['max_atoms_from_anom_map'] = par( desc='max. number of atoms to be picked in "one go" from anomalous diff. map', typ=int )
  supported_params['rms_threshold_anom_map'] = par( desc='RMS threshold for anom. atom picking from map', typ=float )
  supported_params['cycles_period_anom_map'] = par( desc='do anom. atom picking each N building cycles', typ=float )
  supported_params['fast_build'] = par( desc='use fast building option (default: alternate till solved); Buccaneer only)', typ=bool, share=True )
  supported_params['start_shelxe'] = par( desc='Use SHELXE in initial model building cycles then use Buccaneer', typ=bool, share=True )
  supported_params['num_parallel'] = par( desc='Number of parallel building processes', typ=int )
  supported_params['handdet'] = par( desc='Run with both hands (not set by default unless prev.steps failed to determine hand)', typ=bool )


  def TreatInOutPar(self, set_all_par=False):
    # assign the child processes to specific class attributes for convenience
    if self.GetVirtPar('target') not in ('SAD','SIRAS'):
      common.Error('Sorry, the "combined" algorithm is only available for SAD as of now ({0} not supported).'.format(
        self.GetVirtPar('target')))
    self.dmbr=self.GetOrAddProcess('dmfull')
    self.ph = self.dmbr.GetOrAddProcess('ref')
    self.dm = self.dmbr.GetOrAddProcess('dm')
    self.mb = self.GetOrAddProcess('mb')
    # check and set some defaults
    if not self.GetParam('minbigcyc'):
      self.SetParam('minbigcyc',5)
      if self.GetCrankParent() and 'phdmmb' in [p.nick for p in self.GetCrankParent().processes]:
        self.SetParam('minbigcyc',3)
    if not self.minbigcyc_init:
     self.minbigcyc_init = self.GetParam('minbigcyc')
    if not self.GetParam('maxbigcyc'):
      self.SetParam('maxbigcyc',50)
    # we can set skip_initial_build automatically if no mapcoef was inputted
    # only if this is the first subprocess though - otherwise this could not be used for -d
    if not self.GetParam('start_dm') and self.GetParam('skip_initial_build') is None and not self.inp.Get('mapcoef') and \
       self.GetCrankParent() and self in self.GetCrankParent().processes and self.GetCrankParent().processes.index(self)==0:
      self.SetParam('skip_initial_build',True)
    if not self.dmbr.GetParam('dmcyc'):
      self.dmbr.SetParam('dmcyc',7)
    if not self.dmbr.GetParam('biascyc'):
      self.dmbr.SetParam('biascyc',4)
    #if self.GetVirtPar('low_res_par') is not None and self.ph.GetVirtPar('low_res_par') is None:
    #  self.ph.SetVirtPar('low_res_par', self.GetVirtPar('low_res_par'))
    self.dmbr.SetRunDir()
    self.dmbr.TreatInOutPar()
    if not self.IsInputtedParam('ncs_det_ha'):
      self.SetParam('ncs_det_ha', self.dmbr.GetParam('ncs_det_ha'))
    if not self.IsInputtedParam('ncs_det') and self.GetParam('ncs_det_ha') or self.GetParam('ncs_det_mr'):
      self.SetParam('ncs_det', True)
    self.mb.SetRunDir()
    self.mb.TreatInOutPar()
    #if os.environ['CCP4_LIB'] and os.path.isfile(os.path.join(os.environ['CCP4_LIB'],'data','reference_structures','reference-1ajr.pdb')):
    #  self.mb.GetProg().SetKey('pdbin-ref',os.path.join(os.environ['CCP4_LIB'],'data','reference_structures','reference-1ajr.pdb'))
    #  self.mb.GetProg().SetKey('mtzin-ref',os.path.join(os.environ['CCP4_LIB'],'data','reference_structures','reference-1ajr.mtz'))
    if not self.GetParam('refcyc_finish'):
      # using the low_res_par of ref (self.low_res_par may be None)
      if self.ph.GetVirtPar('low_res_par'):
        self.SetParam('refcyc_finish',12)
      else:
        self.SetParam('refcyc_finish',8)
    if self.mb.GetProg(name='buccaneer'):
      # we start with fast and switch to best after first cycle unless in final mode already
      if self.IsTrueOrNoneParam('fast_build') and not self.mb.GetProg(name='buccaneer').IsKey('fast'):
        self.mb.GetProg(name='buccaneer').SetKey('fast')
    if not self.mb.GetProg('shelxe') and not self.IsInputtedParam('start_shelxe') and \
       self.mb.resol and self.mb.resol<2.15 and program.from_name('shelxe',None).CheckBinary(silent=True):
      self.SetParam('start_shelxe',True)
    if not self.mb.GetProg('shelxe') and self.GetParam('start_shelxe'):
      self.mb.AddProg('shelxe')
    self.dm_cyc_init=self.dmbr.GetParam('dmcyc')
    self.finish_mode = False
    process.TreatInOutPar(self,set_all_par)
    if self.GetParam('R_threshold'):
      self.R_tresh=self.GetParam('R_threshold')
    if self.GetParam('fom_threshold'):
      self.fom_tresh_max=self.GetParam('fom_threshold')
    if self.GetParam('always_exclude_free') is None:
      self.SetParam('always_exclude_free', True)
    if self.GetParam('max_atoms_from_anom_map') in (False,0):
      self.find_new = False
    if self.IsTrueOrNoneParam('rms_threshold_anom_map'):
      if self.inp.Get('model',has_atomtypes=True) and self.inp.Get('model',has_atomtypes=True).GetAtomType() in ('SE','S'):
        self.SetParam('rms_threshold_anom_map', 5.0)
      else:
        self.SetParam('rms_threshold_anom_map', 6.0)
    if self.IsTrueOrNoneParam('cycles_period_anom_map'):
      self.SetParam('cycles_period_anom_map', self.bias_reest_cyc*4)
    if self.ph.GetProg('refmac') and not self.ph.GetProg('refmac').IsKey('make'):
      self.ph.GetProg('refmac').SetKey('make',('hydrogens','no'))
      self.hydrog_no_set=True
    if not self.GetParam('num_parallel'):
      self.SetParam('num_parallel',1)

  def RunPreprocess(self,*args,**kwargs):
    process.RunPreprocess(self,*args,**kwargs)
    if not self.ph.inp.Get('model', typ=('substr','partial+substr'), has_atomtypes=True):
      common.Error('No substructure model or atomtype(s) inputted for {0}.'.format(self.name))
    # remove exclude set if it should be only used for the last cycle
    if not self.GetVirtPar('always_exclude_free') and self.inp.Get('exclude',typ='freeR'):
      self.dmbr.inp.Delete(self.inp.Get('exclude',typ='freeR'))
      self.mb.inp.Delete(self.inp.Get('exclude',typ='freeR'))
    # prepare mlhl-comb object
    self.ref_mlhl=self.AddProcess('ref')
    self.ref_mlhl.SetParam('target', 'MLHL')
    self.ref_mlhl.SetParam('cycles', 15)
    if self.ph.GetVirtPar('low_res_par'):
      self.ref_mlhl.SetParam('cycles',25)
    self.ref_mlhl.TreatInOutPar()
    if self.hydrog_no_set:
      self.ref_mlhl.GetProg('refmac').SetKey('make',('hydrogens','no'))
    # run shelxc to get the ins file for shelxe
    if self.mb.GetProg('shelxe'):
      shelxc=self.AddProg('shelxc')
      shelxc.Run()
      self.inp.Add(shelxc.out.Get('datafile',filetype='ins'))
    # set f'' for S to 0 if f'' of the "main" atom is large
    # this is to prevent issues with the built S getting close to heavy atoms and decreasing their occupancy to 0 (removing S atoms close to heavy atoms could be used as an alternative; bucaneer's known-structure would not help as it only excludes building mainchain)
    if self.ph.GetProg(supported=True).nick=='refmac' and (not self.ph.GetProg(supported=True).IsKey('ANOM') or \
        not 'FORM S' in self.ph.GetProg(supported=True).GetKey('ANOM') ):
      fp,fpp,dn,att=self.ph.inp.Get('model', typ=('substr','partial+substr'), has_atomtypes=True).Getfpfpp()
      if fpp is not None and fpp>3 and att!='S':
        self.ph.GetProg(supported=True).SetKey('ANOM','FORM S 0.0 0.0')
    # various other preparations
    self.SetParam('minbigcyc', self.GetParam('minbigcyc')*2)
    self.dmbr.separate_models=self.dmbr.AddProcess('sepsubstrprot',propagate_inp=False)
    self.ref_mlhl.separate_models=self.ref_mlhl.AddProcess('sepsubstrprot',propagate_inp=False)
    self.dmbr.logfilehandle=self.logfilehandle
    self.best_pdb_type='partial+substr'
    self.best_pdb,self.best_pdb_2=os.path.join(self.rundir,"best.pdb"),os.path.join(self.rundir,"best2.pdb")
    self.best_mtz,self.best_mtz_2=os.path.join(self.rundir,"best.mtz"),os.path.join(self.rundir,"best2.mtz")
    self.num_cyc_dm=[self.dmbr.GetParam('dmcyc'),self.dmbr.GetParam('dmcyc'),self.dmbr.GetParam('dmcyc')]
    if self.stop_file and self.GetCrankParent() and self.GetCrankParent().rundir:
      self.stop_file = os.path.join(self.GetCrankParent().rundir,self.stop_file)

  def UpdateResults(self, update, skip_output=False):
    # defining variables just for the convenience
    mb_stat=self.mb.GetProg(supported=True).GetStat
    ph_stat=self.ph.GetProg(supported=True).GetStat
    if update=='mb':
      if self.mb.GetProg(supported=True).nick=='shelxe':
        self.mb_res_all.append( (self.mb.cyc, mb_stat('res_built'), len(list(filter(None,mb_stat('res_built_per_frag')))), 0, 0) )
      else:
        self.mb_res_all.append( (self.mb.cyc, mb_stat('res_built'), mb_stat('frag_built'), 
                                 mb_stat('compl_chain'), mb_stat('compl_res')) )
    elif update=='ph':
      self.ref_res_all.append( (ph_stat('rfact')[-1], ph_stat('fom')[-1], (ph_stat('rfree',accept_none=True) or [None])[-1]) )
    if not skip_output:
      self.PrintActualLogGraph(update)

  def CopyToParallelBuild(self):
    self.mb_parallel = []
    for i in range(self.GetParam('num_parallel')):
      self.mb_parallel.append( self.AddProcessCopy(self.mb, deeper_copy=True) )
      self.mb_parallel[i].SetRunDir( os.path.join(self.rundir,self.mb.nick)+'-'+str(i), reset_subtree=True )
      # buccaneer or shelxe only supported as of now
      if self.mb_parallel[i].GetProg().nick=='buccaneer':
        self.mb_parallel[i].GetProg().SetKey('model-index',i, keep_previous=False)
        #self.mb_parallel[i].GetProg().SetKey('sequence-reliability',0.995-(i*0.005), keep_previous=False)
        #self.mb_parallel[i].GetProg().SetKey('model-filter-sigma',float(i+1),keep_previous=False)
        if os.path.isfile(self.best_mtz) and self.mb.GetProg().GetKey('no-correlation-mode') and i in (1,2):
          modelcopy=self.mb_parallel[i].GetProg().inp.AddCopy( self.mb.GetProg().inp.Get('model',typ=('partial','partial+substr')) )
          mapcopy=self.mb_parallel[i].inp.AddCopy(self.mb.inp.Get('mapcoef',typ='combined'))
          if i==1:
            #print 'trying the "best" map thus far!'
            self.mb_parallel[i].inp.SetFileToChild(mapcopy,self.best_mtz,'mtz')
            #if os.path.getsize(self.best_pdb)==os.path.getsize(modelcopy.GetFileName()) and os.path.isfile(self.best_pdb_2):
            #  self.mb_parallel[i].inp.SetFileToChild(mapcopy,self.best_mtz_2,'mtz')
            if os.path.isfile(self.best_mtz_2):
              self.mb_parallel[i].inp.SetFileToChild(mapcopy,self.best_mtz_2,'mtz')
          if i==2:
            #print 'trying no mr mode!'
            #self.mb_parallel[i].GetProg().inp.Clear('model')
            #self.mb_parallel[i].GetProg().SetKey('mr-model-seed',False,keep_previous=False)
            #print 'trying the "best" model thus far!'
            self.mb_parallel[i].GetProg().inp.SetFileToChild(modelcopy,self.best_pdb,'pdb')
            if os.path.getsize(self.best_pdb)==os.path.getsize(modelcopy.GetFileName()) and os.path.isfile(self.best_pdb_2):
              self.mb_parallel[i].GetProg().inp.SetFileToChild(modelcopy,self.best_pdb_2,'pdb')
            modelcopy.SetType(self.best_pdb_type)
          #if i==5:
        #    self.mb_parallel[i].GetProg().inp.SetFileToChild(modelcopy,self.best_pdb,'pdb')
        #    if os.path.isfile(self.best_mtz_2):
        #      self.mb_parallel[i].GetProg().inp.SetFileToChild(mapcopy,self.best_mtz_2,'mtz')
        #  if i==6:
        #    self.mb_parallel[i].GetProg().inp.SetFileToChild(modelcopy,self.best_pdb_2,'pdb')
        #    self.mb_parallel[i].GetProg().inp.SetFileToChild(mapcopy,self.best_mtz,'mtz')
        if i==3 and self.mb.cyc>1 and hasattr(self,'dmbr_parallel'):
          #print 'trying second best map+model from previous cycle!'
          self.mb_parallel[i].inp.Set(self.dmbr_parallel[self.sumr_i_sort[1]].out.mapcoef)
          self.mb_parallel[i].GetProg().inp.Set(self.dmbr_parallel[self.sumr_i_sort[1]].out.model)
        ### disabled for now due to strange crashes and no buccaneer model used with this option.
        #if i==4 and self.dmbr.out.Get('fsigf',typ='average-derived'):
          #print 'trying average-derived data!'
        #  self.mb_parallel[i].inp.Set(self.dmbr.out.Get('fsigf',typ='average-derived'))
        if i==1 and not self.mb.GetProg().GetKey('no-correlation-mode'):
          self.mb_parallel[i].SetParam('from_weighted_map', not self.mb_parallel[i].GetParam('from_weighted_map'))
        if i==2 and not self.mb.GetProg().GetKey('no-correlation-mode') and \
           self.mb_parallel[i].GetProg().GetKey('sequence-reliability'):
          self.mb_parallel[i].orig_seq_rel=self.mb_parallel[i].GetProg().GetKey('sequence-reliability')
          self.mb_parallel[i].GetProg().SetKey('sequence-reliability',1.-2.0*(1.-self.mb_parallel[i].orig_seq_rel), keep_previous=False)
        #if self.mb_parallel[i].inp.Get('mapcoef',typ='combined'):
        #  print 'map used:',i,self.mb_parallel[i].inp.Get('mapcoef',typ='combined').GetFileName()
        #if self.mb_parallel[i].GetProg().inp.Get('model',typ=('partial','partial+substr')):
        #  print 'model used:',i,self.mb_parallel[i].GetProg().inp.Get('model',typ=('partial','partial+substr')).GetFileName()
      else:
        self.mb_parallel[i].GetProg().SetArg('G',0.7+float(i)*0.05, keep_previous=False)
      aver=self.mb_parallel[i].inp.Get( 'fsigf',typ='average',col='f' )
      if not aver:
        aver=self.mb_parallel[i].inp.Get( 'fsigf',typ='average-derived',col='f' )
      mtzproc = self.mb_parallel[i].AddProcess('mergemtz',propagate_inp=False,propagate_out=False)
      mtzproc.inp.Add(aver)
      mtzproc.Run()
      avercopy=self.mb_parallel[i].inp.AddCopy(aver)
      self.mb_parallel[i].inp.SetFileToChild(avercopy,mtzproc.out.Get('datafile').GetFileName())
      self.mb_parallel[i].inp.Delete(aver)
      self.mb_parallel[i].processes.remove(mtzproc)

  def CopyToParallelDM(self):
    self.dmbr_parallel = []
    #solvent_diff=0.015
    for i in range(self.GetParam('num_parallel')):
      # create copies and assign ph/dm accordingly
      self.dmbr_parallel.append( self.AddProcessCopy(self.dmbr, deeper_copy=True) )
      if self.dmbr.nick == 'dmfull':
        self.dmbr_parallel[i].dm = self.dmbr_parallel[i].GetProcess('dm')
        self.dmbr_parallel[i].ph = self.dmbr_parallel[i].GetProcess('ref')
        # try slight variations in solvent content
        #solvent=self.dmbr.GetParam('solvent_content')
        #if not solvent and self.inp.Get(has_solvent_content=True):
        #  solvent=self.inp.Get(has_solvent_content=True).solvent_content
        #if solvent:
        #  self.dmbr_parallel[i].SetParam('solvent_content', solvent + (float(int(i+1)/2)*pow(-1,i)*solvent_diff))
        #  print 'fdsfdsfs',i,solvent + (float(int(i+1)/2)*pow(-1,i)*solvent_diff)
      self.dmbr_parallel[i].SetRunDir( os.path.join(self.rundir,self.dmbr.nick)+'-'+str(i), reset_subtree=True )
      # make sure that data is copied over to this directory (preventing rewrite issues)
      free=self.dmbr_parallel[i].inp.Get( 'exclude',typ='freeR' )
      plus=self.dmbr_parallel[i].inp.Get( 'fsigf',typ='plus',col='f' )
      minu=self.dmbr_parallel[i].inp.Get( 'fsigf',typ='minus',col='f' )
      aver=self.dmbr_parallel[i].inp.Get( 'fsigf',typ='average',col='f' )
      mtzproc = self.dmbr_parallel[i].AddProcess('mergemtz',propagate_inp=False,propagate_out=False)
      if free:
        mtzproc.inp.Set(free)
      mtzproc.inp.Set(plus), mtzproc.inp.Add(minu), mtzproc.inp.Add(aver)
      mtzproc.Run()
      if free:
        freecopy=self.dmbr_parallel[i].inp.AddCopy(free)
        self.dmbr_parallel[i].inp.SetFileToChild(freecopy,mtzproc.out.Get('datafile').GetFileName())
        self.dmbr_parallel[i].inp.Delete(free)
      pluscopy, minucopy = self.dmbr_parallel[i].inp.AddCopy(plus), self.dmbr_parallel[i].inp.AddCopy(minu)
      self.dmbr_parallel[i].inp.SetFileToChild(pluscopy,mtzproc.out.Get('datafile').GetFileName())
      self.dmbr_parallel[i].inp.SetFileToChild(minucopy,mtzproc.out.Get('datafile').GetFileName())
      avercopy=self.dmbr_parallel[i].inp.AddCopy(aver)
      self.dmbr_parallel[i].inp.SetFileToChild(avercopy,mtzproc.out.Get('datafile').GetFileName())
      self.dmbr_parallel[i].inp.Delete(aver),self.dmbr_parallel[i].inp.Delete(plus),self.dmbr_parallel[i].inp.Delete(minu)
      self.dmbr_parallel[i].processes.remove(mtzproc)
      # input the corresponding model from parallel building
      for mod in self.mb_parallel[i].out.model:
        self.dmbr_parallel[i].inp.Add(mod)

  def ParallelBuild(self, num_proc, output):
    try:
      self.mb = self.mb_parallel[num_proc]
      self.mb.Run()
      self.UpdateResults('mb', skip_output=True)
      self.mb.res_all = self.mb_res_all
      # the parent process is to be removed and attached after this routine to prevent too much copying
      self.mb.parent_process=None
      aver=self.mb.inp.Get( 'fsigf',typ='average',col='f' )
      if not aver and self.mb_parallel[0].inp.Get( 'fsigf',typ='average',col='f' ):
        self.mb.inp.Add( self.mb_parallel[0].inp.Get( 'fsigf',typ='average',col='f' ) )
    except Exception as e:
      error='Parallel building process {0} failed with error: {1}'.format(num_proc, str(e))
      output.put((-1,error))
    else:
      output.put((num_proc,self.mb))
    return True

  def ParallelDM(self, num_proc, output):
    try:
      self.dmbr = self.dmbr_parallel[num_proc]
      if self.dmbr.nick == 'dmfull':
        self.ph = self.dmbr.ph
        self.dm = self.dmbr.dm
      else:
        self.ph = self.dmbr
      self.dmbr.Run()
      if self.dmbr.nick == 'ref':
        self.UpdateResults('ph', skip_output=True)
      self.dmbr.res_all = self.ref_res_all
      # the parent process is removed and reattached after the parallel dm
      self.dmbr.parent_process, self.dmbr.logfilehandle = None, None
    except Exception as e:
      error='Parallel density modification process {0} failed with error: {1}'.format(num_proc, str(e))
      output.put((-1,error))
    else:
      output.put((num_proc,self.dmbr))

  def MergeParallel(self):
    #using a simple approach of picking the lowest R's as of now;  more advanced merges may follow.
    ignorefree=any((dmbr.res_all[-1][0]>dmbr.res_all[-1][2] for dmbr in self.dmbr_parallel))
    sumr=[]
    for i,dmbr in enumerate(self.dmbr_parallel):
      print('sum Rs:',i,dmbr.res_all[-1][0] if ignorefree else dmbr.res_all[-1][0]+dmbr.res_all[-1][2])
      sumr.append(dmbr.res_all[-1][0] if ignorefree else dmbr.res_all[-1][0]+dmbr.res_all[-1][2])
      if hasattr(self.mb_parallel[i],'orig_seq_rel') and self.mb_parallel[i].orig_seq_rel:
        self.mb_parallel[i].GetProg().SetKey('sequence-reliability',self.mb_parallel[i].orig_seq_rel)
        self.mb_parallel[i].orig_seq_rel=None
    self.sumr_i_sort=sorted(range(len(self.dmbr_parallel)), key=lambda p: sumr[p])
    min_i, min_i2 = self.sumr_i_sort[0], self.sumr_i_sort[1]
    print('min_i,min_r (min_i2,minr2)',min_i, sumr[min_i], min_i2,sumr[min_i2])
    self.processes.remove(self.dmbr)
    self.dmbr = self.AddProcessCopy( self.dmbr_parallel[min_i], deeper_copy=True )
    if self.dmbr.nick == 'dmfull':
      self.ph, self.dm = self.dmbr.ph, self.dmbr.dm
    else:
      self.ph = self.dmbr
    for o in self.dmbr.out.GetAll(stored_order=True):
      self.out.AddCopy(o)
    self.processes.remove(self.mb)
    self.mb = self.AddProcessCopy( self.mb_parallel[min_i], deeper_copy=True )
    self.ref_res_all.extend( self.dmbr.res_all[-(self.dmbr.GetParam('dmcyc') if self.dmbr.nick=='dmfull' else 1):] )
    self.mb_res_all.append( self.mb.res_all[-1] )

  def save_reset(self):
    if os.name=='nt':  # windows cannot spawn instances with non-picklable attributes
      if hasattr(self,'ccp4i2job'):  i2job,self.ccp4i2job=self.ccp4i2job,None
      else:  i2job=None
      logfh,logfhdmbr,logfhdm,logfhph,logfhmb,parent=self.logfilehandle,self.dmbr.logfilehandle,self.dm.logfilehandle,self.ph.logfilehandle,self.mb.logfilehandle,self.parent_process
      self.logfilehandle,self.dmbr.logfilehandle,self.dm.logfilehandle,self.ph.logfilehandle,self.mb.logfilehandle,self.parent_process=None,None,None,None,None,None
      return i2job,logfh,logfhdmbr,logfhdm,logfhph,logfhmb,parent
    else:
      return None,None,None,None,None,None,None

  def restore(self,data):
    i2job,logfh,logfhdmbr,logfhdm,logfhph,logfhmb,parent=data[0],data[1],data[2],data[3],data[4],data[5],data[6]
    if os.name=='nt':  # reattaching the atrributes for windows
      if hasattr(self,'ccp4i2job'):  self.ccp4i2job=i2job
      self.logfilehandle,self.dmbr.logfilehandle,self.dm.logfilehandle,self.ph.logfilehandle,self.mb.logfilehandle,self.parent_process=logfh,logfhdmbr,logfhdm,logfhph,logfhmb,parent

  def RunBody(self,*args,**kwargs):
    if not self.GetParam('handdet'):
      self.RunComb()
      self.PrintActualLogGraph(update='endparal') # rewrites the graph for paral
    else:
      windata=self.save_reset() # windows cannot spawn instances with non-picklable attributes
      import multiprocessing
      manager = multiprocessing.Manager()
      queue, queue2 = manager.Queue(), manager.Queue()
      comb_hand,self.cmb_hand=[],[]
      for i in range(2):
        comb_hand.append( multiprocessing.Process(target=self.RunCombHands, args=(i,queue,queue2)) )
        if self.GetParam('num_parallel')<=1:  # daemonic process is not allowed to have children. this means that parallel building with both hands cannot kill other processes if one fails.
          comb_hand[i].daemon=True
        comb_hand[i].start()
        self.cmb_hand.append(None)
      self.restore(windata)
      num_finished, stop_other = 0, None
      while num_finished<len(comb_hand):
        j,cmb_now,logdict = queue.get()
        if j<0:
          num_finished+=1
          if num_finished<len(comb_hand):
            common.Warning(cmb_now)
          else:
            common.Error(cmb_now)
        else:
          self.ph, self.mb, self.dmbr = cmb_now.ph, cmb_now.mb, cmb_now.dmbr
          self.ref_res_all, self.mb_res_all, self.ncs_oper = cmb_now.ref_res_all, cmb_now.mb_res_all, cmb_now.ncs_oper
          for o in cmb_now.out.GetAll(stored_order=True):  self.out.Add(o)
          if hasattr(cmb_now,'mb_parallel'):     self.mb_parallel = cmb_now.mb_parallel
          if hasattr(cmb_now,'dmbr_parallel'):   self.dmbr_parallel = cmb_now.dmbr_parallel
          self.PrintActualLogGraph(hand=j+1,**logdict)
          # stop the other hand when this hand succeeds
          if cmb_now.finish_mode and stop_other is None:
            stop_other = int(not j)
          if stop_other is j:
            queue2.put((j,'stop'))
          else:
            queue2.put((j,'continue'))
        if 'update' in logdict and logdict['update']=='end':
          num_finished+=1
          self.cmb_hand[j]=cmb_now
          self.cmb_hand[j].report_R,self.cmb_hand[j].report_Rfree=self.report_R,self.report_Rfree
          if stop_other is j:
            self.stop_message = cmb_now.stop_message
            self.PrintActualLogGraph(update='end',hand=j+1)
            self.stop_message = None
      chosen_hand,chosen_cmb = 0,self.cmb_hand[0]
      if self.cmb_hand[0] is None or (self.cmb_hand[1] and self.cmb_hand[1].R_best<self.cmb_hand[0].R_best):
        chosen_hand,chosen_cmb = 1,self.cmb_hand[1]
      self.ncs_oper,self.report_R,self.report_Rfree=chosen_cmb.ncs_oper,chosen_cmb.report_R,chosen_cmb.report_Rfree
      for o in chosen_cmb.out.GetAll(stored_order=True):  self.out.Add(o)
      if self.cmb_hand[not chosen_hand]:
        self.out2=inout.input_output(is_output=True,parent=self)
        self.out2=self.cmb_hand[not chosen_hand].out
      self.ph, self.mb, self.dmbr = chosen_cmb.ph, chosen_cmb.mb, chosen_cmb.dmbr
      self.ref_res_all, self.mb_res_all = chosen_cmb.ref_res_all, chosen_cmb.mb_res_all
      self.stop_message,self.result_str = chosen_cmb.stop_message, chosen_cmb.result_str
      self.spacegroup_change=False
      if chosen_cmb.out.Get('mapcoef').GetSpacegroupNumber(self)!=self.inp.Get('fsigf').GetSpacegroupNumber(self):
        self.spacegroup_change=chosen_cmb.out.Get('mapcoef').GetSpacegroup(self)
        self.result_str+='  Spacegroup changed from {0} to {1}.'.format(self.inp.Get('fsigf').GetSpacegroup(self),self.spacegroup_change)
      self.GetParam('refcyc_finish'),self.GetParam('maxbigcyc'),self.GetParam('cycles_period_anom_map'),self.GetParam('rms_threshold_anom_map')
      if self.cmb_hand[not chosen_hand]:
        self.ref_res_all, self.mb_res_all = self.cmb_hand[not chosen_hand].ref_res_all, self.cmb_hand[not chosen_hand].mb_res_all
        self.PrintActualLogGraph(update='endparal',hand=int(not chosen_hand)+1) # rewrites the graphs for paral other hand
        self.ref_res_all, self.mb_res_all = self.cmb_hand[chosen_hand].ref_res_all, self.cmb_hand[chosen_hand].mb_res_all
      self.PrintActualLogGraph(update='endparal2',hand=int(chosen_hand)+1) # rewrites the graph for paral chosen hand
      self.PrintActualLogGraph(update='end')  # needed here since the messages are only now in the main thread


  def RunCombHands(self,num_proc=-1,queue=None,queue2=None):
    self.num_proc=num_proc
    self.queue,self.queue2=queue,queue2
    if self.logfilehandle:
      self.logfilehandle.close()
    try:
      if num_proc>0:
        self.SetRunDir(self.rundir+'_hand'+str(num_proc+1), reset_subtree=True, change_cwd=True)
        self.ReplaceInpOut('inp','inp2')
      # the other hand is not kept for each hand subprocess.  If needed then at the very least AddProcessCopy() needs to support inp2.
      self.RemoveInpOut('inp2')
      # this might not be needed?
      #self.parent_process.Initialize(self)
      self.TreatInOutPar()
      self.RunPreprocess()
      self.RunComb()
    except Exception as e:
      error='Parallel combined building process {0} failed with error: {1}'.format(num_proc, str(e))
      queue.put((-1,error,'err'))

  def RunComb(self):
    # prepare phases for initial model building
    if self.GetVirtPar('start_dm'):
      phcmb=self.dmbr.AddProcess('phcomb')
      phcmb.SetParam('target',self.ph.GetParam('target'))
      dmbr_backup=self.dmbr.BackupAnyPars()
      self.dmbr.SetParam('dmcyc',25)
      self.dmbr.Run()
      self.dmbr.RestoreAnyPars(dmbr_backup)
      self.ph.out.Set(phcmb.out.mapcoef)
      self.ph.out.Set(phcmb.out.model)
      self.dmbr.processes.remove(phcmb)
    else:
      self.ph.out.Set(self.dmbr.inp.mapcoef)
      self.ph.out.Set(self.dmbr.inp.model)
    # separate substructure/partial if needed
    obj,N=self.ph.GetProg().inp.GetTargetObjects(self.ph.GetParam('target'),modtyp='pdb')
    if obj and obj['mod'].GetType()=='partial+substr':
      self.dmbr.separate_models.inp.Set(obj['mod'])
      self.dmbr.separate_models.Run()
    self.ph.SetParam('cycles', 1)
    # if shelxe should be used then start with it
    if self.mb.GetProg('shelxe') and not self.mb.GetProg(supported=True).nick=='shelxe':
      self.mb.AddProgCopy( self.mb.programs.pop(self.mb.programs.index(self.mb.GetProg('shelxe'))), ind=0, deeper_copy=True )
    # go with model building!
    for self.mb.cyc in range(1, self.GetParam('maxbigcyc')+1):
      # cleaning all unless we are in the very last cycle(s) (so that eg anom. diff data remain in out)
      if self.ph != self.ref_mlhl or self.GetParam('minbigcyc')>1:
        self.out.ClearAll()
      # build model
      paral_build=[]
      if not self.GetVirtPar("skip_initial_build") or self.mb.cyc>1:
        self.Info("\nStarting model building cycle {0}".format(self.mb.cyc))
        self.mb.inp.Set(self.dmbr.out.mapcoef)
        self.mb.GetProg().runname=self.mb.GetProg().name+'_cyc'+str(self.mb.cyc)
        if self.GetParam('num_parallel')>1:
          #self.lock_paral = multiprocessing.Lock()
          self.CopyToParallelBuild()
          windata=self.save_reset_Win() # windows cannot spawn instances with non-picklable attributes
          output = multiprocessing.Queue()
          for i in range(self.GetParam('num_parallel')):
            paral_build.append( multiprocessing.Process(target=self.ParallelBuild, args=(i,output)) )
            paral_build[i].daemon=True
            paral_build[i].start()
          self.restore_Win(windata)
          for p in paral_build:
            j,mb_now = output.get()
            if j<0:
              common.Error(mb_now)
            self.mb_parallel[j] = mb_now
          for j,p in enumerate(paral_build):
            self.PrintActualLogGraph('mb', j)
        else:
          self.mb.Run()
          self.UpdateResults('mb')
      else:
        self.mb.out.Set(self.dmbr.out.Get('model',typ='partial'))
      # setup and run density modification recycling (with bias reduction at specified mb cycles)
      self.dmbr.inp.Set(self.inp.fsigf)
      if self.dmbr.out.Get('model',typ='substr',filetype='pdb') or self.inp.GetAll('model',custom='ncs'):
        self.dmbr.inp.Set([self.dmbr.out.Get('model',typ='substr',filetype='pdb')]+self.inp.GetAll('model',custom='ncs'))
      else:
        self.dmbr.inp.Clear('model')
      # DM specific parameters
      if self.dmbr.nick == 'dmfull':
        if not self.GetParam('no_bias_est'):
          self.dmbr.SetParam( 'no_bias_est', bool((self.mb.cyc-1)%self.bias_reest_cyc) )
        # ncs parameters will be ignored for solomon
        if self.GetParam('ncs_det') is not False:
          self.dmbr.SetParam( 'ncs_det', self.GetParam('ncs_det') if not (self.mb.cyc-2)%self.bias_reest_cyc else False)
      # use biased or unbiased phases for mlhl.  this needs more testing?
      if self.dmbr.nick=='ref' and self.inp.Get('mapcoef',typ='best'): 
        if self.GetParam('minbigcyc')>1 and self.mb.cyc>self.cyc_mlhl+1:
          self.ph.out.AddCopy( self.inp.Get('mapcoef',typ='best') ).custom.append('mlhl')
        else:
          self.ph.out.AddCopy( self.ph.out.Get('mapcoef',typ='combined') ).custom.append('mlhl')
      self.dmbr.inp.Set(self.ph.out.Get('mapcoef',typ=('best','combined')))
      if paral_build:
        paral_dm=[]
        self.CopyToParallelDM()
        windata=self.save_reset_Win() # windows cannot spawn instances with non-picklable attributes
        output = multiprocessing.Queue()
        for i in range(self.GetParam('num_parallel')):
          paral_dm.append( multiprocessing.Process(target=self.ParallelDM, args=(i,output)) )
          paral_dm[i].daemon=True  # daemon makes sure that quitting main process kills its children processes
          paral_dm[i].start()
        self.restore_Win(windata)
        for p in paral_dm:
          j,dmbr_now = output.get()
          self.dmbr_parallel[j].logfilehandle,self.dmbr_parallel[j].parent_process=self.logfilehandle,self
          if j<0:
            common.Error(dmbr_now)
          self.dmbr_parallel[j] = dmbr_now
        # a separate loop to prevent hanging if error occurs in the printactualloggraph function
        for j,p in enumerate(paral_dm):
          self.PrintActualLogGraph('ph', j)
        self.MergeParallel()
      else:
        for mod in self.mb.out.model:
          self.dmbr.inp.Add(mod)
        self.dmbr.Run()
        # update results in case of mlhl-comb ref (dmfull updates results from itself on the fly)
        if self.dmbr.nick == 'ref':
          self.UpdateResults('ph')
      if self.finish_mode and self.dmbr.nick == 'ref' and self.ph.R<self.R_tresh-0.1:
        self.Info('Going to add waters and re-refine.')
        addwat=self.GetOrAddProcess('addwaters',propagate_out=False)
        addwat.inp.Set(self.dmbr.out.Get('model',typ=('partial','partial+substr')))
        addwat.inp.Set(self.dmbr.out.Get('mapcoef',typ='weighted'))
        try:
          addwat.Run()
        except (Exception) as e:
          common.Warning('Waters could not be added due to error: {}'.format(str(e)))
          if hasattr(sys,'exc_clear'): sys.exc_clear()
        else:
          if os.path.isfile(addwat.out.Get('model').GetFileName()):
            self.dmbr.inp.Add(addwat.out.Get('model'))
            dmc=self.dmbr.GetParam('cycles')
            self.dmbr.SetParam('cycles',3)
            self.dmbr.Run()
            self.dmbr.SetParam('cycles',dmc)
            self.UpdateResults('ph')
          else:
            common.Warning('Waters could not be added.')
      # simple (mainly FOM based) decision making
      self.decide()
      if self.GetParam('minbigcyc')<=0:
        break
      # try some simple extra actions before new bias estimation if we are not successful yet
      ###if self.mb.cyc%self.bias_reest_cyc==0:
      self.try_something()
    # this file cleanup stuff should be properly organized. at the moment just ad hoc like this.
    if os.path.isfile(self.best_mtz):
      os.remove(self.best_mtz),os.remove(self.best_pdb)
    if os.path.isfile(self.best_mtz_2):
      os.remove(self.best_mtz_2)
    if os.path.isfile(self.best_pdb_2):
      os.remove(self.best_pdb_2)
    #self.Clean(self.dm_be.mtzinitwork)
    self.WriteXMLToContinue()
    # rough FOM based result estimation and printing
    self.FOMresult()
    self.UpdateResults('end')

  def WriteXMLToContinue(self):
    continue_copy = self.from_xml(self.init_xmlelem, dummy=True)
    if not continue_copy.GetVirtPar('minbigcyc') or continue_copy.GetVirtPar('minbigcyc')<10:
      continue_copy.SetVirtPar('minbigcyc',10)
    continue_copy.inp.Add( self.ph.out.Get('mapcoef',typ=('best','combined')), propagate=False )
    if self.ph.out.Get('model',typ=('partial+substr','partial')):
      continue_copy.inp.Add( self.ph.out.Get('model',typ=('partial+substr','partial')), propagate=False )
    if self.GetCrankParent():
      crank_p=copy.copy(self.GetCrankParent())
    else:
      crank_p=process.from_name('crank',None)
    crank_p.processes=[continue_copy,]
    common.WriteXML(crank_p, os.path.join(self.rundir, self.nick+'_continue.xml'))


  def __init__(self,*args,**kwargs):
    super(comb_phdmmb,self).__init__(*args,**kwargs)
    self.fom_tresh_min = 0.38
    self.fom_tresh_min_init = self.fom_tresh_min
    self.fom_tresh_max=0.58
    self.R_tresh=0.45
    self.bias_reest_cyc=3
    #self.find_rem_heavy=find_remove_refine_heavy(ph)
    self.mb_correlcyc=2
    self.try_dm=0
    self.R_best=1.0
    self.Rfree_best=1.0
    self.report_R,self.report_Rfree=1.0,1.0
    self.cyc_best=0
    self.cyc_mlhl=0
    self.mb_res_all, self.ref_res_all = [], []
    self.find_new=1
    self.extra_cyc=0
    self.stop_message=None
    self.stop_file='stop_file'
    self.result_str=None
    self.ncs_oper=None
    self.hydrog_no_set=None
    self.minbigcyc_init = None
    #if debug:
    #  self.debug=debug
    #  self.ph.debug=debug
    #  self.dm.debug=debug
    #  self.mb.debug=debug
    # queue and 2 processes if both hands are run in parallel
    self.queue,self.queue2=None,None
    self.num_proc=-1


  def decide(self):
    """FOM and Rcomb based decision making"""
    # for convenience
    ph, dmbr, mb = self.ph, self.dmbr, self.mb
    if (hasattr(self,'ccp4i2job') and self.ccp4i2job and self.ccp4i2job.testForInterrupt()) or os.path.isfile(self.stop_file):
      self.SetParam('minbigcyc',0)
      self.stop_message='Stopping early on user request!'
      self.Info(self.stop_message)
    if self.dmbr.nick == 'dmfull':
      # increase weight of refinement as we are getting close
      if ph.R<self.R_tresh+0.07 and ph.fom>=self.fom_tresh_max-0.17 and not self.finish_mode:
        dmbr.SetParam('dmcyc', self.dm_cyc_init-1  if self.dm_cyc_init>1  else 1)
        ph.SetParam('cycles', 2)
      if ph.R<self.R_tresh+0.07 and ph.fom>=self.fom_tresh_max-0.07 and not self.finish_mode:
        dmbr.SetParam('dmcyc', self.dm_cyc_init-2  if self.dm_cyc_init>2  else 1)
        if ph.R<self.R_tresh+0.05:
          ph.SetParam('cycles', 3)
          self.mb_correlcyc=2
          if ph.R<self.R_tresh+0.025:
            ph.SetParam('cycles', 5)
            self.mb_correlcyc=3
      # switch to finishing stage
      if not self.finish_mode:
        if ph.R<self.R_tresh+0.05 and ph.fom>=self.fom_tresh_max:
          dmbr.SetParam('dmcyc', self.dm_cyc_init-3  if self.dm_cyc_init>3  else 1)
        if (ph.fom>=self.fom_tresh_max and ph.R<self.R_tresh) or (ph.fom>=self.fom_tresh_max-0.05 and ph.R<self.R_tresh-0.03):
          ph.SetParam('cycles', self.GetParam('refcyc_finish'))
          ph.SetParam('beta',1.0)
          dmbr.SetParam('biascyc', 0)
          self.mb_correlcyc=3
          # the number of finishing cycles depends on the number of cycles it took to get to the finishing stage
          self.SetParam( 'minbigcyc', min(self.minbigcyc_init+mb.cyc-1, 2*self.minbigcyc_init) )
          self.finish_mode = 1
          self.Info("Switching to finishing mode ({0} cycles).".format(self.GetParam('minbigcyc')))
          if self.hydrog_no_set:
            self.ph.GetProg('refmac').SetKey('make',False,keep_previous=False)
            self.ref_mlhl.GetProg('refmac').SetKey('make',False,keep_previous=False)
    if self.fom_tresh_min<0.5:
      self.fom_tresh_min+=0.005
    self.Info("FOM after {0}. model building cycle is {1}".format(mb.cyc,ph.fom))
    # decrease number of min cycles if we are doing very well or very badly
    if (ph.fom<self.fom_tresh_min or self.finish_mode) and not (ph.R<self.R_best-0.01 and self.GetParam('minbigcyc')<2):
      self.SetParam('minbigcyc', self.GetParam('minbigcyc')-1)
    cyc2go = min(self.GetParam('minbigcyc'), self.GetParam('maxbigcyc')-mb.cyc)
    # add extra cycles if R factor remains high in finishing mode
    if self.finish_mode and cyc2go==2 and min(self.R_best,ph.R)>0.36 and not self.extra_cyc:
      self.extra_cyc=min(self.minbigcyc_init,self.GetParam('maxbigcyc')-mb.cyc-cyc2go,10)
      if self.extra_cyc:
        self.SetParam('minbigcyc', self.GetParam('minbigcyc')+self.extra_cyc)
        self.Info('R-factor high, extra {0} cycles will be added.'.format(self.extra_cyc))
        cyc2go+=self.extra_cyc
    switch_mlhl = True if self.finish_mode and (cyc2go<=1 or mb.cyc==self.cyc_mlhl or (ph.R<self.R_tresh-0.1 and cyc2go<=4) or (ph.R<self.R_tresh-0.08 and cyc2go<=2)) and self.GetParam('maxbigcyc')>1 and self.dmbr.nick!='ref'  else False
    # "anneal" from time to time - DISABLED!
    #if mb.cyc%(self.bias_reest_cyc*4)==0 and \
    #   cyc2go>self.mb_correlcyc+(cyc2go!=(self.GetParam('maxbigcyc')-mb.cyc)) and \
    #   (ph.R>(self.R_tresh-0.08) or cyc2go>5):
    #  self.try_dm=1
    #  self.fom_tresh_min=self.fom_tresh_min_init+0.05
    # take the best model to the last cycle(s)
    if ((cyc2go==2 and switch_mlhl) or cyc2go<=1) and \
       (self.R_best<ph.R-0.01 or (self.R_best<ph.R and self.cyc_mlhl and self.cyc_best and self.cyc_best<=self.cyc_mlhl)):
      self.Info("Taking the results from cycle {0} into the last cycle(s).".format(self.cyc_best))
      shutil.copy(self.best_mtz, dmbr.out.Get('mapcoef',typ='combined').GetFileName('mtz'))
      shutil.copy(self.best_pdb, dmbr.out.Get('model',typ=self.best_pdb_type).GetFileName('pdb'))
      # run a few cycles of ref. with prev. model if something wrong happens in the very last cycle
      if self.finish_mode and cyc2go<=0 and dmbr.inp.Get('model',typ=self.best_pdb_type):
        shutil.copy(self.best_pdb, dmbr.inp.Get('model',typ=self.best_pdb_type).GetFileName('pdb'))
        dmbr.SetParam('cycles',8)
        dmbr.Run(), self.UpdateResults('ph')
    # separate substructure/partial if available
    if dmbr.out.Get('model',typ='partial+substr'):
      dmbr.separate_models.inp.Set(dmbr.out.Get('model',typ='partial+substr'))
      dmbr.separate_models.Run()
    # if only substr. was refined then copy it otherwise it might get rewritten by partial+susbtr in the next cycles
    elif dmbr.out.Get('model',typ='substr'):
      newfile=os.path.join(ph.rundir,'substr.pdb')
      shutil.copy(dmbr.out.Get('model',typ='substr').GetFileName('pdb'),newfile)
      dmbr.out.SetFileToChild(dmbr.out.Get('model',typ='substr'),newfile,'pdb')
    # switch between correlation building and standard building
    if mb.GetProg('buccaneer'):
      if mb.GetProg('buccaneer').inp.Get('model',typ=('partial','partial+substr')):
        mb.GetProg('buccaneer').inp.model.remove(mb.GetProg('buccaneer').inp.Get('model',typ=('partial','partial+substr')))
      if not mb.IsInputtedParam('from_weighted_map'):
        mb.SetVirtPar('from_weighted_map',False)
      mb.GetProg('buccaneer').SetKey('cycles',3,keep_previous=False)
      mb.GetProg('buccaneer').SetKey('no-correlation-mode',True,keep_previous=False)
      if self.out.Get('model',typ='partial'):
        mb.GetProg('buccaneer').inp.Add( self.out.Get('model',typ='partial') )
        mb.GetProg('buccaneer').SetKey('mr-model-seed',True,keep_previous=False)
        mb.SetVirtPar('from_mr_model',True)
      if  self.GetParam('fast_build',ignore_error=True) is None: #and not self.finish_mode:
        mb.GetProg('buccaneer').SetKey('fast',False,keep_previous=False)
      if self.GetParam('rebuild_only') or   ( mb.cyc>1 and \
         (mb.cyc%self.mb_correlcyc!=0 or (cyc2go<self.mb_correlcyc and self.finish_mode) or \
          (ph.R<self.R_best-0.01-0.02*int(self.cyc_mlhl==mb.cyc-1) and self.GetParam('num_parallel')<=1) ) ):
        if self.out.Get('model',typ='partial'):
          mb.GetProg('buccaneer').SetKey('cycles',2,keep_previous=False)
          mb.SetVirtPar('from_mr_model',False)
          mb.GetProg('buccaneer').SetKey('no-correlation-mode',False,keep_previous=False)
          mb.GetProg('buccaneer').SetKey('mr-model-seed',False,keep_previous=False)
          if self.GetParam('fast_build',ignore_error=True) is None and (mb.cyc%self.mb_correlcyc!=1 or self.finish_mode):
            mb.GetProg('buccaneer').SetKey('fast',True,keep_previous=False)
        if self.finish_mode and mb.cyc%self.mb_correlcyc>1 and not mb.IsInputtedParam('from_weighted_map'):
          mb.SetVirtPar('from_weighted_map')
    # if we have the best model then save it and take it as result at the end
    # usually the last model is the best but time to time this helps to keep a bit more complete model
    if ph.R<self.R_best:
      self.R_best, self.Rfree_best = ph.R, ph.Rfree
      self.cyc_best = mb.cyc
      #if self.finish_mode:
      if ph.out.Get('model',typ='partial+substr'):
        #if os.path.isfile(self.best_mtz):
        #  shutil.copy(self.best_mtz, self.best_mtz_2), shutil.copy(self.best_pdb, self.best_pdb_2)
        if self.mb.inp.Get('mapcoef',typ='combined') and self.mb.inp.Get('mapcoef',typ='combined').GetFileName('mtz')!=self.best_mtz_2:
          shutil.copy(self.mb.inp.Get('mapcoef',typ='combined').GetFileName('mtz'), self.best_mtz_2)
        if self.mb.inp.Get('model',typ=self.best_pdb_type) and self.mb.inp.Get('model',typ=self.best_pdb_type).GetFileName()!=self.best_pdb_2:
          shutil.copy(self.mb.inp.Get('model',typ=self.best_pdb_type).GetFileName('pdb'), self.best_pdb_2)
        shutil.copy(ph.out.Get('mapcoef',typ='combined').GetFileName('mtz'), self.best_mtz)
        shutil.copy(ph.out.Get('model',typ=self.best_pdb_type).GetFileName('pdb'), self.best_pdb)
    # use MLHL-comb (or return to dmfull from MLHL-comb)
    #or (mb.cyc%self.bias_reest_cyc==2 and (self.finish_mode or self.mb.cyc>self.bias_reest_cyc) and mb.cyc%(self.bias_reest_cyc*4)!=(self.bias_reest_cyc*4-1)):
    #######!!!!!!!
    #if self.finish_mode and cyc2go<=1 and self.GetParam('maxbigcyc')>1 and self.dmbr.nick!='ref':
    #if self.finish_mode and self.GetParam('maxbigcyc')>1 and self.dmbr.nick!='ref':
    if switch_mlhl:
      if cyc2go>3 and mb.cyc!=self.cyc_mlhl:
        self.cyc_mlhl = mb.cyc+3 if not self.cyc_mlhl else self.cyc_mlhl
      else:
        if self.inp.Get('exclude',typ='freeR') and not self.GetParam('always_exclude_free'):
          self.ref_mlhl.inp.Set(self.inp.Get('exclude',typ='freeR'))
        #self.mb.inp.Set(self.inp.Get('exclude',typ='freeR'))
        self.ref_mlhl.out.Add(dmbr.out.Get('model',typ='substr'))
        self.ref_mlhl.out.Add(dmbr.out.Get('model',typ='partial'))
        self.ref_mlhl.out.Add(dmbr.out.Get('mapcoef',typ='combined'))
        self.ph, self.dmbr = self.ref_mlhl, self.ref_mlhl
        self.cyc_mlhl = mb.cyc
        # run at most 3 more cycles with mlhl
        self.SetParam( 'minbigcyc', min(3, cyc2go) )
        self.Info('Switching to last {0} cycles of MLHL refinement.'.format(self.GetParam('minbigcyc')))
    ##elif self.dmbr.nick != 'dmfull':
    ##  self.dmbr = self.GetProcess('dmfull')
    ##  self.ph = self.dmbr.GetProcess('ref')
    ##  self.dmbr.out.Add(self.ref_mlhl.out.Get('model',typ='substr'))
    ##  self.dmbr.out.Add(self.ref_mlhl.out.Get('model',typ='partial'))
    ##  self.ph.out.Set(self.ref_mlhl.out.Get('mapcoef',typ=('combined','best')))
    # switch to buccaneer from shelxe 
    if self.mb.GetProg('buccaneer') and not self.mb.GetProg(supported=True).nick=='buccaneer':
      if self.mb.cyc>9 or (self.mb.cyc>4 and (ph.R<self.R_tresh+0.04 or self.mb.GetProg('shelxe').GetStat('cc')>19.)) or ph.R<self.R_tresh+0.025 or self.mb.GetProg('shelxe').GetStat('cc')>25.:
        self.mb.AddProgCopy( self.mb.programs.pop(self.mb.programs.index(self.mb.GetProg('buccaneer'))), ind=0, deeper_copy=True )
        mb.GetProg('buccaneer').inp.Add( self.out.Get('model',typ='partial') )
        mb.GetProg('buccaneer').SetKey('no-correlation-mode',False,keep_previous=False)
        self.Info('Switching from ShelxE building to Buccaneer building.')
    # just warning suppressing.
    if mb.GetProg(supported=True).nick=='shelxe':
      mb.GetVirtPar('from_weighted_map'), mb.GetVirtPar('from_mr_model')
    # minbigcyc might have changed in this routine but we want to force stop
    if (hasattr(self,'ccp4i2job') and self.ccp4i2job and self.ccp4i2job.testForInterrupt()) or os.path.isfile(self.stop_file):
      self.SetParam('minbigcyc',0)


  def FOMresult(self):
    """rough FOM/R based result estimation and its printing"""
    if self.finish_mode:
      self.result_str="Majority of model should be successfully built!"
    else:
      self.GetParam('refcyc_finish')
      if self.ph.fom<self.fom_tresh_min-0.05 or self.ph.R>self.R_tresh+0.06:
        if not self.mb_res_all[-1][1] or not self.ph.out.Get('model',typ=('partial+substr','partial')):
          self.result_str = "No residues could be traced, map too noisy. Structure solution was unsuccessful."
          common.Error(self.result_str, nosuccess=True)
        self.result_str="Wrong substructure or very weak phases suspected: the model building seems unsuccessful."
      else:
        self.result_str="A partial model was built."
    self.Info('\n'+self.result_str)


  def CreatePlotsParal(self, hand=''):
    setattr(self,'rv_res_paral'+hand,[]), setattr(self,'rv_fr_paral'+hand,[])
    setattr(self,'rv_fom_paral'+hand,[]), setattr(self,'rv_R_paral'+hand,[]), setattr(self,'rv_Rfree_paral'+hand,[])
    for i in range(self.GetParam('num_parallel')):
      #self.rv_res_paral.append( self.rv_plotmb.PlotLine(["x"],["residues paral. "+str(i)], color='#{0:02X}{1:02X}{2:02X}'.format(min(75+6*(i+3),255),min(178+6*(i+3),255),min(197+6*(i+3),255))) )
      getattr(self,'rv_res_paral'+hand).append( self.rv_plotmb.PlotLine(["x"],["residues "+str(i),""], color='#{0:02X}{1:02X}{2:02X}'.format(min(75+6*(i+3),255),min(178+6*(i+3),255),min(197+6*(i+3),255))) )
      getattr(self,'rv_fom_paral'+hand).append( self.rv_plotref.PlotLine(["x"],["FOM "+str(i)], color='#{0:02X}{1:02X}{2:02X}'.format(min(75+6*(i+3),255),min(178+6*(i+3),255),min(197+6*(i+3),255))) )
    for i in range(self.GetParam('num_parallel')):
      getattr(self,'rv_fr_paral'+hand).append( self.rv_plotmb.PlotLine(["x"],["fragments "+str(i)], color='#{0:02X}{1:02X}{2:02X}'.format(min(250+6*(i+3),255),min(162+6*(i+3),255),min(40+6*(i+3),255))) )
      getattr(self,'rv_R_paral'+hand).append( self.rv_plotref.PlotLine(["x"],["Rcomb "+str(i)], color='#{0:02X}{1:02X}{2:02X}'.format(min(250+6*(i+3),255),min(162+6*(i+3),255),min(40+6*(i+3),255))) )
    for i in range(self.GetParam('num_parallel')):
      getattr(self,'rv_Rfree_paral'+hand).append( self.rv_plotref.PlotLine(["x"],["Rfree "+str(i)], color='#{0:02X}{1:02X}{2:02X}'.format(min(153+6*(i+3),255),min(50+6*(i+3),255),min(204+6*(i+3),255))) )

  def CreatePlots(self, paral=None, hand=''):
    if hasattr(self,'rv_plotref') and hasattr(self,'rv_plotref_prev') and self.rv_plotref and self.rv_plotref_prev:
      # if both hands' plots exist already then do nothing
      return
    handstr = '' if not hand else ' - hand 2'
    if not hasattr(self,'rv_plotref') or self.rv_plotref is None: # hand 1
      self.rv_plotref = self.rv_report.Plot( 'Rcomb and FOM vs refinement cycle', \
             "Cycle", "Rcomb-factor and FOM", block="Refinement statistics", legendloc='nw', intx=True )
      self.rv_fom = self.rv_plotref.PlotLine(["x"],["FOM"])
      self.rv_R = self.rv_plotref.PlotLine(["x"],["Rcomb"])
      self.rv_Rfree = self.rv_plotref.PlotLine(["x"],["Rfree"], color='darkorchid')
      self.rv_plotmb = self.rv_plotref.parent.parent.Plot( 'Residues and fragments vs build cycle', "Cycle", \
            "Number of residues and fragments", block="Model building statistics", legendloc='nw', intx=True )
      self.rv_res = self.rv_plotmb.PlotLine(["x"],["residues"])
      self.rv_fr = self.rv_plotmb.PlotLine(["x"],["fragments"])
      if self.GetParam('num_parallel')>1 and paral is not False:
        self.CreatePlotsParal('')
    if handstr: # hand 2
      self.rv_plotref_prev = self.rv_plotref
      self.rv_plotref = self.rv_plotref_prev.parent.Plot( 'Rcomb and FOM vs refinement cycle'+handstr, \
             "Cycle", "Rcomb-factor and FOM", legendloc='nw', intx=True )
      self.rv_fom2 = self.rv_plotref.PlotLine(["x"],["FOM"])
      self.rv_R2 = self.rv_plotref.PlotLine(["x"],["Rcomb"])
      self.rv_Rfree2 = self.rv_plotref.PlotLine(["x"],["Rfree"], color='darkorchid')
      self.rv_plotmb_prev = self.rv_plotmb
      self.rv_plotmb = self.rv_plotmb_prev.parent.Plot( 'Residues and fragments vs build cycle'+handstr, "Cycle", \
            "Number of residues and fragments", legendloc='nw', intx=True )
      self.rv_res2 = self.rv_plotmb.PlotLine(["x"],["residues"])
      self.rv_fr2 = self.rv_plotmb.PlotLine(["x"],["fragments"])
      if self.GetParam('num_parallel')>1 and paral is not False:
        self.CreatePlotsParal(hand)

  def PrintActualLogGraph(self, update=None, paral_num=None, hand=-1):
    # not nedeed: rvapi and conversions are not threadsafe and called from parallel jobs, thus lock is used.
    #if hasattr(self,'lock_paral'):
    #  self.lock_paral.acquire()
    if self.queue and self.num_proc>=0:
      # in case of both hands, passing to/from the main comb process
      # the queues, parent process and logfilehandles are removed and reattached after putting to the queue
      queue,queue2,parent = self.queue,self.queue2,self.parent_process
      if hasattr(self,'ccp4i2job'):  i2job,self.ccp4i2job=self.ccp4i2job,None
      logfh,phlogfh,dmlogfh,dmbrlogfh,mblogfh,dmflogfh = self.logfilehandle,self.ph.logfilehandle,self.dm.logfilehandle,self.dmbr.logfilehandle,self.mb.logfilehandle,self.GetProcess('dmfull').logfilehandle
      mapseg = self.GetProcess('dmfull').GetProcess('segmentmap') 
      if mapseg:
        self.GetProcess('dmfull').processes.remove(mapseg)
      processes=self.processes[:]
      self.queue,self.queue2,self.parent_process = None,None,None
      self.logfilehandle,self.ph.logfilehandle,self.dm.logfilehandle,self.dmbr.logfilehandle,self.mb.logfilehandle,self.GetProcess('dmfull').logfilehandle = None,None,None,None,None,None
      self.processes=[p for p in processes if p in (self.dmbr,self.ph,self.mb,self.dm)]
      queue.put((self.num_proc,self,{'update':update,'paral_num':paral_num}))
      self.processes=processes
      if mapseg:
        self.GetProcess('dmfull').processes.append(mapseg)
      self.queue,self.queue2,self.parent_process = queue,queue2,parent
      if hasattr(self,'ccp4i2job'):  self.ccp4i2job=i2job
      self.logfilehandle,self.ph.logfilehandle,self.dm.logfilehandle,self.dmbr.logfilehandle,self.mb.logfilehandle,self.GetProcess('dmfull').logfilehandle = logfh,phlogfh,dmlogfh,dmbrlogfh,mblogfh,dmflogfh
      numproc,check_stop=self.queue2.get()
      if check_stop=='stop' and self.num_proc==numproc and not self.stop_message:
        self.SetParam('minbigcyc',0)
        self.stop_message = "Hand {} building stopped since the other hand succeeded.".format(numproc+1)
        self.Info(self.stop_message)
      return
    mb_res = self.mb_res_all if paral_num is None  else self.mb_parallel[paral_num].res_all
    ph_res = None
    if update=='ph':
      ph_res = self.ref_res_all if paral_num is None  else self.dmbr_parallel[paral_num].res_all
    if self.opened_loggraph and hand<2 and self.GetLogGraphHandle():  # loggraph for hand 2 not printed - to be removed completely at some point?
      self.GetLogGraphHandle().seek(0,0)
      self.LGInfo(self.GetCCP4Header())
      self.LGInfo(self.GetLogGraphPrograms())
      if ph_res:
        self.LGInfo('\n $TABLE : Rcomb-factor and FOM per cycle:')
        self.LGInfo('$GRAPHS    : Rcomb-factor and FOM per cycle :A:1,2,3:\n$$')
        self.LGInfo('Cycle  Rcomb-factor   FOM  $$ $$')
        for cyc,(R,fom,Rfree) in enumerate(ph_res):
          self.LGInfo('{0:3}  {1:8.3f}  {2:8.3f}'.format(cyc+1,R,fom))
        self.LGInfo('$$\n')
      if mb_res:
        self.LGInfo('\n $TABLE : Built per building cycle:')
        self.LGInfo('$GRAPHS    : Residues and fragments per cycle :A:1,2,3:\n$$')
        self.LGInfo('Cycle Residues Fragments  $$ $$')
        for cyc,res,fr,c1,c2 in mb_res:
          self.LGInfo('{0:3}  {1:8}  {2:8}'.format(cyc,res,fr))
        self.LGInfo('$$\n')
      if update=='end' and self.result_str:
        self.LGInfo('\n$TEXT:Result: $$ Final result $$')
        self.LGInfo("{0}\n$$\n".format(self.result_str))
    if self.rv_report is not None:
      if not hasattr(self,'rv_res') or (hand==2 and not hasattr(self,'rv_res2')):
        self.CreatePlots(hand='2' if hand==2 else '')
      if hand<=1:
        rv_res, rv_fr, rv_fom, rv_R, rv_Rfree = self.rv_res, self.rv_fr, self.rv_fom, self.rv_R, self.rv_Rfree
      else:
        rv_res, rv_fr, rv_fom, rv_R, rv_Rfree = self.rv_res2, self.rv_fr2, self.rv_fom2, self.rv_R2, self.rv_Rfree2
      if paral_num is not None:
        handstr='2' if hand==2 else ''
        rv_res_paral, rv_fr_paral = getattr(self,'rv_res_paral'+handstr), getattr(self,'rv_fr_paral'+handstr)
        rv_fom_paral, rv_R_paral, rv_Rfree_paral = getattr(self,'rv_fom_paral'+handstr), getattr(self,'rv_R_paral'+handstr), getattr(self,'rv_Rfree_paral'+handstr)
      if mb_res and update=='mb':
        if paral_num is not None:
          rv_res_paral[paral_num].Reset(),  rv_fr_paral[paral_num].Reset()
          cyc,res,fr,c1,c2 = mb_res[-1]
          rv_res_paral[paral_num].Data( cyc, res, int )
          rv_fr_paral[paral_num].Data( cyc, fr, int )
        if self.mb_res_all:
          cyc,res,fr,c1,c2 = self.mb_res_all[-1]
          if not paral_num: # ie always the first one
            rv_res.Data( cyc, res, int ),  rv_fr.Data( cyc, fr, int )
          if paral_num is not None:
            rv_res_paral[paral_num].Data( cyc, res, int )
            rv_fr_paral[paral_num].Data( cyc, fr, int )
        if paral_num is None or paral_num==self.GetParam('num_parallel')-1: # ie the last one
          crvapi.Flush()
      if ph_res and update=='ph':
        if self.out.Get('model',typ=('partial','partial+substr')):
          if paral_num is None:
            # this displays the actual pdb+mtz updated cycle by cycle
            # the files are continuously rewritten which might cause issues - be careful
            out_pdb = self.out.Get('model',typ=('partial','partial+substr')).GetFileName()
            if not out_pdb:  # this can happen if skipping the initial building cycle
              out_pdb = self.out.Get('model').GetFileName()
            out_mtz = crvapi.SetMtzFWT(self,self.out.Get('mapcoef',typ='combined'))
            # conversions to map are not threadsafe and called from parallel jobs, thus lock is used.
            #if hasattr(self,'lock_paral'):
            #  self.lock_paral.acquire()
            #out_map = self.out.Get('mapcoef',typ='weighted',filetype='map').GetFileName('map')
            #out_dmap = self.out.Get('mapcoef',typ='diff',filetype='map',conv_opts=['no_rewrite']).GetFileName('map')
            #out_dmap = self.out.Get('mapcoef',typ='diff',filetype='map',conv_opts=['outfilename','diff.map']).GetFileName('map')
            #if hasattr(self,'lock_paral'):
            #  self.lock_paral.release()
          if (not self.parent_process or not self.parent_process.ccp4i2) and paral_num is None and \
              (not hasattr(self,'rv_files') or (hand==2 and not hasattr(self,'rv_files2'))):
            handtext = ' hand {}'.format(hand) if hand>0 else ''
            handstr = '' if hand<=1 else '2'
            self.rv_report.Text("<BR><i><small>Output{} (updating after each refinement cycle):</small></i>".format(handtext))
            setattr(self,'rv_files'+handstr, self.rv_report.DataFile(out_pdb, "xyz", "<small>Built model and map</small>") )
            getattr(self,'rv_files'+handstr).DataFile(out_mtz, "hkl:map", flush=True)
            #getattr(self,'rv_files'+handstr).DataFile(out_map, "hkl:ccp4_map")
            #getattr(self,'rv_files'+handstr).DataFile(out_dmap, "hkl:ccp4_dmap", flush=True)
        if paral_num is not None:
          rv_fom_paral[paral_num].Reset(), rv_R_paral[paral_num].Reset(), rv_Rfree_paral[paral_num].Reset()
          #self.num_cyc_dm = self.dmbr.GetParam('dmcyc') if self.dmbr.nick=='dmfull' else 1
          # display 'best' for the previous cycles (once)
          if paral_num==0 and self.ref_res_all:
            for i in reversed(range(self.num_cyc_dm[hand])):
              R,fom,Rfree = self.ref_res_all[-1-i]
              if Rfree:
                rv_Rfree.Data( len(self.ref_res_all)-i, Rfree, int )
              rv_fom.Data( len(self.ref_res_all)-i, fom, int )
              rv_R.Data( len(self.ref_res_all)-i, R, int )
          self.num_cyc_dm[hand] = self.dmbr.GetParam('dmcyc') if self.dmbr.nick=='dmfull' else 1
          # display the starting point for paral as the last 'best' from the previous cycle
          if self.ref_res_all:
            R,fom,Rfree = self.ref_res_all[-1]
            if Rfree:
              rv_Rfree_paral[paral_num].Data( len(self.ref_res_all), Rfree, int )
            rv_fom_paral[paral_num].Data( len(self.ref_res_all), fom, int )
            rv_R_paral[paral_num].Data( len(self.ref_res_all), R, int )
          # display paral for the actual cycles
          for i in reversed(range(self.num_cyc_dm[hand])):
            R,fom,Rfree = ph_res[-1-i]
            if Rfree:
              rv_Rfree_paral[paral_num].Data( len(ph_res)-i, Rfree, int )
            rv_fom_paral[paral_num].Data( len(ph_res)-i, fom, int )
            rv_R_paral[paral_num].Data( len(ph_res)-i, R, int )
          if paral_num==self.GetParam('num_parallel')-1: # ie the last one
            crvapi.Flush()
        else:
          R,fom,Rfree = self.ref_res_all[-1]
          if Rfree:
            rv_Rfree.Data( len(self.ref_res_all), Rfree, int )
          rv_fom.Data( len(self.ref_res_all), fom, int )
          rv_R.Data( len(self.ref_res_all), R, int, flush=True )
      if update=='ncs':
        self.rv_report.Text("<BR><i><small>NCS detection:</small></i>")
        handstr = '' if hand<=0 else ' for hand {}'.format(hand)
        self.rv_report.Text("&emsp;{0} NCS operators detected{1}.".format(self.ncs_oper,handstr), flush=True)
      if update.startswith('endparal') and self.GetParam('num_parallel')>1:
        #self.rv_plotref.parent.Remove()
        if update=='endparal':
          self.rv_plotref.parent.parent.Remove()
          self.rv_plotref,self.rv_plotref_prev,self.rv_plotmb,self.rv_plotmb_prev=None,None,None,None
        self.CreatePlots(paral=False,hand='2' if hand==2 else '')
        if hand<=1:
          rv_res, rv_fr, rv_fom, rv_R, rv_Rfree = self.rv_res, self.rv_fr, self.rv_fom, self.rv_R, self.rv_Rfree
        else:
          rv_res, rv_fr, rv_fom, rv_R, rv_Rfree = self.rv_res2, self.rv_fr2, self.rv_fom2, self.rv_R2, self.rv_Rfree2
        for cyc,res,fr,c1,c2 in self.mb_res_all:
          rv_res.Data( cyc, res, int ),  rv_fr.Data( cyc, fr, int, flush=True )
        for cyc,(R,fom,Rfree) in enumerate(self.ref_res_all):
          if Rfree:
            rv_Rfree.Data( cyc+1, Rfree, int )
          rv_fom.Data( cyc+1, fom, int ),  rv_R.Data( cyc+1, R, int )
      if update=='end':
        self.report_R = self.R_best  if self.R_best<1.0  else self.ref_res_all[-1][0]
        self.report_Rfree = self.Rfree_best  if self.Rfree_best and self.Rfree_best<1.0  else self.ref_res_all[-1][2]
        if self.stop_message:
          self.rv_report.Text("<i><small>Stop:</small></i>")
          self.rv_report.Text('&emsp;'+self.stop_message)
        if self.result_str:
          self.rv_report.Text("<i><small>Result:</small></i>")
          self.rv_report.Text('&emsp;'+self.result_str, flush=True)
        #self.rv_report.Text("<BR><i><small>Output:</small></i>")
        #self.rv_files = self.rv_report.DataFile(out_pdb, "xyz", "<small>Built model and map</small>")
        #self.rv_files.DataFile(out_mtz, "hkl:map", flush=True)
    #if hasattr(self,'lock_paral'):
    #  self.lock_paral.release()

  def try_something(self):
    # try searching for new atoms and removing low occupancy ones 
    # restricted to a few times in the building due to too many related fake atoms could freeze NCS detection by Parrot etc
    if self.find_new and (self.mb.cyc<=self.bias_reest_cyc or self.mb.cyc%(self.GetParam('cycles_period_anom_map'))==0) and \
       self.ph.GetProg().out.Get('mapcoef',typ='anom-diff'):
      # find new atoms
      atomsfrommap=self.AddProcess('atomsfrommap')
      atomsfrommap.SetParam('occupancy', 0.15)
      if self.GetParam('rms_threshold_anom_map'):
        atomsfrommap.SetParam('rms_threshold', self.GetParam('rms_threshold_anom_map'))
      if self.GetParam('max_atoms_from_anom_map'):
        atomsfrommap.SetParam('max_new_atoms', self.GetParam('max_atoms_from_anom_map'))
      atomsfrommap.inp.Add(self.ph.GetProg().out.Get('mapcoef',typ='anom-diff'))
      atomsfrommap.inp.Add(self.dmbr.out.Get('model',typ='substr',filetype='pdb'))
      atomsfrommap.Run()
    # remove low occ atoms
    if self.out.Get('model',typ='substr',filetype='pdb'):
      pdbcur=self.AddProg('pdbcur')
      pdbcur.runname=pdbcur.name+'_removelowocc'
      pdbcur.inp.Add(self.out.Get('model',typ='substr',filetype='pdb'))
      with open(self.out.Get('model',typ='substr',filetype='pdb').GetFileName('pdb')) as f:
        occ = pdbcur.GetStatGrep('occup',from_str=f.read())
      maxocc = max(occ) if occ else 0.5
      pdbcur.SetKey('cutocc',0.14*maxocc)
      pdbcur.Run()
      self.ph.GetProg().out.Add(pdbcur.out.model[-1])
      self.programs.remove(pdbcur)
    #self.dmbr.out.Add(self.out.Get('model',typ='substr',filetype='pdb'))
    #if not self.finish_mode and self.mb.cyc==2*self.bias_reest_cyc:
    #  # try increasing the DM filtering radius
    #  if hasattr(self.dmbr,'dm') and self.dmbr.dm.GetProg().nick=='parrot':
    #    self.dmbr.dm.GetProg().SetKey( 'solvent-mask-filter-radius', self.dmbr.dm.GetProg().GetStat('radius_auto')*2 )
    # try including standard DM if DM2 might got stuck - DISABLED!!!!
    if self.try_dm:
      self.Info("Phasing-DM annealing will follow.")
      # first run 10 cycles of current model refinement without phase combination
      ph_backup=self.ph.BackupAnyPars()
      self.ph.SetParam('cycles',10)
      self.ph.SetParam('phcomb',False)
      self.ph.SetParam('beta',1.0)
      self.ph.GetProg().runname=self.ph.GetProg().nick+'_ref'
      self.ph.Run()
      self.ph.RestoreAnyPars(ph_backup)
      # now run standard DM
      phcomb=self.dmbr.AddProcess('phcomb')
      phcomb.SetParam('target',self.ph.GetParam('target'))
      phcomb.SetParam('cycles',0)
      phcomb.AddProg(self.ph.GetProg().nick)
      dmbr_backup=self.dmbr.BackupAnyPars()
      self.dmbr.SetParam('dmcyc',15)
      self.dmbr.SetParam('biascyc',4)
      self.dmbr.SetParam('no_bias_est',False)
      # try solvent optimization - unless the user inputted solvent content or number of monomers
      mon_obj=self.dmbr.inp.Get(has_monomers_asym=True)
      solv_obj=self.dmbr.inp.Get(has_solvent_content=True)
      if self.mb.cyc<=5*self.bias_reest_cyc and self.GetParam('optimize_solvent') is not False and \
         not self.dmbr.IsInputtedParam('solvent_content') and \
         (not mon_obj or not mon_obj.IsInputtedAttr('monomers_asym')) and \
         (not solv_obj or not solv_obj.IsInputtedAttr('solvent_contest')):
        self.dmbr.SetParam('optimize_solvent')
        self.dmbr.SetParam('threshold_stop',1.)
      #self.dmbr.SetParam('ncs_det', None)
      self.dmbr.Run()
      self.ph.out.Set(phcomb.out.mapcoef,propagate=False)
      self.dmbr.processes.remove(phcomb)
      self.dmbr.RestoreAnyPars(dmbr_backup)
      self.try_dm=0
