#!/usr/bin/python

import multiprocessing
import os
import shutil

import gemmi
import numpy as np

from .. import common, crvapi, inout
from ..process import process


par=common.parameter

class dmfull(process):
  name="density modification with Fourier recycling"
  short_name="density modification"
  supported_progs=["parrot", "shelxe"]
  supported_procs=["dm", "phcomb", "ref"]
  supported_params = {}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True, share=True )
  supported_params['dmcyc'] = par( desc='Number of density modification iterations', typ=int )
  supported_params['biascyc'] = par( desc='Number of iterations used for bias estimation', typ=int )
  supported_params['beta'] = par( desc='Bias parameter beta', typ=float, share=True )
  supported_params['solvent_content'] = par( desc='(Expected) solvent fraction of crystal', typ=float, share=True )
  supported_params['no_bias_est'] = par( desc='Skip bias estimation', typ=bool )
  supported_params['ncs_det'] = par( desc='Determine NCS from heavy atoms or partial model (Parrot only)', typ=bool )
  supported_params['ncs_det_ha'] = par( desc='Determine NCS from heavy atoms (Parrot only)', typ=bool, share=True )
  supported_params['ncs_det_mr'] = par( desc='Determine NCS from partial model (Parrot only)', typ=bool, share=True )
  supported_params['solventmask_radius'] = par( desc='Use the specified solvent mask radius (Parrot only)', typ=float, share=True )
  supported_params['optimize_solvent'] = par( desc='Try to optimize solvent content fraction', typ=bool )
  supported_params['incr_radius_threshold_fom'] = par( desc='Filtering radius is increased if FOM after cycle 5 is below this value (Parrot only)', typ=float )
  supported_params['threshold_stop'] = par( desc='Stop earlier if this FOM threshold is reached', typ=(float,bool) )
  supported_params['handdet'] = par( desc='Run with both hands (not set by default unless prev.steps failed to determine hand)', typ=bool )
  supported_params['solvent_perturb'] = par( desc='Use solvent perturbation by the specified solvent fraction difference', typ=(float,bool), share=True )
  supported_params['map_segmentation'] = par( desc='Use map solvent segmentation from deep learning', typ=bool, share=True )


  def Init(self):
    self.stopped_early=False
    self.report_fom=0.0
    # queue and 2 processes if both hands are run in parallel
    self.queue=None
    self.num_proc=-1
    self.phfile_end_sep='_'
    self.guess=None
    self.ph=None
    self.dm=None
    self.dmcyc_base=25
    self.ncs_opers,self.ncs_oper_correls,self.ncs_cyc=[],[],-1

  def TreatInOutPar(self, set_all_par=False):
    # defaults if no programs and no subprocesses are specified
    if not self.GetProg(supported=True) and not self.GetProcess(supported=True):
      if self.GetVirtPar('target') in ('SAD','SIRAS','MLHL'):
        self.AddProcess('dm')
        self.AddProcess('phcomb')
      else:
        self.AddProg('parrot')
    if self.GetProg(supported=True):
      # any complete DM with recycling within a single program stuff shall go here
      if not self.IsInputtedParam('no_bias_est'):
        self.SetParam('no_bias_est', True)
    else:
      # parameters specific for our multivariate and/or bias reduced density modification recycling
      # many of these parameters are our DM specific, so here we also set their defaults
      # now first assign the child processes to specific class attributes for convenience
      # and set defaults if an incomplete subprocesses list is specified:
      self.dm = self.GetOrAddProcess('dm')
      self.ph = self.GetProcess('phcomb')
      if not self.ph:
        self.ph = self.GetProcess('ref')
      if not self.ph:
        self.ph = self.AddProcess('phcomb')
      if self.ph.nick=='ref':
        self.ph.SetVirtPar('phcomb')
      # check and set some defaults
      solv_cont = self.GetParam('solvent_content')
      if not solv_cont:
          solv_obj = self.inp.Get(has_solvent_content=True)
          if solv_obj:
            solv_cont=solv_obj.solvent_content
      if self.GetParam('dmcyc') is None:
        self.SetParam('dmcyc',self.dmcyc_base)
        if self.GetParam('target') in ('MLHL','MAD'):
          self.SetParam('dmcyc',15)
        self.SetDMcycSolv(solv_cont)  # this can change with mapseg after its solv. det.
        if self.parent_process and self.parent_process.nick=='handdet':
          self.SetParam('dmcyc',10)
      # bias estimation parameters
      if self.GetParam('biascyc') is None:
        self.SetParam('biascyc',5)
      if self.GetParam('no_bias_est') is None and (self.GetParam('beta') \
         or (self.parent_process and self.parent_process.nick in ('handdet','mbref'))):
        self.SetParam('no_bias_est',True)
      if not self.IsInputtedParam('ncs_det_ha'):
        if self.GetParam('ncs_det') is not False and not self.GetParam('ncs_det_mr'):
          self.SetParam('ncs_det_ha', True)
        # real time decision making - disabled due to problems with how to present in the GUI...
        #if self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True,filetype='pdb'):
        #  # determine the actual number of atoms in the pdb
        #  pdbcur=self.AddProg('pdbcur', propagate_out=False)
        #  for at in self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True).GetAtomTypes():
        #    pdbcur.SetKey('delatom', "{0}[{0}]:*".format(at))
        #  pdbcur.Run()
        #  if pdbcur.GetStat('atoms_deleted')>35:
        #    self.SetParam('ncs_det_ha', False)
        #  self.programs.remove(pdbcur)
      # disable ncs detection from heavy atoms by default if it would likely take too long
      # this should be a reasonable solution in the case of substructure detection...
        if self.GetParam('ncs_det') is not True and self.inp.Get(typ='substr',has_num_atoms=True) and \
           self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms>60:
          if self.GetParam('ncs_det_ha') is not False and (not self.dm.GetProg(supported=True) or self.dm.GetProg(supported=True).nick=='parrot'):
            if self.logfilehandle and not self.logfilehandle.closed:
              self.Info('Many expected substructure atoms: NCS detection by Parrot would be likely slow and wont be performed.', stdout=False)
            self.SetParam('ncs_det_ha', False)
      if self.GetParam('incr_radius_threshold_fom') is None:
        self.SetParam('incr_radius_threshold_fom', 0.49)
      if self.GetParam('threshold_stop') is None:
        self.SetParam('threshold_stop', 0.56)
      if self.GetParam('solvent_perturb') is None:
        #if self.GetParam('map_segmentation'):
        #  self.SetParam('solvent_perturb', 0.0)
        #elif solv_cont:
          self.SetParam('solvent_perturb', 0.05)#-0.1*abs(solv_cont-0.5))
      self.ph.TreatInOutPar(set_all_par)
      self.dm.TreatInOutPar(set_all_par)
      # this will have to be enabled for free lunch or other extension schemes
      #self.dm.GetProg().keep_nodata=True
    process.TreatInOutPar(self,set_all_par)

  def SetDMcycSolv(self,solv_cont):
    self.SetParam('dmcyc',int(self.dmcyc_base))
    if solv_cont and solv_cont<0.45:
      self.SetParam('dmcyc',int(self.GetParam('dmcyc')*4/5))
    if solv_cont and solv_cont<0.4:
      self.SetParam('dmcyc',int(self.GetParam('dmcyc')*4/5))
    if solv_cont and solv_cont<0.35:
      self.SetParam('dmcyc',int(self.GetParam('dmcyc')*4/5))

  def RunBody(self,*args,**kwargs):
    # if program is to be run directly (eg "original" parrot, shelxe)
    if any(p for p in self.programs if p.nick in self.supported_progs):
      process.RunBody(self)
    # dm program and phase comb program combination
    else:
      # make sure there is container with input data mtz file in self.inp 
      # (copy it from program input if needed) - ie the dmfull input data mtz is saved here and used later
      if self.dm.GetProg().inp.Get('fsigf',typ='average',col='f'):
        fsf = self.inp.AddCopy( self.dm.GetProg().inp.Get('fsigf',typ='average',col='f'), propagate=False )
      else:
        common.Error('(Average) F/SIGF is not available for input to {0}'.format(self.name))
      if self.dm.GetProg().nick=='prasa':
        ecalc=self.AddProg('ecalc',propagate_out=False)
        ecalc.inp.Set(fsf)
        ecalc.Run()
        self.inp.AddCopy(ecalc.out.Get('fsigf',typ='average'))
      self.inp.SetCopy( self.dm.GetProg().inp.GetAll('mapcoef',stored_order=True), propagate=False )
      self.CreateWorkingSetMTZ(new_seed=True)
      if self.GetVirtPar('optimize_solvent'):
        dmf_copy=self.GetProcessCopy()

      if not self.GetParam('handdet'):
        self.RunComb()
      else:
        if os.name=='nt':  # windows cannot spawn instances with non-picklable attributes
          if hasattr(self,'ccp4i2job'):  i2job,self.ccp4i2job=self.ccp4i2job,None
          logfh,parent=self.logfilehandle,self.parent_process
          self.logfilehandle,self.parent_process=None,None
        manager = multiprocessing.Manager()
        queue = manager.Queue()
        dm_hand,self.dmf=[],[]
        for i in range(2):
          dm_hand.append( multiprocessing.Process(target=self.RunCombHands, args=(i,queue)) )
          dm_hand[i].daemon=True
          dm_hand[i].start()
          self.dmf.append(None)
        if os.name=='nt':  # reattaching the atrributes for windows
          if hasattr(self,'ccp4i2job'):  self.ccp4i2job=i2job
          self.logfilehandle,self.parent_process=logfh,parent
        num_finished=0
        while num_finished<len(dm_hand):
          j,dm_now,logdict = queue.get()
          if j<0:
            num_finished+=1
            common.Error(dm_now)
          else:
            self.ph=dm_now.ph
            self.out=dm_now.out
            self.PrintActualLogGraph(hand=j+1,**logdict)
          if 'finished' in logdict and logdict['finished']:
            num_finished+=1
            self.dmf[j]=dm_now
            self.dmf[j].report_fom=self.report_fom
        self.out=self.dmf[0].out
        self.out2=inout.input_output(is_output=True,parent=self)
        self.out2=self.dmf[1].out
        self.GetParam('biascyc'),self.GetParam('threshold_stop'),self.GetParam('incr_radius_threshold_fom'),self.GetParam('dmcyc')
        self.report_fom = max(self.dmf[0].report_fom,self.dmf[1].report_fom)

      if self.GetVirtPar('optimize_solvent'):
        self.OptimizeSolvent(dmf_copy)
      #self.Clean(mtz_in_saved)
      #if dm_be is not None:
      #  self.Clean(dm_be.mtzinitwork)

  def OptimizeSolvent(self, dmf_copy):
    #tries to optimize solvent by FOM based multiple jobs with different solvent content
    #requires a dmfull object on input which is assumed to be run already (thus having self.fom*)
    # first make sure all is available and consistent
    if not self.dm.inp.Get('sequence'):
      common.Warning('Solvent optimization will be skipped as sequence not known.')
      return
    if self.dm.inp.Get(has_monomers_asym=True):
      num_mon = self.dm.inp.Get(has_monomers_asym=True).monomers_asym
    else:
      matthews=dmf_copy.GetProcess(self.dm.nick).GetOrAddProcess('matthews')
      matthews.out.never_propagate=True
      matthews.Run()
      # this assumes gcx; matthews_coef may be added
      num_mon = matthews.GetProg().GetStat('monomers_asym')
      if self.dm.GetParam('solvent_content') and \
        abs(matthews.GetProg().GetStat('solvent_content'),self.dm.GetParam('solvent_content'))>0.1:
          common.Warning('Solvent optimization will be skipped as inputted solvent content differs too much from the estimation.')
          return
    # now try the optimization
    num_mon_start=num_mon
    fom_best=sum(self.ph.fom_all[5:10])
    dmf_best=self
    num_mon_change=-1
    if num_mon<=1:
      num_mon_change=1
    while num_mon_change!=0:
      num_mon = num_mon+num_mon_change
      dmf_new = dmf_copy.GetProcessCopy()
      dmf_new.Info("\nTesting solvent content corresponding to {0} monomers.".format(num_mon))
      dmf_new.out.never_propagate=True
      dmf_new.SetRunDir(dmf_new.rundir+'_'+str(num_mon))
      dmf_new.ph=dmf_new.GetProcess(self.ph.nick)
      dmf_new.dm=dmf_new.GetProcess(self.dm.nick)
      dmf_new.sft=dmf_new.GetProcess('createfree').GetProg('sftools')
      dmf_new.sft.seed=self.sft.seed
      dmf_new.sft.SetRunDir(os.path.join(dmf_new.rundir,'betafromfree'))
      dmf_new.dm.inp.Get('sequence').monomers_asym=num_mon
      dmf_new.inp.Set(dmf_new.dm.inp.Get('sequence'))
      matthews=dmf_new.dm.GetOrAddProcess('matthews')
      try:
        matthews.Run()
      except (matthews.GetProg().ProgramRunError,):
        self.Info("Matthews coeficients for {0} monomers too unlikely.".format(num_mon))
        fom=0.0
      else:
        dmf_new.dm.inp.Get('sequence').solvent_content=matthews.GetProg(supported=True).GetStat('solvent_content')
        dmf_new.RunComb()
        fom=sum(dmf_new.ph.fom_all[5:10])
      if abs(fom-fom_best)<0.02:
        self.Info("Very weak discrimination - the solvent content estimation may be unreliable.")
      if fom>fom_best:
        fom_best=fom
        dmf_best=dmf_new
        if num_mon<=1:
          num_mon_change=0
      elif num_mon_change<0 and dmf_best is self:
        num_mon = num_mon_start
        num_mon_change=1
      else:
        num_mon_change=0
    # evaluate and copy the best results
    if dmf_best is self:
      self.Info("More likely solvent content not found.")
    else:
      self.dm.inp.Set( dmf_best.dm.inp.Get('sequence') )
      self.dm.out.Add( dmf_best.dm.inp.Get('sequence') )
      self.dm.out.Add( dmf_best.out.Get('mapcoef',typ='densmod') )
      self.ph.out.Add( dmf_best.out.Get('mapcoef',typ='combined') )
      self.Info("Solvent content of {0} corresponding to {1} monomers will be used.".format(
        self.dm.inp.Get('sequence').solvent_content, self.dm.inp.Get('sequence').monomers_asym))

  def MergeBiasedPhases(self):
    if not self.inp.Get(is_native=True) or not self.inp.Get(is_native=False) or self.GetParam('dmcyc')<=0 or \
       self.inp.Get(is_native=True).GetResolution(self)>=self.inp.Get(is_native=False).GetResolution(self):
      return
    self.sft.inp.ClearAll()
    self.sft.ClearAnyParams()
    self.sft.runname=self.sft.name+'_mergebiased'
    biased_ph=self.dm.out.Get('mapcoef',typ='densmod',filetype='mtz',col='ph')
    comb_ph=self.ph.out.Get('mapcoef',typ='combined',filetype='mtz',col='ph')
    self.sft.inp.Add(biased_ph)
    self.sft.inp.Add(comb_ph)
    self.sft.MergeMTZ(*self.sft.inp.GetAll(), force_merge=True)
    merged_ph=self.sft.out.AddCopy(biased_ph)
    merged_ph.SetLabel( ['ph','fom','f'], ['PHCOMB_MERGE','FOMCOMB_MERGE', 'FMOD_USED'], ignore_prefix=True )
    merged_ph.SetType('combined')
    self.sft.out.SetFileToChild(merged_ph, self.sft.inp.Get(filetype='mtz').GetFileName('mtz')[:-4]+'_merged.mtz', 'mtz')
    self.sft.SetKey( 'read', ('"'+self.sft.inp.Get(filetype='mtz').GetFileName('mtz')+'"') )
    self.sft.SetKey( 'checkhkl' )
    self.sft.SetKey( 'calc', ('P col', merged_ph.GetLabel('ph')+' = col '+comb_ph.GetLabel('ph')) )
    self.sft.SetKey( 'calc', ('W col', merged_ph.GetLabel('fom')+' = col '+comb_ph.GetLabel('fom')) )
    self.sft.SetKey( 'set', ('label', 'rename', 'col', biased_ph.GetLabel('f')) )
    self.sft.SetKey( merged_ph.GetLabel('f') )
    self.sft.SetKey( 'select', ('col', comb_ph.GetLabel('ph'), 'absent') )
    self.sft.SetKey( 'select', ('col', biased_ph.GetLabel('ph')) )
    self.sft.SetKey( 'calc', ('col', merged_ph.GetLabel('ph')+' = col '+biased_ph.GetLabel('ph')) )
    self.sft.SetKey( 'calc', ('col', merged_ph.GetLabel('fom')+' = 0.1') )
    self.sft.SetKey( 'select', 'all' )
    self.sft.SetKey( 'checkhkl' )
    self.sft.SetKey( 'write', (merged_ph.GetFileName('mtz'), '\nY') )
    self.sft.SetKey( 'quit\nY' )
    self.sft.Run(restore=False, clear_out=False)
    self.ph.out.Add(merged_ph)

  def RunPostprocess(self,restore=True,*args,**kwargs):
    #if self.dmf and self.dmf[0]:
      #self=self.dmf[0]
    if self.GetProcess('createfree'):
      self.processes.remove(self.createfree)
      delattr(self,'createfree')
    ###self.PrintActualLogGraph(finished=True)
    # a workaround for multicomb replacing the residue type (for example IOD always becomes I)
    # this can (and is better to) be removed after this is fixed in multicomb
    if self.GetProcesses(supported=True) and self.ph.GetProg().nick=='multicomb' and \
       self.ph.out.Get('model',typ='substr',filetype='pdb'):
      fixsub=self.ph.AddProcess('fixsubstrpdb',propagate_inp=False)
      fixsub.inp.Add(self.ph.out.Get('model',typ='substr',filetype='pdb'))
      fixsub.AddProg('pdbcur')
      fixsub.GetProg().outfilename['pdb'] = "fixed_"+self.ph.nick+".pdb"
      fixsub.Run()
      self.ph.processes.remove(fixsub)
    # pass larger solv rad to comb
    comb=self.GetCrankParent().prep.GetProcess('comb_phdmmb') if self.GetCrankParent() else None
    if comb and (not comb.GetProcess('dmfull') or not comb.GetProcess('dmfull').GetProcess('dm') or not comb.GetProcess('dmfull').GetProcess('dm').GetProg() or comb.GetProcess('dmfull').GetProcess('dm').GetProg('parrot')):
      if self.dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius'):
        comb.SetParam('solventmask_radius', self.dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius'))
      if self.GetParam('handdet') and self.dmf[0].dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius') and self.dmf[1].dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius'):
        comb.SetParam('solventmask_radius', min(self.dmf[0].dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius'), self.dmf[1].dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius')))
    # take determined solvent cont. to comb
      if self.GetParam('map_segmentation'):
        if self.GetParam('handdet') and self.dmf[0].GetParam('solvent_content') and self.dmf[1].GetParam('solvent_content'):
          comb.SetParam('solvent_content', min(self.dmf[0].GetParam('solvent_content'), self.dmf[1].GetParam('solvent_content')))
        else:
          comb.SetParam('solvent_content',self.GetParam('solvent_content'))
    #self.MergeBiasedPhases()
    self.guess=1
    if self.ph and self.ph.fom_all and self.ph.fom_all[-1]<0.35:
      self.guess=0
    elif self.ph and self.ph.fom_all and self.ph.fom_all[-1]<0.5:
      self.guess=2
    mapseg=self.GetProcess('segmentmap')
    if mapseg:
      self.processes.remove(mapseg)  # useful for comb - the s.cont. will be estimated again after the next building
    process.RunPostprocess(self,restore,*args,**kwargs)

  def PrintActualLogGraph(self, finished=False, beta=False, ncs=False, hand=-1, clear=False):
    if self.queue and self.num_proc>=0:
      queue,parent,logfilehandle,ph_lfh,dm_lfh = self.queue,self.parent_process,self.logfilehandle,self.ph.logfilehandle,self.dm.logfilehandle
      if hasattr(self,'ccp4i2job'):  i2job,self.ccp4i2job=self.ccp4i2job,None
      # the queue and parent process are removed and reattached after putting to the queue
      self.queue,self.parent_process,self.logfilehandle,self.ph.logfilehandle,self.dm.logfilehandle = None,None,None,None,None
      mapseg=self.GetProcess('segmentmap')
      if mapseg:
        self.processes.remove(mapseg)
      queue.put((self.num_proc,self,{'finished':finished,'beta':beta,'ncs':ncs,'clear':clear}))
      if mapseg:
        self.processes.append(mapseg)
      self.queue,self.parent_process,self.logfilehandle,self.ph.logfilehandle,self.dm.logfilehandle = queue,parent,logfilehandle,ph_lfh,dm_lfh
      if hasattr(self,'ccp4i2job'):  self.ccp4i2job=i2job
      return
    if self.opened_loggraph and hand<2:  # loggraph for hand 2 not printed, may be removed completely at some point.
      self.GetLogGraphHandle().seek(0,0)
      self.LGInfo(self.GetCCP4Header())
      self.LGInfo(self.GetLogGraphPrograms())
      if not self.GetVirtPar('no_bias_est'):
        self.LGInfo('\n$TEXT: Beta estimation result: $$ Beta estimation result $$')
        self.LGInfo('Beta bias parameter estimated as {0}\n$$'.format(self.ph.GetParam('beta')))
      if self.GetProcesses(supported=True) and self.ph.fom_all:
        self.LGInfo('\n $TABLE : Fom per cycle:')
        self.LGInfo('$GRAPHS    : FOM per cycle :A:1,2,3:\n$$')
        self.LGInfo('Cycle FOM $$ $$')
        for i,n in enumerate(zip(self.ph.fom_all)):
          self.LGInfo('{0:3}  {1:8.3f}'.format(i+1,n[0]))
        self.LGInfo('$$\n')
        if finished:
          self.LGInfo('\n$TEXT:Result: $$ Final result $$')
          self.LGInfo('The final FOM is {0}\n$$'.format(self.ph.fom_all[-1]))
    if self.rv_report is not None:
      if clear:
        if hasattr(self,'rv_plot'):
          self.rv_plot.parent.parent.Remove()
          del self.rv_plot
        if hasattr(self,'rv_fom'):  del self.rv_fom
        if hasattr(self,'rv_fom2'): del self.rv_fom2
      elif beta:
        handtext = ' hand {}'.format(hand) if hand>0 else ''
        self.rv_report.Text("<BR><i><small>Beta (bias) parameter estimation result{}:</small></i>".format(handtext))
        button=crvapi.init_meta["help_btn_template"].replace('.html','.html#term-what-is-the-parameter-in-density-modification') if "help_btn_template" in crvapi.init_meta else ""
        self.rv_report.Text("&emsp;Beta parameter"+button+" estimated as {0}".format(self.ph.GetParam('beta')), flush=True)
      elif ncs is not False:
        self.rv_report.Text("<BR><i><small>NCS detection result:</small></i>")
        if ncs>=0:
          self.rv_report.Text("&emsp;{0} NCS operators detected.".format(ncs), flush=True)
        else:
          self.rv_report.Text("&emsp;NCS detection disabled.", flush=True)
      elif self.GetProcesses(supported=True) and self.ph.fom_all:
        fom_thr=0.35
        fom_ymax=max([fom_thr]+self.ph.fom_all)+0.01
        if not hasattr(self,'rv_plot'):
          self.rv_plot = self.rv_report.Plot( 'FOM vs cycle', "Cycle", "Figure of merit", legendloc='se', ymin=0.15, ymax=fom_ymax, intx=True )
        if hand<=1 and not hasattr(self,'rv_fom'):
          self.rv_fom = self.rv_plot.PlotLine( ["x"], ["FOM"], color=self.rv_plot.def_color[0] )
          if hand<0:
            self.rv_shadow = self.rv_plot.PlotLine( ["x",(0,)], ["Very weak phases",(fom_thr,)], fill=1 )
          for i,p in enumerate(self.ph.fom_all[:-1]):
            self.rv_fom.Data( i+1, p, x_type=int )
        if hand>1 and not hasattr(self,'rv_fom2'):
          self.rv_fom2 = self.rv_plot.PlotLine( ["x"], ["FOM hand 2"], color=self.rv_plot.def_color[1] )
          self.rv_shadow = self.rv_plot.PlotLine( ["x",(0,)], ["Very weak phases",(fom_thr,)], fill=1 )
          for i,p in enumerate(self.ph.fom_all[:-1]):
            self.rv_fom2.Data( i+1, p, x_type=int )
        if hand>1:
          self.rv_fom2.Data( len(self.ph.fom_all), self.ph.fom_all[-1], x_type=int )
        else:
          self.rv_fom.Data( len(self.ph.fom_all), self.ph.fom_all[-1], x_type=int )
        self.rv_plot.SetProperties(ymax=fom_ymax)
        if hasattr(self,'rv_shadow'):
          self.rv_shadow.Data( len(self.ph.fom_all)+1, fom_thr, x_type=int, flush=True )
        if finished:
          handtext = ' hand {}'.format(hand) if hand>0 else ''
          self.report_fom = self.ph.fom_all[-1]
          self.rv_report.Text("<i><small>Result{}:</small></i>".format(handtext)) # title
          if self.stopped_early:
            self.rv_report.Text('&emsp;'+'Stopping early on FOM higher than the threshold ({0})'.format(self.GetParam('threshold_stop')))
          self.rv_report.Text('&emsp;'+'The final FOM is {0}'.format(self.ph.fom_all[-1]))
          if not self.parent_process or not self.parent_process.ccp4i2:
            self.rv_report.Text("<BR><i><small>Output{}:</small></i>".format(handtext))
            usemtz=crvapi.SetMtzFWT(self,self.out.Get('mapcoef',typ='combined'))
            if self.out.Get('model'):
              self.rv_files = self.rv_report.DataFile(self.out.Get('model').GetFileName(), "xyz", "<small>DM improved map</small>")
              self.rv_files.DataFile(usemtz, "hkl:map")
            else:
              self.rv_files = self.rv_report.DataFile(usemtz, "hkl:map", "<small>DM improved map</small>")
            #self.rv_files.DataFile(self.out.Get('mapcoef',typ='combined',filetype='map').GetFileName('map'), "hkl:ccp4_map")
            #anomdiff=self.out.Get('mapcoef',typ='anom-diff',filetype='map',conv_opts=['no_rewrite'])
            #if anomdiff:
            #  self.rv_files.DataFile(anomdiff.GetFileName('map'), "hkl:ccp4_dmap")
          crvapi.Flush()


  def RunCombHands(self,num_proc=-1,queue=None):
    self.num_proc=num_proc
    self.queue=queue
    if self.logfilehandle:
      self.logfilehandle.close()
    try:
      if num_proc>0:
        self.SetRunDir(self.rundir+'_hand'+str(num_proc+1), reset_subtree=True, change_cwd=True)
        #self.inp=self.inp2
        self.ReplaceInpOut('inp','inp2')
        shutil.copy(self.sft.out.Get('exclude',filetype='mtz').GetFileName(), os.path.join(self.sft.rundir,'workonly_copy.mtz'))
        self.sft.out.SetFileToChild( self.sft.out.Get('exclude',filetype='mtz'), 'workonly_copy.mtz', 'mtz' )
      self.TreatInOutPar()
      self.RunPreprocess()
      self.RunComb()
      # conversions have to be done here rather than in printactuallog since the parent dmfull process is not running, just reporting...
      self.out.Get('mapcoef',typ='combined',filetype='map')
      self.out.Get('mapcoef',typ='anom-diff',filetype='map',conv_opts=['no_rewrite'])
      self.PrintActualLogGraph(finished=True)
    except Exception as e:
      error='Parallel density modification process {0} failed with error: {1}'.format(num_proc, str(e))
      queue.put((-1,error,'err'))
    #else:
      #self.queue=None
      #queue.put((num_proc,self,'end'))

  def TestAdjFeedSolvmask(self,mult):
      # if low fom is estimated then rerun with solv. filter radius adjusted
      # use estimate of fom with bias applied
      # (assumptions: change of fom between 00 and 0 cycle is similar to change betwen 0 and 1; bias hampers fom with square of beta)
      beta=self.ph.GetParam('beta')
      #twice0min1=2*self.ph.fom_all[0]-self.ph.fom_all[1]
      twice0min1=2.25*self.ph.fom_all[0]-1.25*self.ph.fom_all[1]
#      fomest=twice0min1+beta*beta*(self.ph.fom-twice0min1)+max(0.,(self.ph.fom-self.ph.fom_all[-3])*beta*2.5-0.05)
      fomest=twice0min1+beta*beta*(self.ph.fom-twice0min1)+max(0.,0.1*beta-0.05)
      #print('debug',fomest,twice0min1,max(0.,0.1*beta-0.05),max(0.,(self.ph.fom-self.ph.fom_all[-3])*beta*2.-0.05),self.GetParam('incr_radius_threshold_fom'))
      if self.dm.GetProg('parrot') and (not self.dm.GetProg().GetKey('solvent-mask-filter-radius') or mult>1) and \
         fomest<self.GetParam('incr_radius_threshold_fom') and (not self.parent_process or self.parent_process.nick!='comb_phdmmb'):
#         beta*(2*self.ph.fom-twice0min1)/2.<self.GetParam('incr_radius_threshold_fom') and (not self.parent_process or self.parent_process.nick!='comb_phdmmb'):
#        mult=1.5 if twice0min1+beta*(self.ph.fom-twice0min1)<self.GetParam('incr_radius_threshold_fom')-0.025 else 1.25
        mult=min(1.0+20.*(self.GetParam('incr_radius_threshold_fom')-fomest),1.5)
        #print('debug2',mult)
        self.AdjSolvMask(mult)
        #self.ph.out.Get('mapcoef',typ='weighted').custom.append('dmmsk')
        #self.inp.Add(self.ph.out.Get('mapcoef',typ='weighted',custom='dmmsk'))
      return mult

  def RunComb(self):
    # contin is used so that the process can be repeated if needed
    num_failed_bias, feedback, contin, mult, beta_aver = 0, 0, 1, 1, 0.
    while contin:
      #if hasattr(self,'rv_plot') and num_failed_bias:
      if num_failed_bias:
        #self.rv_plot.parent.parent.Remove()
        #del self.rv_plot
        self.PrintActualLogGraph(clear=True)
        beta_aver+=self.ph.GetParam('beta')
      # estimate density modification bias
      extra_corr=False
      if not self.GetParam('no_bias_est') and self.GetParam('biascyc')>0:
        if num_failed_bias>3:
          self.ph.SetParam('beta',beta_aver/num_failed_bias)
          self.Info('Beta estimation seems unstable and will be skipped, with beta set to {}.\n'.format(self.ph.GetParam('beta')))
          self.SetParam('no_bias_est', True)
        else:
          # possibly restarting due to detected NCS estimation issues
          if self.RunWithGivenBeta(bias_est=True):
            continue
          self.EstBias()
          mult = self.TestAdjFeedSolvmask(mult)
#          if mult>1.001 and not feedback and not self.GetParam('map_segmentation') and self.ph.GetProg(supported=True).nick=='refmac' and self.dm.GetProg(supported=True).nick=='parrot':
#            feedback=1
#            continue
          # if beta is low then we should estimate the correction separately - ie more cycles needed...
          if self.ph.GetParam('beta')<0.7 or self.GetParam('optimize_solvent'):
            extra_corr=True
      # run density modification + check bias estimation, apply correction (if bias is estimated)
      #self.ph.out.Get('mapcoef',typ='weighted').custom.append('dmmsk')
      #self.inp.Add(self.ph.out.Get('mapcoef',typ='weighted',custom='dmmsk'))
      if self.RunWithGivenBeta(extra_corr=extra_corr, feedback=feedback):
        # if bias estimation problem then let's reseed and repeat bias estimation
        self.CreateWorkingSetMTZ()
        num_failed_bias+=1
      else:
        contin=0


  def ReInitPhDM(self):
    self.dm.inp.SetCopy(self.inp.fsigf)
    self.dm.inp.SetCopy(self.inp.mapcoef)
    self.dm.inp.SetCopy(self.inp.model)
    # this will have to be disabled for free lunch or other extension schemes!
    # done here as a lot of extra nodata can hinder the dm performance.
    self.dm.GetProg().keep_nodata=False
    self.ph.inp.SetCopy(self.inp.fsigf)
    self.ph.inp.SetCopy(self.inp.model)
    self.ph.inp.SetCopy(self.inp.mapcoef)
    self.stopped_early=False

  def SetPhOutFileName(self,ncyc,bias_est,before_run):
    # adding end_ prefix to the output files to make sure they won't be rewritten when used as an input in a potential following dmfull
    # alternate between two different names needed for parallel comb jobs
    if self.cyc==ncyc and not bias_est and hasattr(self.ph.GetProg(),'outfilename'):
      for i,(key,ofn) in enumerate(self.ph.GetProg().outfilename.items()):
        if before_run:
          if i==0:
            self.phfile_end_sep='_' if self.phfile_end_sep=='.' else '.'
          self.ph.GetProg().outfilename[key] = 'end'+self.phfile_end_sep+ofn
        else:
          self.ph.GetProg().outfilename[key] = ofn[4:]

  def AdjustResol(self,o1,o2):
    # adjusts resolution of mtz from o1 to that of o2
    mtz = gemmi.read_mtz_file(o1.GetFileName('mtz'))
    mtz_resol = gemmi.read_mtz_file(o2.GetFileName('mtz'))
    all_data = np.array(mtz, copy=False)
    mtz.set_data(all_data[mtz.make_d_array() <= mtz_resol.resolution_low()])
    mtz.write_to_file(o1.GetFileName())

  def RunWithGivenBeta(self, bias_est=False, extra_corr=False, no_corr=False, feedback=False, verbose=1):
    """Simple DM+phase combination recycling"""
    self.ph.R_all, self.ph.corr_all, self.ph.fom_all, self.dm.contr_inv_all = [], [], [], []
    beta_saved = 0
    self.ReInitPhDM()
    if bias_est or extra_corr:
      ncyc = self.GetParam('biascyc')
      if extra_corr:
        beta_saved=self.ph.GetParam('beta')
      self.ph.SetParam('beta',1.0)
    elif feedback:
      ncyc = 12 if not self.GetParam('map_segmentation') else 5
      if self.GetParam('map_segmentation'):
        beta_saved=self.ph.GetParam('beta')
        self.ph.SetParam('beta',1.0)
    else:
      ncyc = int(self.GetParam('dmcyc'))
    if verbose and ncyc:
      if bias_est:
        self.Info("Starting {0} bias estimation cycles\n".format(ncyc))
      elif extra_corr:
        self.Info("Starting {0} bias correction cycles\n".format(ncyc))
      elif feedback:
        self.Info("Starting {0} feedback cycles\n".format(ncyc))
      else:
        #refi="and refinement" if mb_mod else ""
        self.Info("Starting {0} density modification cycles\n".format(ncyc))
    if not self.inp.Get('mapcoef',typ=('best','combined')):
      ph_backup=self.ph.BackupAnyPars()
      self.cyc=0
      self.Info(" Cycle 0")
      self.ph.SetParam('comb_with_typ','no_phcomb')
      self.ph.Run()
      self.GetActualStats(verbose, update_loggraph=False)
      self.ph.RestoreAnyPars(ph_backup)
      self.dm.inp.AddCopy(self.ph.out.Get('mapcoef',typ=('best','combined')))
    self.ph.out.Set(self.ph.inp.model, propagate=False)
    init_solmask = self.dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius')
    init_mapseg_solv,prev_solv,max_solv_move=None,None,0.125  # solvent related variables used with map_segmentation
    # the recycling starts here
    for self.cyc in range(1, ncyc+1):
      if self.IsNonFalseParam('solvent_perturb'):# and not bias_est and not extra_corr and not feedback:
#        self.dm.SetParam('solvent_perturb', (self.cyc%3-1.0)*self.GetParam('solvent_perturb'))#*int(ncyc!=self.cyc))
        self.dm.SetParam('solvent_perturb', (self.cyc%3-0.5)*self.GetParam('solvent_perturb') *int(ncyc!=self.cyc and (self.cyc!=2 or self.GetParam('ncs_det') is False)))
#        self.dm.SetParam('solvent_perturb', min(0.0,-0.25+0.015*self.cyc))
      if bias_est or extra_corr or feedback:
        self.dm.SetParam('solvent_perturb',0.0)
      inp_solmask = self.dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius')
      if feedback:
#        self.dm.GetProg(supported=True).SetKey('solvent-mask-filter-radius',min(max(0.75,self.dm.GetProg().GetStat('radius_auto')-2.0)+self.cyc*0.5,init_solmask-0.5),keep_previous=False)
        self.dm.GetProg(supported=True).SetKey('solvent-mask-filter-radius',min(max(0.75,init_solmask-3.0)+self.cyc*0.5,init_solmask),keep_previous=False)
#        self.dm.GetProg(supported=True).SetKey('solvent-mask-filter-radius',min(1.0+self.cyc*0.5,init_solmask),keep_previous=False)
      if verbose:
        self.Info(" Cycle {0}".format(self.cyc))
      self.out.ClearAll(propagate=False)
      self.SetNCSParams(feedback)  # parrot only (skip ncs det. in feedback)
      self.ref_use=[]
      fom=self.ph.fom_all
      if self.GetParam('map_segmentation') and ((self.cyc-3)%5!=0 or (self.cyc>3 and self.cyc>=ncyc-1) or (self.cyc>6 and fom and fom[-1]>self.GetParam('threshold_stop')-0.03)):# or self.cyc==6):# and self.cyc!=4):
        inthemiddle,mapseg=True,self.GetProcess('segmentmap')
        if not mapseg:
          inthemiddle=False
          mapseg=self.AddProcess('segmentmap')
          matthews=self.GetOrAddProcess('matthews')
          matthews.out.never_propagate=True
          matthews.Run()
          matthews_probs = [float(prob) for prob in matthews.GetProg().GetStat('probab_matth_all') if prob!='']
          matthews_solvs = [float(sol) for sol in matthews.GetProg().GetStat('solvent_content_all') if sol]
        if self.dm.inp.Get('mapcoef',typ='weighted'):
          mapseg.inp.Set(self.dm.inp.Get('mapcoef',typ='weighted'))
        else:
          mapseg.inp.Set(self.dm.inp.Get('mapcoef',typ=('combined','best')))
        mapseg.Run()
        mapseg.out.Get('mapcoef',typ='mask').custom.append('dmmsk')
        if not self.IsInputtedParam('solvent_content') and self.parent_process.nick!='comb_phdmmb':# and self.parent_process.nick!='comb_phdmmb':# and self.cyc==1 and not inthemiddle:
          prev_solv=self.GetParam('solvent_content') if self.GetParam('solvent_content') else self.inp.Get(has_solvent_content=True).solvent_content
          if self.cyc==1 and not inthemiddle:
            logit_estsolv = mapseg.GetLogitFromSolv(mapseg.est_solvcont)
            matthews_logits = [mapseg.GetLogitFromSolv(sol) for sol in matthews_solvs]
            matthews_probs2 = [ mapseg.GetProbFromLogitDiff(abs(logit_diff-logit_estsolv) if logit_diff is not None else None) for logit_diff in matthews_logits ]
            matthews_probs_joint = [p2*matthews_probs[i] for i,p2 in enumerate(matthews_probs2)]
            init_solv = matthews_solvs [ matthews_probs_joint.index( max(matthews_probs_joint) ) ]
            if max(matthews_probs_joint)<0.25 or abs(mapseg.est_solvcont-init_solv)>0.1:
              if max(matthews_probs_joint)<0.005:
                init_solv = mapseg.est_solvcont
              else:
                matthews_probs = [0]+matthews_probs+[0]
                matthews_solvs = [1.]+matthews_solvs+[0.]
                init_solv=mapseg.GetExpSolvFromCombDist(matthews_solvs,matthews_probs)
            init_solv=round(init_solv,3)
            matthews_solv=self.inp.Get(has_solvent_content=True).solvent_content
            #if matthews_solv and matthews_solv>mapseg.est_solvcont_min and matthews_solv<mapseg.est_solvcont_max and abs(matthews_solv-init_solv)<0.1: 
             # if mapseg.est_solvcont_max-mapseg.est_solvcont_min>0.2:
            #    init_solv=matthews_solv
             # elif mapseg.est_solvcont_max-mapseg.est_solvcont_min>0.12:
            #    init_solv=(mapseg.est_solvcont+matthews_solv)/2.0
          #self.SetParam('solvent_content',init_solv+max(-max_solv_move,min(max_solv_move,(mapseg.est_solvcont+prev_solv)/2.-init_solv))) # do not go further than max_solv_move from the initial estimate
            self.SetParam('solvent_content',init_solv)
            if abs(matthews_solv-self.GetParam('solvent_content'))>0.001:
              self.Info("Using solvent content {} estimated using map segmentation instead of the Matthews estimate {}".format(self.GetParam('solvent_content'),matthews_solv))
            if not self.IsInputtedParam('dmcyc'):
              init_dmcyc=self.GetParam('dmcyc')
              self.SetDMcycSolv(init_solv)
              if init_dmcyc!=self.GetParam('dmcyc'):
                self.Info('Number of DM cycles to run adjusted to {}'.format(self.GetParam('dmcyc')))
        self.dm.inp.Add(mapseg.out.Get('mapcoef',typ='mask',custom='dmmsk'))
      self.dm.GetProg().runname=self.dm.GetProg().name+'_cyc'+str(self.cyc)
      self.dm.Run()
      if self.GetParam('map_segmentation'):
        self.AdjustResol(self.dm.out.Get('mapcoef'),self.dm.inp.Get('fsigf'))
      if self.dm.GetProg('shelxe'):
        self.dm.out.Get('mapcoef',typ='densmod',filetype='mtz',inp_cont=self.inp.Get('fsigf',typ='average',col='f'))
      if feedback:
        if inp_solmask is None and self.dm.GetProg(supported=True).GetKey('solvent-mask-filter-radius'):
          self.dm.GetProg(supported=True).UnsetParam('solvent-mask-filter-radius',is_key=True)
        elif inp_solmask is not None:
          self.dm.GetProg(supported=True).SetKey('solvent-mask-filter-radius',inp_solmask,keep_previous=False)
      self.ph.inp.AddCopy(self.dm.out.Get('mapcoef',typ='densmod'))
      if self.ph.out.model: # take the latest model if available (not available for MLHL)
        self.ph.inp.Set(self.ph.out.model)
      self.ph.GetProg().runname=self.ph.GetProg().name+'_cyc'+str(self.cyc)
      self.SetPhOutFileName(ncyc,bias_est,before_run=True)
      self.ph.Run()
      self.SetPhOutFileName(ncyc,bias_est,before_run=False)
      if self.ph.GetProg().nick=='refmac': 
        self.ph.inp.Clear('datafile')
        if self.cyc>1 and self.cyc<ncyc:
          self.ph.inp.Set(self.ph.out.Get('datafile',typ='dluz'))
        # comb_phdmmb specific - mainly to optimize the initial cycle for good starting maps
        if self.parent_process and self.parent_process.nick=='comb_phdmmb' and self.cyc==1 and bias_est:
          if self.ph.R<self.parent_process.R_tresh+0.06 and self.GetParam('dmcyc')>=self.parent_process.dm_cyc_init:
            self.SetParam('dmcyc', self.parent_process.dm_cyc_init-1  if self.parent_process.dm_cyc_init>1  else 1)
            self.ph.SetParam('cycles', 3)
          if not self.parent_process.IsInputtedParam('no_bias_est') and self.ph.R<self.parent_process.R_tresh-0.05:
            self.SetParam('no_bias_est',1)
            self.Info('Skipping bias estimation as the model looks good. Please restart with the no_bias_est parameter set to False to supress this choice.')
            return 1
      if not self.dm.GetProg('parrot'):
        self.dm.inp.Clear('mapcoef')
      if self.dm.GetProg('shelxe'):
        self.dm.inp.AddCopy(self.ph.out.Get('mapcoef',typ='combined'))
      else:
        self.dm.inp.AddCopy(self.ph.out.Get('mapcoef',typ='weighted'))
      # models can be used for ncs detection parrot (otherwise ignored)
      ncs_mod = self.inp.GetAll('model',custom='ncs') + self.out.GetAll('model',typ='substr')
      if ncs_mod:
        self.dm.inp.Set( ncs_mod )
      if bias_est and self.cyc>0:
        self.dm.GetProg().keep_nodata=True
        self.MakeFreeDMInpAbsent()
      self.GetActualStats(verbose, update_loggraph=not bias_est and not extra_corr and not feedback)
      if self.GetNCSParams(update_loggraph=not bias_est and not extra_corr and not feedback):  # parrot only
        self.ReInitPhDM()
        return 1
      # automatically increase the filtering radius for weak phases (only for parrot as of now)
      #if self.dm.GetProg('parrot') and not self.dm.GetProg().GetKey('solvent-mask-filter-radius') \
        # and self.cyc==5 and self.ph.fom<self.GetParam('incr_radius_threshold_fom'):
        #self.AdjSolvMask(1.25)
        #self.ReInitPhDM()
        #return 1
      if not bias_est and not self.GetParam('no_bias_est') and not no_corr and \
         self.cyc==self.GetParam('biascyc') and self.EstBias(correction=True,extra_correction=beta_saved):
        self.ReInitPhDM()
        return 1
      # stop early
      if not bias_est and self.ph.fom>self.GetParam('threshold_stop') and self.cyc>=8:
        self.stopped_early=True
        break
    if extra_corr:
      if self.TestAdjFeedSolvmask(1)>1.001 and not feedback:
        feedback=1
      if self.RunWithGivenBeta(no_corr=True,feedback=feedback):
        self.ReInitPhDM()
        return 1
    elif feedback:
      if not self.GetParam('map_segmentation'): ###!!!
        self.ph.out.Get('mapcoef',typ='weighted').custom.append('dmmsk')
        self.inp.Set(self.ph.out.Get('mapcoef',typ='weighted',custom='dmmsk'))
      #feedback+=1 #if we wanted more feedback recycling
      feedback=0 #if not self.GetParam('map_segmentation') or feedback>2 else feedback
      if self.GetParam('map_segmentation'): # return to the determined beta
        self.ph.SetParam('beta',beta_saved)
      if self.RunWithGivenBeta(no_corr=True,feedback=feedback):
        self.ReInitPhDM()
        return 1
    self.ReInitPhDM()
    return 0

  def AdjSolvMask(self,mult,maxim=8.):
    self.dm.GetProg().SetKey( 'solvent-mask-filter-radius', min(self.dm.GetProg().GetStat('radius_auto')*mult,maxim), keep_previous=False )

  def GetActualStats(self,verbose,update_loggraph):
    self.ph.fom = self.ph.GetProg().GetStat('fom')[-1]
    self.ph.fom_all.append(self.ph.fom)
    if verbose:
      self.Info("  Overall MEAN FOM is {0}".format(self.ph.fom))
    if self.ph.nick=='ref':
      self.ph.R = self.ph.GetProg().GetStat('rfact')[-1]
      self.ph.R_all.append(self.ph.R)
      if verbose:
        self.Info("  Overall Rcomb-factor is {0}".format(self.ph.R))
    else:
      self.ph.corr = self.ph.GetProg().GetStat('correl')
      self.ph.corr_all.append(self.ph.corr)
      if verbose and (self.ph.nick=='ref' or self.ph.GetProg().nick!='refmac'):
        self.Info("  Overall Correlation is {0}".format(self.ph.corr))
      if self.dm.GetProg('solomon'):
        self.dm.contr_inv=self.dm.GetProg('solomon').GetStat('contrast_inv')
        self.dm.contr_inv_all.append(self.dm.contr_inv)
        if verbose:
          self.Info("  Inverse contrast is {0}".format(self.dm.contr_inv))
    if verbose and update_loggraph:
      self.PrintActualLogGraph()
      if self.parent_process and self.parent_process.nick=='comb_phdmmb' and self.ph.nick=='ref':
        self.parent_process.UpdateResults('ph', skip_output=self.parent_process.GetParam('num_parallel')>1)

  def SetNCSParams(self,skip_ncs=False):
    if not self.dm.GetProg('parrot') or self.GetParam('ncs_det') is False:
      self.dm.SetParam('ncs_det', False)
      return
    # we do not want to run ncs determination every cycle
    if ((self.cyc-2)%10==0 or (self.cyc-4)%13==0) and not skip_ncs:
      self.ncs_cyc = self.cyc
      self.dm.SetParam('ncs_det', self.GetParam('ncs_det'))
      self.dm.GetProg('parrot').SetKey('ncs-operator',False,keep_previous=False)
    else:
      self.dm.SetParam('ncs_det', False)
    if self.cyc==1:
      self.dm.GetProg('parrot').SetKey('ncs-operator',False,keep_previous=False)

  def GetNCSParams(self, update_loggraph=False):
    #if not self.dm.GetProg('parrot') or self.GetParam('ncs_det') is False:
    if not self.dm.GetProg('parrot'):
      return 0
    #if self.dm.GetParam('ncs_det') is not False:
    ncs_opers=self.dm.GetProg('parrot').GetStat('ncs_operator')
    ncs_oper_correls=[float(c) for c in self.dm.GetProg('parrot').GetStat('ncs_operator_correl') if c!=''] if ncs_opers else []
    self.dm.GetProg('parrot').SetKey('ncs-operator',False,keep_previous=False)
    if self.cyc==self.ncs_cyc:
      for ncs in ncs_opers:
        self.dm.GetProg('parrot').SetKey('ncs-operator',ncs)
    elif ncs_opers:
      for i,(ncs,ncsp,c,cp) in enumerate(zip(ncs_opers,self.ncs_opers,ncs_oper_correls,self.ncs_oper_correls)):
        if c>cp:
          self.dm.GetProg('parrot').SetKey('ncs-operator',ncs)
        else:
          self.dm.GetProg('parrot').SetKey('ncs-operator',ncsp)
          ncs_opers[i]=ncsp
          ncs_oper_correls[i]=cp
    self.ncs_opers,self.ncs_oper_correls=ncs_opers,ncs_oper_correls
    if self.dm.GetParam('ncs_det') is not False and (self.dm.GetParam('ncs_det_ha') or self.dm.GetParam('ncs_det_mr')):
      self.Info('  {0} NCS operators detected by {1}.'.format(len(ncs_opers), self.dm.GetProg('parrot').name))
      if update_loggraph and self.cyc>3 and self.cyc<10:
        self.PrintActualLogGraph(ncs=len(ncs_opers))
        if self.parent_process and self.parent_process.nick=='comb_phdmmb' and self.parent_process.ncs_oper!=len(ncs_opers):
          self.parent_process.ncs_oper=len(ncs_opers)
          self.parent_process.PrintActualLogGraph(update='ncs')
    # if NCS detection seems to go wrong disable it
    # also this assumes the phases should improve (ie starting phases are from phasing)!
    if self.GetParam('ncs_det') is None and self.dm.GetParam('ncs_det') is not False and \
       self.ph.nick=='phcomb' and not self.GetParam('optimize_solvent') and \
       self.cyc==self.ncs_cyc+1 and self.ph.fom_all[self.cyc-1]+0.003<self.ph.fom_all[self.cyc-3]:
      self.Info("Disabling automatic NCS detection and restarting DM.\n")
      self.SetParam('ncs_det', False)
      self.dm.GetProg('parrot').SetKey('ncs-operator',False,keep_previous=False)
      self.PrintActualLogGraph(ncs=-1)
      return 1
    return 0

  def EstBias(self, correction=False, extra_correction=False): 
    """Estimate DM bias from performed bias estimation jobs"""
    # run sftools to calculate the correlations
    self.sft.ClearAnyParams()
    self.sft.runname=self.sft.name+'_correl'
    if not correction and not extra_correction:
      self.sft.runname+='_betaest'
    if extra_correction:
      self.ph.SetParam('beta',extra_correction)
    #self.sft.MergeMTZ(*self.sft.inp.GetAll())
    exclude_sft=self.sft.out.Get('exclude',col='free')
    mapc_dm_out=self.dm.out.Get('mapcoef',typ='densmod',filetype='mtz',col='f')
    mapc_obs_in=self.inp.Get('fsigf',typ='average',filetype='mtz',col='f')
    self.sft.SetKey( 'read', ('"'+exclude_sft.GetFile('mtz').name+'"', 'col', "\"{0}\"".format(exclude_sft.GetLabel('free'))) )
    self.sft.SetKey( 'read', ('"'+mapc_dm_out.GetFile('mtz').name+'"', 'col', "\"{0}\"".format(mapc_dm_out.GetLabel('f'))) )
    self.sft.SetKey( 'read', ('"'+mapc_obs_in.GetFile('mtz').name+'"', 'col', "\"{0}\"".format(mapc_obs_in.GetLabel('f'))) )
    self.sft.SetKey( 'select', ('col', exclude_sft.GetLabel('free'), '=1') )
    self.sft.SetKey( 'checkhkl' )
    self.sft.SetKey( 'correl', ('col', mapc_dm_out.GetLabel('f')+' '+mapc_obs_in.GetLabel('f')) )
    self.sft.SetKey( 'select', 'invert' )
    self.sft.SetKey( 'checkhkl' )
    self.sft.SetKey( 'correl', ('col', mapc_dm_out.GetLabel('f')+' '+mapc_obs_in.GetLabel('f')) )
    self.sft.SetKey( 'quit\nY' )
    self.sft.Run(restore=False, clear_out=False)
    # get the correlations from the log
    corr=self.sft.GetStat('correl',multiple=True)
    if len(corr)!=2:
      common.Error( "Unexpected number of matches in {0}{1}: {2}".format(
        self.sft.runname,self.sft.log_suffix,len(corr)) )
    corr=(float(corr[0]), float(corr[1]))
    # now calculate the estimation or correction to the estimation
    thres=0.15
    if not correction:
      self.ph.SetParam('beta',round(corr[1]/corr[0],3))
    else:
      perc=corr[1]/corr[0]-1.0
      if (perc>thres or perc<-thres) and (corr[0]>2*thres or perc>2*thres or perc<-2*thres):
        self.Info('Too unrepresentative free set. Will need to generate new free set and start from beginning again.\n')
        return 1
      else:
        self.ph.SetParam('beta', round(self.ph.GetParam('beta')*corr[0]/corr[1],3))
    if self.ph.GetParam('beta')<0.2:
      self.ph.SetParam('beta',0.2)
    elif self.ph.GetParam('beta')>1.0:
      self.ph.SetParam('beta',1.0)
    if correction:
      self.Info("Beta corrected to {0}\n".format(self.ph.GetParam('beta')))
    else:
      self.Info("Beta estimated as {0}\n".format(self.ph.GetParam('beta')))
      self.PrintActualLogGraph(beta=True)
    return 0


  def CreateWorkingSetMTZ(self,new_seed=False):
    """Create MTZ containing only working set (for bias estimation) reflections"""
    self.createfree=self.GetOrAddProcess('createfree')
    self.createfree.SetVirtPar('freetype','freebias')
    self.createfree.inp.Set(self.inp.Get('fsigf',typ='average',filetype='mtz'))
    self.createfree.SetRunDir('betafromfree')
    self.sft=self.createfree.GetOrAddProg('sftools')
    self.sft.SetRunDir()
    # this makes sure that nodata are purged
    self.sft.MergeMTZ( self.createfree.inp.Get(filetype='mtz'), force_merge=True )
    self.sft.ClearAnyParams()
    self.sft.runname=self.sft.name+'_createfree'
    if new_seed:
      self.sft.seed=0
    else:
      self.sft.seed+=1
    self.sft.SetKey( 'calc', ('seed', self.sft.seed), keep_previous=False )
    self.createfree.Run()

  def MakeFreeDMInpAbsent(self):
    exclude_sft=self.sft.out.Get('exclude',col='free')
    self.sft.ClearAnyParams()
    self.sft.inp.ClearAll()
    self.sft.runname=self.sft.name+'_workonly'
    # catch all objects needed for input of the DM program
    self.dm.GetProg().BackupAnyPars()
    for o in self.dm.GetProg().CatchMTZObjects():
      self.sft.inp.Add(o)
    # merge and treat parameters and input
    sft_rundir=self.sft.rundir
    self.sft.SetRunDir(os.path.join(sft_rundir,'workingonly'),change_cwd=True)
    self.sft.MergeMTZ( *self.sft.inp.GetAll(), force_merge=True )
    for obj in self.sft.inp.GetAll(filetype='mtz'):
      excl=self.sft.out.AddCopy(obj,propagate=False)
      excl.AddCustomTag('workonly')
      self.sft.out.SetFileToChild( excl, 'workonly.mtz', 'mtz' )
    self.sft.SetKey( 'read', '"'+self.sft.inp.Get(filetype='mtz').GetFileName('mtz')+'"' )
    self.sft.SetKey( 'read', ('"'+exclude_sft.GetFileName('mtz')+'"', 'col', "\"{0}\"".format(exclude_sft.GetLabel('free'))) )
    for obj in self.sft.inp.GetAll(filetype='mtz'):
      for c,l in obj.GetAllLabels():
        self.sft.SetKey( 'absent', ('col', l, 'if col', exclude_sft.GetLabel('free'), '= 0') )
    self.sft.SetKey( 'write', (excl.GetFileName('mtz'), '\nY') )
    self.sft.SetKey( 'quit\nY' )
    self.sft.Run(restore=False, clear_out=False)
    self.sft.SetRunDir(os.path.join(sft_rundir))
    self.dm.inp.Set(self.sft.out.fsigf)
    self.dm.inp.Set(self.sft.out.mapcoef)
