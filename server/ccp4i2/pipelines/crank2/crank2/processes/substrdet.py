#!/usr/bin/python
import heapq
import math
import os
import shutil
import sys

from .. import common
from .. import data
from ..process import crvapi, process
from ..program import program

par=common.parameter

class substrdet(process):
  name="substructure determination"
  short_name=name
  supported_progs=["crunch2","shelxd","prasa"]
  supported_params = {}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True )
  supported_params['threshold_stop'] = par( desc='Stop when this threshold is reached (CFOM for shelx, FOM for crunch)', typ=(float,bool) )
  supported_params['catch_output'] = par( desc='Catch program output continously (always True if threshold_stop)', typ=bool )
  supported_params['num_atoms'] = par( desc='Number of atoms (peaks) to search for', typ=int )
  supported_params['num_trials'] = par( desc='(Maximal) Number of substructure detection trials', typ=(int,bool) )
  supported_params['patterson_seed'] = par( desc='Seed from Patterson (shelxd/crunch2 only)', typ=bool )
  supported_params['num_threads'] = par( desc='Number of CPU threads (shelxd/prasa only)', typ=(int,bool) )
  supported_params['high_res_cutoff'] = par( desc='High resolution cutoff (in A)', typ=(float,bool), share=False )
  supported_params['high_res_cutoff_radius'] = par( desc='Radius of high resolution cutoff search (in A)', typ=(float,bool) )
  supported_params['high_res_cutoff_step'] = par( desc='Step of high resolution cutoff search (in A)', typ=(float,bool) )
  supported_params['high_res_cutoff_cchalf'] = par( desc='Determine high resolution cutoff from CCanom1/2 (if avail.)', typ=bool )
  supported_params['min_dist_atoms'] = par( desc='Minimum distance between atoms (shelx only)', typ=(float,bool) )
  supported_params['min_dist_symm_atoms'] = par( desc='Minimum distance between symmetry related atoms (shelx only)', typ=(float,bool) )
  supported_params['num_dsul'] = par( desc='Number of disulfides to be searched for as 1 atom (shelx only)', typ=(float,bool) )
  supported_params['threshold_weak'] = par( desc='Phases are considered weak if below this threshold (CFOM for shelx, FOM for crunch)', typ=(float,bool) )
  supported_params['optimize_sol'] = par( desc='Optimize the found solutions (prasa only; 0/False=disabled, 1=only the best, 2=all candidates;  default: 2 if below threshold_weak, otherwise 1)', typ=(int,bool) )
  supported_params['afroprasa'] = par( desc='Use default parameters for the Afro-Prasa workflow (prasa only)', typ=bool )

  def Init(self):
    self.min_score=1000.
    self.score=-1000.
    self.score_cut=-1000.
    self.score_adj=-1000.
    self.fom_unopt=-1000.
    self.prog=None
    self.refsol=False
    self.trial=0
    self.trial_adj=0
    self.trial_temp=0
    self.stop_message=""
    self.stop_sent=False
    self.stats_table=[]
    self.check_solution=(-1,-1)
    self.resel_cycle=2000
    self.cutoff2present=None
    self.solcheck_skipped=False
    #self.resel_cycle=20
    self.resel_init_half=self.resel_cycle/2
    self.hist_bin_size, self.hist_num_bins, self.hist_min = 1.0, 200, -10.0
    self.cfom_hist = [ [(float(i)+self.hist_min-0.5)*self.hist_bin_size, 0]  for i in range(self.hist_num_bins) ]
    self.thr_up, self.thr_down = 0., 0.
    self.fom_unopt_all = [-1000.,]
    self.guess=None # a crude classification guess: 0=probably not,1=who knows,2=probably yes
    self.stop_file='stop_file'

  def TreatInOutPar(self, set_all_par=False, skip_cutoffest=False):
    cp=self.GetCrankParent()
    if not self.GetProg(supported=True):
      if program.from_name('shelxd',None).CheckBinary(silent=True) and (not cp or not cp.jscofe or self.GetParam('target')!='SAD' or self.GetParam('afroprasa') is False):
        self.AddProg('shelxd')
      else:
        self.AddProg('prasa')
        if cp and cp.GetProcess('faest') and cp.GetProcess('faest').GetProg(supported=True).nick=='afro' and not self.IsParam('afroprasa'):
          self.SetParam('afroprasa',True)
    prog=self.GetProg(supported=True)
    at = prog.GetKey('sfac')
    if not at and self.inp.Get(has_atomtypes=True):
      at = self.inp.Get(has_atomtypes=True).GetAtomType()
    # catch output continously (eg to stop early if certain stat is reached) - shelxd only as of now.
    if self.GetParam('catch_output') is False:
      self.Info('Early stop on score (cc/fom) threshold disabled.')
    else:
      prog.interact_output=True
      if self.IsTrueOrNoneParam('threshold_stop'):
        if prog.nick=='shelxd':
          self.SetParam('threshold_stop',75.0)
          if self.GetParam('target')=='MAD':
            self.SetParam('threshold_stop',120.0)
        elif prog.nick=='crunch2':
          self.SetParam('threshold_stop',1.0)
        elif prog.nick=='prasa':
          if at in ('S','SE'):
            self.SetParam('threshold_stop',26.0)
          else:
            self.SetParam('threshold_stop',31.0)
      if self.IsTrueOrNoneParam('threshold_weak'):
        if prog.nick=='shelxd':
          self.SetParam('threshold_weak',45.0)
          if self.GetParam('target')=='MAD':
            self.SetParam('threshold_weak',65.0)
        elif prog.nick=='prasa':
          self.SetParam('threshold_weak',17.5)
    # we are changing the shelx/crunch2 program defaults here!
    # the original defaults (if any) will be used when setting threshold_stop, num_trials, min_dist_atoms to False.
    if self.IsTrueOrNoneParam('num_trials'):
      if prog.nick=='shelxd':
        if self.GetParam('target')=='MAD':
          self.SetParam('num_trials',500)
        else:
          self.SetParam('num_trials',2000)
      elif prog.nick=='crunch2':
        self.SetParam('num_trials',20)
      elif prog.nick=='prasa':
        if not cp or not cp.jscofe:
          self.SetParam('num_trials',2000)
        else:
          self.SetParam('num_trials',1000)
    # use cutoff from cc1/2 if unmerged SAD data inputted and cc1/2 known from shelxc
    # (cc1/2 can be later also calculated by other means)
    if self.GetParam('high_res_cutoff_cchalf') is None and not self.IsInputtedParam('high_res_cutoff') and \
       cp.GetParam('target')!='MAD' and self.inp.Get('fsigf',typ='plus',filetype=('HKL','sca'),try_convert=False) and \
       (cp and self in cp.processes and cp.processes.index(self)-1>=0 and \
       cp.processes[cp.processes.index(self)-1].nick=='faest' and \
       'shelxc' in [p.nick for p in cp.processes[cp.processes.index(self)-1].programs]):
      self.SetParam('high_res_cutoff_cchalf',True)
    if self.IsTrueOrNoneParam('high_res_cutoff') and (self.GetParam('target') in ('SAD','SIRAS') or prog.nick=='prasa'):
      # automatic cutoff unless shelx pipeline is used
      if self.GetParam('high_res_cutoff') is None and prog.nick in ('shelxd','prasa') and cp and not cp.GetProcess('phdmmb'):
        self.SetParam('high_res_cutoff',True)
      if prog.nick=='crunch2' or not self.GetParam('high_res_cutoff'):
        self.SetParam('high_res_cutoff', self.EstimateCutOff(halfang=True))
      elif not skip_cutoffest:
        self.SetParam('high_res_cutoff', self.EstimateCutOff())
    if prog.nick=='prasa' and self.IsParam('high_res_cutoff'):
      if self.IsTrueOrNoneParam('high_res_cutoff_radius'):
        self.SetParam('high_res_cutoff_radius',0.5)
      if self.IsTrueOrNoneParam('high_res_cutoff_step'):
        if self.GetParam('high_res_cutoff_radius'):
          self.SetParam('high_res_cutoff_step', self.GetParam('high_res_cutoff_radius')/2.);
        else:
          self.SetParam('high_res_cutoff_step',0.25)
    # changing SHELXD DSUL defaults!
    seq_obj=self.inp.Get('sequence', filetype=data.sequence.supported_filetypes, try_convert=False, has_monomers_asym=True)
    if prog.nick=='shelxd' and seq_obj and not prog.IsKey('DSUL') and not self.IsParam('num_dsul') and at=='S':
      anom_res = prog.GetKey('SHEL')[1] if len(prog.GetKey('SHEL',as_list=True))>1 else self.GetParam('high_res_cutoff')
      # this may be bit off for MAD... reading SHEL from shelxc's ins would be the best but might be slow for defaults generation
      if not anom_res or anom_res is True:
        for o in self.inp.GetAll(typ='plus'):
          if o.GetCrystalName()!='native' and not 'native' in o.custom and self.inp.Get(typ='average',xname=o.GetCrystalName()):
            anom_res = self.inp.Get(typ='average',xname=o.GetCrystalName()).GetResolution(self,accept_none=True)
            if anom_res:
              anom_res = anom_res+0.5
              break
      if anom_res and anom_res>=2.1:
        if not seq_obj.seqstr:
          seq_obj.GetSequenceString()
        self.SetParam('num_dsul', int(seq_obj.seqstr.count('C')*seq_obj.monomers_asym/2))
      else:
        self.SetParam('num_dsul', 0)
    # changing SHELXD MIND defaults!
    if self.IsTrueOrNoneParam('min_dist_atoms') and prog.nick=='shelxd':
      # this should also include resolution check (only to do it when better than 2A or so)
      if at=='S' and not prog.GetKey('DSUL') and not self.GetParam('num_dsul'):
        self.SetParam('min_dist_atoms',-1.5)
      else:
        self.SetParam('min_dist_atoms',-3.0)
    if self.IsTrueOrNoneParam('min_dist_symm_atoms') and prog.nick in ('shelxd','prasa'):
      if at in ('SE','S'):
        self.SetParam('min_dist_symm_atoms',3)
      else:
        self.SetParam('min_dist_symm_atoms',-0.1)
    # changing SHELXD NTPR and PATS defaults for larger structures!  on George's request.
    if prog.nick=='shelxd':
      num_at = prog.GetKey('FIND') if prog.IsKey('FIND') else self.GetParam('num_atoms')
      if not num_at and self.inp.Get(typ='substr',has_num_atoms=True) and self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms:
        num_at = self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms
      if num_at and num_at >= 25:
        # there seems to be no way to turn patterson seeding off other than manual deletion from ins... 
        # which is done at the preprocessing and thus the keyword has to be set here rather than in the shelxd wrapper...
        if not prog.IsKey('PATS') and not self.IsParam('patterson_seed'):
          self.SetParam('patterson_seed', False)
      # George advises 8*FIND for NTPR.  Implemented so that the current minimum of 100 remains for smaller substr.
      if num_at and num_at >= 13:
        if not prog.IsKey('NTPR'):
          prog.SetKey('NTPR', num_at*8)
    if self.GetParam('patterson_seed') is None and prog.nick=='crunch2':
      self.SetParam('patterson_seed')
    ### use Afro with chargeflip 2 and tail histogram matching.
    if self.GetParam('afroprasa'):
      if not prog.IsKey('chargeflip'):
        prog.SetKey('chargeflip',2)
      if not prog.IsKey('histmatch'):
        prog.SetKey('histmatch',2)
    process.TreatInOutPar(self,set_all_par)

  def EstimateCutOff(self,halfang=False,use_faest_08=False,verbose=True):
    cut08=1000.
    # do not try to determine if called only to check binaries...
    if self.GetCrankParent() and hasattr(self.GetCrankParent(),'check_binaries'):  return None
    cutoff=None
    if halfang:
      mtz=[o  for o in self.inp.GetAll(filetype='mtz',typ=('fa','plus'))  if o.GetCrystalName()!='native' and not 'native' in o.custom]
      if mtz:
        cutoff = mtz[0].GetResolution(self)+0.5
        return cutoff
      else:
        common.Warning('Resolution for cutoff could not be determined.')
    # only do this if ecalc was used in the previous step - it does not work for E's from shelxc...
    cp, fa_obj = self.GetCrankParent(), None
    if cp and self in cp.processes and cp.processes.index(self)-1>=0:
      faest = next(p for p in cp.processes if p.nick=='faest')
      if faest and faest.nick!='shelxc':
        fa_obj = faest.inp.Get('fsigf',typ=('delta-anom'),col='e')
        if not fa_obj:
          fa_obj = self.inp.Get('fsigf',typ=('delta-anom'),col='e')
      if faest and hasattr(faest,'rescut08') and use_faest_08:
        cut08=faest.rescut08
    # this should always happens in emulation mode (to get the default cutoff from GUI) but fa's/delta's needed...
    # it may currently happen that the cutoff in emulation mode is different from the real run if different fa's are produced!
    if not fa_obj:
      ecalc=self.AddProg('ecalc',propagate_out=False)
      if not self.inp.Get('fsigf',typ=('delta-anom'),filetype='mtz'):
        mtzmadmod=self.AddProg('mtzMADmod',propagate_out=False)
        try:
          mtzmadmod.Run()
        except:
          if hasattr(sys,'exc_clear'): sys.exc_clear()
        else:
          ecalc.inp.Set(mtzmadmod.out.Get(typ='delta-anom'))
      # currently we prefer delta's over fa's as they are available at the defaults generation, ie to 
      # make sure the defaults match the used params
      if self.inp.Get('fsigf',typ=('delta-anom'),filetype='mtz'):
        ecalc.inp.Set(self.inp.Get('fsigf',typ=('delta-anom'),filetype='mtz'))
      try:
        ecalc.Run()
      except:
        if hasattr(sys,'exc_clear'): sys.exc_clear()
      else:
        fa_obj=ecalc.out.Get('fsigf',typ=('fa','delta-anom'),col='e')
        #self.inp.Add(fa_obj)
      self.programs.remove(ecalc)
    if fa_obj:
        sft=self.AddProg('sftools')
        sft.runname=sft.name+'_compl'
        sft.SetKey('read', '"'+fa_obj.GetFileName('mtz')+'"')
        sft.SetKey('reduce')
        sft.SetKey('select', ('col', fa_obj.GetLabel('f'), 'present'))
        sft.SetKey('select', ('col', fa_obj.GetLabel('sigf'), '>', '0'))
        sft.SetKey('purge\nY')
        sft.SetKey('checkhkl')
        sft.SetKey('complete', ('table', 'shells', '22', 'col', fa_obj.GetLabel('f')))
        sft.SetKey('quit\nY')
        sft.Run()
        filename=os.path.join(sft.rundir,'COMPLETE.TAB')
        filename2=os.path.join(sft.rundir,'COMPLETE2.TAB')
        shutil.move(filename, filename2)
        sft.ClearAllKeys()
        sft.runname=sft.name+'_sig'
        sft.SetKey('read', '"'+fa_obj.GetFileName('mtz')+'"')
        sft.SetKey('reduce')
        sft.SetKey('select', ('col', fa_obj.GetLabel('f'), 'present'))
        sft.SetKey('select', ('col', fa_obj.GetLabel('sigf'), '>', '0'))
        sft.SetKey('purge\nY')
        sft.SetKey('calc', ('R','col','_delpm_','=','col','"'+fa_obj.GetLabel('f')+'"','abs'))
        sft.SetKey('calc', ('R','col','_delsigratio_','=','col','_delpm_','col','"'+fa_obj.GetLabel('sigf')+'"','/') )
        fp_obj,fm_obj=self.inp.Get('fsigf',typ='plus',filetype='mtz'),self.inp.Get('fsigf',typ='minus',filetype='mtz')
        #if fp_obj and fm_obj:
        #  sft.SetKey('calc', ('R','col','_sigpratio_','=','col','"'+fp_obj.GetLabel('f')+'"','col','"'+fp_obj.GetLabel('sigf')+'"','/'))
        #  sft.SetKey('calc', ('R','col','_sigmratio_','=','col','"'+fm_obj.GetLabel('f')+'"','col','"'+fm_obj.GetLabel('sigf')+'"','/'))
        #  sft.SetKey('absent', ('col', fa_obj.GetLabel('f'), 'if','col','_sigpratio_','<','5'))
        #  sft.SetKey('absent', ('col', fa_obj.GetLabel('f'), 'if','col','_sigmratio_','<','5'))
#        sft.SetKey('plot', ('col', '_delsigratio_', 'resol'))
        sft.SetKey('absent', ('col', fa_obj.GetLabel('f'), 'if','col','_delsigratio_','<','2.75'))
        sft.SetKey('complete', ('table', 'shells', '22', 'col', fa_obj.GetLabel('f')))
        sft.SetKey('quit\nY')
        sft.Run()
        self.programs.remove(sft)
        if os.path.isfile(filename) and os.path.isfile(filename2):
          with open(filename) as f:
            with open(filename2) as g:
              resbin = []
              for line in f:
                resbin.append([float(v) for v in line.split()])
              for i,line in enumerate(g):
                if float(line.split()[-1])>0.01:
                  resbin[i][-1] /= (float(line.split()[-1])/100.0)
                else:
                  resbin[i][-1] = -1
          i=0
          for r in resbin[:]:
            if r[-1]<0:
              resbin.pop(i)
            else:
              i+=1
          # smooth
          s=[]
          for i,r in enumerate(resbin):
            #print( r[0],r[-1] )
            smooth, tot_smooth = 0., 0.
            for j in range(max(0,i-2),min(i+3,len(resbin))):
              smooth += (1.0-0.425*abs(i-j))*resbin[j][-1]
              tot_smooth += 1.0-0.425*abs(i-j)
            s.append(smooth/tot_smooth)
          for i,r in enumerate(resbin):
            r[-1] = s[i]
          # 4.5A derived threshold
          #resder=4.5
          resder = 3.0+self.inp.Get(filetype='mtz').GetResolution(self)
          #print ('ggg',resder)
          thres,tot_w,part = 0., 0., 1.
          for r in resbin:
            if r[1]<resder:
              part = (r[0]-resder)/(r[0]-r[1])
            thres += part*r[2]*r[3]
            tot_w += part*r[2]
            if part<1.0:
              break
          # add bonus for very low res.
          for r in resbin:
            r[-1] *= max(1.0,2*(r[0]+r[1])-2.5*resder)
            #print (r[0],r[1],r[-1])
          thres = thres/tot_w/7.0
          #print ('thres',thres)
          cutoff,cutoff0 = resbin[0][1],resbin[0][0]
          for r in resbin:
            # two subsequent bins above the threshold
            #if r[-1]>thres: 
            #  if abs(r[0]-cutoff0)<1e-5:
            #    cutoff = r[1]
            #  cutoff0 = r[1]
            # two subsequent bins below the threshold
            if r[-1]>thres:
              cutoff = r[0]#(r[1]+r[0])/2
            else:
              if abs(r[0]-cutoff0)<1e-5:
                break
              cutoff0 = r[1]
          else:
            cutoff=r[1]
          for r in resbin:
            # correction if no signal in higher res
            if r[-1]<0.1*thres:
              if (r[0]+0.5)>cutoff:
                cutoff = r[0]+0.5
              break
          if verbose:
            print('Anom. cutoff:',cutoff)
        else:
          common.Warning('Unexpected sftools problem when determining resolution cutoff.')
    elif not halfang:
      common.Warning('Resolution cutoff could not be determined.')
    if use_faest_08:
      cutoff=min(cutoff,cut08)
    return cutoff

  def RunPreprocess(self,*args,**kwargs):
    process.RunPreprocess(self,*args,**kwargs)
    prog = self.GetProg(supported=True)
    # lets convert fa to e-values
    if not self.inp.Get('fsigf',typ=('fa','delta-anom'),col='e') and prog.nick!='shelxd':
      for fa in self.inp.GetAll('fsigf',typ='fa'):
        fa_obj = self.inp.Get('fsigf',typ='fa',filetype='mtz', inp_cont=self.inp.Get('fsigf',typ='average',filetype='mtz',xname=fa.xname,try_convert=False))
        if fa_obj:
          ecalc=self.AddProg('ecalc')
          ecalc.inp.Add(fa_obj)
          ecalc.Run()
          self.inp.Add(ecalc.out.Get('fsigf',typ='fa'))
          self.programs.remove(ecalc)
          break
      else:
        common.Warning('Calculation of E-values by ecalc from Fa values failed.')
    # run pmf for crunch2 (no special process defined as of now)
    if self.GetParam('patterson_seed') and prog.nick=='crunch2':
      pmf=self.AddProg('pmf',propagate_out=False)
      pmf.SetKey('TRY',(1,150))
      try:
        pmf.Run()
      except (pmf.ProgramRunError, common.CrankError) as e:
        common.Warning("Will start without seeding as PMF crashed with message: {0}".format(str(e)))
        if hasattr(sys,'exc_clear'): sys.exc_clear()
      else:
        prog.inp.Add(pmf.out.Get('model'))
      self.programs.remove(pmf)
    # a hack needed if automatic res. cutoff is used for shelxd... or if pats is disabled.
    if prog.nick=='shelxd' and ( (self.GetParam('patterson_seed') is False and not prog.IsKey('PATS')) or \
       (self.IsParam('high_res_cutoff') and not self.IsInputtedParam('high_res_cutoff') and not prog.IsKey('SHEL')) ):
      ins = self.inp.Get(filetype='ins').GetFileName()
      if ins:
        with open(ins) as f:
          lines = f.readlines()
        with open(ins,'w') as g:
          for line in lines:
            if line.startswith('SHEL') and self.IsParam('high_res_cutoff') and \
               not self.IsInputtedParam('high_res_cutoff') and not prog.IsKey('SHEL'):
              line='SHEL 999.0 {0}\n'.format(self.GetParam('high_res_cutoff'))
            if not line.startswith('PATS') or self.GetParam('patterson_seed') is not False or prog.IsKey('PATS'):
              g.write(line)
      else:
        common.Error('No ins file inputted to shelxd.')
    #if prog.nick=='prasa':
    #  self.hist_bin_size, self.hist_num_bins, self.hist_min = 0.5, 100, -10.0
    #  self.cfom_hist = [ [(float(i)+self.hist_min-0.5)*self.hist_bin_size, 0]  for i in range(self.hist_num_bins) ]
    self.prog=prog
    if self.stop_file and self.GetCrankParent() and self.GetCrankParent().rundir:
      self.stop_file = os.path.join(self.GetCrankParent().rundir,self.stop_file)

  def Cutoffs(self):
    # multiplied by 100 just so that range can be used
    cut0, radius, step = int(self.GetParam('high_res_cutoff')*100), int(self.GetParam('high_res_cutoff_radius')*100), int(self.GetParam('high_res_cutoff_step')*100)
    reso = self.inp.Get(filetype='mtz').GetResolution(self)*100
    #print(reso,step,reso-step)
    cutoffs = [float(format(cut*0.01,"2.2f")) for cut in range(cut0,cut0+radius+1,step)] + [float(format(cut*0.01,"2.2f")) for cut in range(cut0-step,cut0-radius-1,-step) if cut>=reso-step]
    return cutoffs

  def RunBody(self,*args,**kwargs):
    prog0=self.GetProg(supported=True)
    if self.GetParam('high_res_cutoff_radius') and self.GetParam('high_res_cutoff'):
      cutoffs=self.Cutoffs()
      tot_trials=self.GetParam('num_trials')
      self.Info('The following high resolution cutoffs will be tried: '+', '.join(format(r,"2.2f") for r in cutoffs))
      if prog0.nick=='prasa':
        #prog0.SetKey('statrescut', max(max(cutoffs), reso+2))
        prog0.SetKey('minstatrescut', min(cutoffs))
        prog0.SetKey('maxstatrescut', max(cutoffs))
        prog0.SetKey('ntrials', int(tot_trials/len(cutoffs)))
        # generate reference histogram
#        prasa_refer=self.AddProgCopy(prog0, deeper_copy=True)
#        prasa_refer.SetKey('generaterefdist', 1)
#        prasa_refer.SetKey('histmatch', 2)
#        prasa_refer.SetKey('mtzout', 'refdist.mtz')
#        prasa_refer.SetKey('rescut',max(cutoffs))
#        prasa_refer.Run( rundir=os.path.join(self.rundir,prasa_refer.nick+'_refdist') )
#        prog0.SetKey('refmtzin', os.path.join(self.rundir,prasa_refer.nick+'_refdist','refdist.mtz'))
#        prog0.SetKey('colin-ref-fo', '*/*/[pras.E_sigE.E,pras.E_sigE.sigE]')
#        self.programs.remove(prasa_refer)
      for i,cutoff in enumerate(cutoffs):
        prog=self.AddProgCopy(prog0, deeper_copy=True)
        self.score = -1000.0
        if i==0:  self.prog=prog
        prog.jobnum=i
        if prog.nick=='prasa':
          prog.SetKey('rescut', cutoff, keep_previous=False)
        self.Info('Running with high resolution cutoff of {0}.'.format(cutoff))
        self.PrintActualLogGraph(cutoff=cutoff)
        prog.Run( rundir=os.path.join(self.rundir,prog.nick+'_'+str(cutoff)) )
        if self.stop_message:
          break
      self.programs.remove(prog0)
    else:
      try:
        prog0.Run()#restore=False)
      except program.from_name('shelxd',None).Exception_ShelxD_TooSmall:
        self.Info('ShelxD arrays too small - trying with -L20')
        prog0.ClearAllArgs()
        prog0.SetArg('L20')
        prog0.Run()

  def RunPostprocess(self,restore=True,*args,**kwargs):
    prg=self.GetProg(supported=True)
    outmod=self.out.Get('model',typ='substr')
    if not outmod or not outmod.GetFileName() or not os.path.isfile(outmod.GetFileName()):
      self.stop_message='{} failed - no substructure solution could be found.  See {} for more information.'.format(prg.name,prg.GetLogFileName())
      for mod in self.out.GetAll('model',typ='substr'):
        mod.DetachFile(outmod.GetFileName())
      self.PrintActualLogGraph(finished=True)
      common.Error(self.stop_message, nosuccess=True)
    # for crunch2 we need to adjust the filename!
    if prg.nick=='crunch2':
      if not self.trial:
        common.Warning('Could not pick good model. Detection likely did not work.')
        self.trial = self.fom_unopt_all.index(max(self.fom_unopt_all))
      if self.trial:
        #outmod.SetFile( outmod.GetFileName().replace('trial1a','trial'+str(self.trial)+'a'), 'pdb' )
        outmod.SetFile( self.GetProg(supported=True).GetTrialPdb(self.trial), 'pdb' )
      else:
        common.Error('Could not pick good model. Detection failed.', nosuccess=True)
    # taking the solution with largest adjusted score for prasa
    if prg.nick=='prasa':
      if self.IsTrueOrNoneParam('optimize_sol') or self.GetParam('optimize_sol') in (1,2):
        ref_all_sol = True if self.GetParam('optimize_sol')==2 or \
          (self.GetParam('optimize_sol') is None and self.score_adj*100<self.GetParam('threshold_weak')) else False
        if not ref_all_sol and self.GetParam('optimize_sol') is None: # be careful with solutions with many special pos.
          import gemmi
          struct=gemmi.read_structure(self.prog.out.Get('model').GetFileName())
          spec_pos=[struct.cell.is_special_position(s.atom.pos) for s in struct[0].all() if s.atom.occ>=0.3]
          if spec_pos and sum(spec_pos)/len(spec_pos)>0.35:
            ref_all_sol = True
        if ref_all_sol:
          self.score_cut=-1000
        for ip,prog in enumerate(self.GetProgs('prasa')):
         # optimize for each prog that found some solutions
         if (ref_all_sol or prog is self.prog) and hasattr(prog,'trial'):
          prasa_ref = self.AddProgCopy(prog, deeper_copy=True)
          prasa_ref.runname = prasa_ref.name+'_ref'
          self.refsol = True
          self.score = None
          prasa_ref.SetKey('specialpos', 1) #perhaps disable due to 3og5?
          prasa_ref.SetKey('chargeflip', 2, keep_previous=False)
          #prasa_ref.SetKey('histmatch', 1, keep_previous=False)
          if os.path.isfile(os.path.join(prog.rundir,'stop_prasa')):
            os.remove( os.path.join(prog.rundir,'stop_prasa') )
          if ref_all_sol:
            import glob
            for pdb in glob.glob(os.path.join(prog.rundir,'*.pdb_*')):
              prasa_ref.SetKey('pdbin', pdb)
          else:
            prasa_ref.SetKey('pdbin', self.prog.out.Get('model').GetFileName())
            if prog.trial:
              prasa_ref.SetKey('pdbin', self.prog.GetTrialPdb(prog.trial))
            if hasattr(prog,'trial_stop') and prog.trial!=prog.trial_stop:
              prasa_ref.SetKey('pdbin', self.prog.GetTrialPdb(prog.trial_stop))
            self.prog=prasa_ref
            #prasa_ref.SetKey('chargeflip',0)
          prasa_ref.SetKey('shannon', 2.2) #larger needed eg for 3uot (issues with large for 7bmv?, 6fmw not with 2?)
          prasa_ref.SetKey('statcyc', 20, keep_previous=False) # larger needed eg for 7bmy
          prasa_ref.SetKey('statcycskip', 1), prasa_ref.SetKey('ncycles', 45, keep_previous=False)
          prasa_ref.SetKey('ntrials', False, keep_previous=False)#, prasa_ref.SetKey('addhalfstat', 0)
          if prasa_ref.GetKey('delta'):
            prasa_ref.SetKey('delta', prasa_ref.GetKey('delta')+0.1, keep_previous=False)
          elif prasa_ref.GetKey('chargeflip') and prasa_ref.GetKey('chargeflip')==2:
            prasa_ref.SetKey('delta', 1.81, keep_previous=False) # lower helps 7waa, larger helps uncorrected 2w5o? and 3tx3
          else:
            prasa_ref.SetKey('delta', -0.1, keep_previous=False)
          prasa_ref.outfilename['pdb'] = 'prasa_ref.pdb'
          self.Info('Refining solutions from high resol. cutoff {0}'.format(prasa_ref.GetKey('rescut') or ''))
          self.PrintActualLogGraph(cutoff=prasa_ref.GetKey('rescut'),optim=ip)
          prasa_ref.Run()
      outmod=self.prog.out.Get('model')
      if self.trial_adj:
        outmod=self.prog_adj.out.Get('model')
        outmod.SetFile( self.prog_adj.GetTrialPdb(self.trial_adj), 'pdb' )
        self.Info('Solution from trial {0} (adjusted score {1}) taken.'.format(self.trial_adj,self.score_adj))
      elif self.score_cut>0.0:
        if self.prog.GetKey('rescut'):
          self.Info('Solution from cutoff {2} (cutoff CCrange {1}, trial {0}) taken.'.format(self.trial,self.score_cut,self.prog.GetKey('rescut')))
        else:
          self.Info('Solution with CC {1}, trial {0} taken.'.format(self.trial,self.score))
          self.GetParam('high_res_cutoff_step')
    # "fix" the substructure pdb file
    fixsubstr=self.AddProcess('fixsubstrpdb')
    fixsubstr.inp.Set(outmod)
    fixsubstr.Run()
    # print the final loggraph
    self.PrintActualLogGraph(finished=True)
    # if weak phases and phdmmb follows then set the thorough building mode (perhaps this should be in manager?)
    crank=self.GetCrankParent()
    if crank:
      score_here = self.score_adj*100  if self.GetProg(supported=True).nick=='prasa' else self.score
      if score_here<self.GetParam('threshold_weak') and crank.prep.GetProcess('phdmmb') and not crank.prep.GetProcess('phdmmb').IsParam('thorough_build'):
        crank.prep.GetProcess('phdmmb').SetParam('thorough_build')
      if score_here<self.GetParam('threshold_weak')/2.:
        self.guess=0
      elif self.guess is None:
        self.guess=1
      if crank.prep.GetProcess('refatompick'):
        if score_here>=self.GetParam('threshold_weak'):
          if not crank.prep.GetProcess('refatompick').IsParam('num_iter'):
            crank.prep.GetProcess('refatompick').SetParam('num_iter',1)
        elif self.GetProg(supported=True).nick=='shelxd' and not crank.prep.GetProcess('refatompick').IsParam('bfactor'):
          crank.prep.GetProcess('refatompick').SetParam('bfactor',20)
        elif self.GetProg(supported=True).nick=='prasa' and not crank.prep.GetProcess('refatompick').IsParam('bfactor'):
          crank.prep.GetProcess('refatompick').SetParam('bfactor',30)
        #if not crank.prep.GetProcess('refatompick').IsParam('res_cut'):
         # crank.prep.GetProcess('refatompick').SetParam('res_cut',self.EstimateCutOff(use_faest_08=True))
    self.GetParam('high_res_cutoff_cchalf')
    process.RunPostprocess(self,restore,*args,**kwargs)

  def RedrawHistogram(self):
    self.rv_hist.Reset()
    max_hind = max([i for i,v in enumerate(self.cfom_hist) if v[1]>0])+1 if any(list(zip(*self.cfom_hist))[1]) else 0
    min_hind = min([i for i,v in enumerate(self.cfom_hist) if v[1]>0])   if any(list(zip(*self.cfom_hist))[1]) else 0
    self.rv_hist.Data( [v[0] for v in self.cfom_hist[min_hind:max_hind]], [math.log(v[1]+1,2) for v in self.cfom_hist[min_hind:max_hind]] )

  def PrintActualLogGraph(self, finished=False, cutoff=None, optim=None):
    # print loggraph info for shelxd and prasa/crunch2
    prog=self.GetProg(supported=True)
    if cutoff: self.cutoff2present=cutoff
    if finished and not self.stop_message:
      self.stop_message='Maximum number of trials reached.'
    if prog.nick=='shelxd':
      sc,sc2,scb=('CCweak','CC'),"correlation coefficients",'CFOM'
    elif prog.nick=='crunch2':
      sc,sc2,scb=('FOM','FOM_unoptimized'),"figures of merit",'FOM'
    elif prog.nick=='prasa':
      sc,sc2,scb=('CCrange','CC'),"correlation coefficients",'CCrange'
    # writing out loggraph takes time & causes slowdown.
    # with larger number of cycles, we need to optimize by only writing at certain cycles.
    if self.opened_loggraph and not prog.nick=='prasa' and \
       ( finished or len(self.stats_table)<2000 or len(self.stats_table)%500==0 ):
      self.GetLogGraphHandle().seek(0,0)
      self.LGInfo(self.GetCCP4Header())
      self.LGInfo(self.GetLogGraphPrograms())
      self.LGInfo('\n$TABLE : {0} {1} versus trial:\n$SCATTER'.format(prog.name,sc[0]))
      self.LGInfo(' : {0} per trial :A:1,2,3:\n$$'.format(' and '.join(sc)))
      self.LGInfo('Trial       {0}$$ $$'.format('        '.join(sc)))
      for trial,ccw,cc,fom,fom_unopt in self.stats_table:
        if prog.nick=='shelxd':
          self.LGInfo('{0:5}  {1:10.2f} {2:10.2f}'.format(trial,ccw,cc))
        elif prog.nick=='crunch2':
          self.LGInfo('{0:5}  {1:10.3f} {2:10.3f}'.format(trial,fom,fom_unopt))
      self.LGInfo('$$\n\n')
      if finished:
        self.LGInfo('$TEXT:Result: $$ Final result $$')
        self.LGInfo(self.stop_message)
        self.LGInfo('Maximum {0}: {1}, from trial {2}.\n$$'.format(scb, self.score_cut, self.trial))
    if self.rv_report is not None:
      if not finished:
        trial,ccw,cc,fom,fom_unopt = self.stats_table[-1]  if self.stats_table  else (0,0,0,0,0)
        # for shelxd we restrict the number of points as rvapi is slow or gets stuck dealing with many hundreds of points
        # we do this by re-selecting the points each resel_cycle cycles, thus keeping the max. displayed below roughly 500 (per line)
        if prog.nick in ('shelxd','prasa'):
          if prog.nick=='prasa' and cc and not ccw:  ccw=cc
          tr_half, cfom = trial/self.resel_init_half, cc+ccw if prog.nick=='shelxd' else ccw
        if optim==0:
          self.rv_plot_against = self.rv_plot_against.parent.Plot( '{0} against {1} optimized'.format(sc[1],sc[0]), sc[0]+' [x100]', sc[1]+' [x100]', legendloc='nw', xmin=0., ymin=0. )
        if not hasattr(self,'rv_plot') or (trial>=self.resel_cycle and prog.nick=='shelxd'):
          if trial>=self.resel_cycle and prog.nick in ('shelxd','prasa'):
            self.rv_plot.parent.parent.Remove()
          plot_name = '{1} and {2} vs trial'.format(prog.name,sc[0],sc[1])
          if prog.nick in ('shelxd','prasa'):
            self.rv_plot_against = self.rv_report.Plot( '{0} against {1}'.format(sc[1],sc[0]), sc[0]+' [x100]', sc[1]+' [x100]', block="Correlation coefficients", legendloc='nw', xmin=0., ymin=0. )
            self.rv_plot_hist = self.rv_plot_against.parent.Plot( 'Distribution of {0} from trials'.format(scb), scb+' [x100]', "Number of trials (log scale)", xmin=0., ymin=0., inty=True )
            #max_hind = max([i for i,v in enumerate(self.cfom_hist) if v[1]>0])+1 if any(zip(*self.cfom_hist)[1]) else 0
            #min_hind = min([i for i,v in enumerate(self.cfom_hist) if v[1]>0])   if any(zip(*self.cfom_hist)[1]) else 0
            #self.rv_hist = self.rv_plot_hist.PlotLine( ["x",[v[0] for v in self.cfom_hist[min_hind:max_hind]]], \
              #                  ["CFOM distr.",[v[1] for v in self.cfom_hist[min_hind:max_hind]]], style='bars' )
            self.rv_hist = self.rv_plot_hist.PlotLine( ["x"], ["{0} distr.".format(scb)], style='bars' )
            self.rv_against = self.rv_plot_against.PlotLine( ["x"], ["cutoff {0}".format(self.cutoff2present if self.cutoff2present else "CC")], style="off", size=1.2 )
            if prog.nick=='prasa':
              self.rv_plot = None
            else:
              self.rv_plot = self.rv_plot_against.parent.Plot( plot_name, "Trial", "{0} and {1}".format(sc[0],sc[1]), intx=True )
            #self.rv_agc={}  # to keep the cutoff plotlines
            #if cutoff:   self.rv_agc[cutoff] = self.rv_against
          else:
            self.rv_plot = self.rv_report.Plot( plot_name, "Trial", "{0} and {1} [[x100]]".format(sc[0],sc[1]), block="{0} {1}".format(prog.name,sc2), intx=True )
          if prog.nick not in ('prasa'):
            self.rv_stat1 = self.rv_plot.PlotLine( ["x"], [sc[0]], style='off', size=1.5)
            self.rv_stat2 = self.rv_plot.PlotLine( ["x"], [sc[1]], style='off', size=1.5)
          if trial>=self.resel_cycle and prog.nick in ('shelxd','prasa'):
            # new selection is made and plotted.
            stat_all = [s[1]+s[2] if prog.nick=='shelxd' else s[1] for s in self.stats_table]
            self.thr_up, self.thr_down = heapq.nlargest(20,stat_all)[-1], heapq.nsmallest(20,stat_all)[-1]
            stat_sel = [s for i,s in enumerate(self.stats_table) if i%tr_half==0 or stat_all[i]>=self.thr_up or stat_all[i]<self.thr_down]
            stat_sel = list(zip(*stat_sel))
            if prog.nick != 'prasa':
              self.rv_stat1.Data(stat_sel[0],stat_sel[1]), self.rv_stat2.Data(stat_sel[0],stat_sel[2])
            ##self.rv_against.Data(stat_sel[1],stat_sel[2])
            stat_full = list(zip(*self.stats_table))
            self.rv_against.Data(stat_full[1],stat_full[2])
            ####self.resel_cycle*=2
            self.resel_cycle += self.resel_init_half
        elif cutoff or optim is not None:
          #if cutoff not in self.rv_agc:
            color = None
            if len(self.rv_plot_against.children)==2:  color='darkorchid'
            if len(self.rv_plot_against.children)==3:  color='darkgreen'
            if len(self.rv_plot_against.children)==4:  color='brown'
            self.rv_against = self.rv_plot_against.PlotLine( ["x"], ["cutoff {0}{1}".format(cutoff," optim." if optim is not None else "")], style="off", size=1.2, color=color )
            #self.rv_against = self.rv_agc[cutoff] = self.rv_plot_against.PlotLine( ["x"], ["optim." if optim is not None else "cutoff {0}".format(cutoff)], style="off", size=1.2, color=color, marker="diamond" if optim is not None else "" )
          #else:
          #  self.rv_against = self.rv_agc[cutoff]
        if self.stats_table and (ccw or fom):
          if prog.nick in ('shelxd') and (trial<self.resel_init_half or trial%tr_half==0 or cfom>=self.thr_up or cfom<self.thr_down):
            self.rv_stat1.Data(trial,ccw,x_type=int), self.rv_stat2.Data(trial,cc,x_type=int)
          if prog.nick in ('shelxd','prasa'):
            self.rv_against.Data(ccw,cc)
          elif prog.nick=='crunch2':
            self.rv_stat1.Data(trial,fom,x_type=int), self.rv_stat2.Data(trial,fom_unopt,x_type=int, flush=True)
          if prog.nick in ('shelxd','prasa'):
            hist_i = max(0, int((cfom/self.hist_bin_size-self.hist_min)))
            self.cfom_hist[hist_i][1]+=1
            #self.rv_hist.Data(self.cfom_hist[hist_i][0], self.cfom_hist[hist_i][1])
            if trial%50==0:
              self.RedrawHistogram()
            crvapi.Flush(timing_restr=3.0 if len(self.stats_table)<9000 else len(self.stats_table)/3000)
      else:  # if finished
        if prog.nick in ('shelxd','prasa') and hasattr(self,'rv_hist'):
          self.RedrawHistogram()
          crvapi.Flush()
        occ=None
        if self.out.Get('model',filetype='pdb') and os.path.isfile(self.out.Get('model',filetype='pdb').GetFileName('pdb')):
          with open(self.out.Get('model',filetype='pdb').GetFileName('pdb')) as f:
            occ = program.from_name('pdbcur',None).GetStatGrep('occup',from_str=f.read())
          if occ:
            self.occ_plot_hist = self.rv_report.Plot( 'Occupancies from the best trial', "Atom number", "Estimated occupancy", \
                                                       block="Estimated atom occupancies", ymin=0., ymax=1., xmin=0.5, xmax=len(occ)+0.5, intx=True )
            custom_x_ticks = False if len(occ)>10 else True
            if custom_x_ticks:
              self.occ_plot_hist.CustomTick('x',0,' ')
            self.occ = self.occ_plot_hist.PlotLine( ["x",[o+1 for o in range(len(occ))]], ["Atom occupancies",occ], style='bars', custom_x_tick=custom_x_ticks )
            if custom_x_ticks:
              self.occ_plot_hist.CustomTick('x',len(occ)+1,' ')
            occ_thr = float(occ[0])*0.25
            self.shadow = self.occ_plot_hist.PlotLine( ["x",(0,len(occ)+1)], ["Suspiciously low occupancy area",(occ_thr,occ_thr)], fill=1 )
        self.rv_report.Text("<i><small>Stop:</small></i>")
        button=crvapi.init_meta["help_btn_template"].replace('.html','.html#term-what-does-the-threshold-in-substructure-determination-mean') if "help_btn_template" in crvapi.init_meta else ""
        self.rv_report.Text('&emsp;'+self.stop_message+button)
        self.rv_report.Text("<i><small>Result:</small></i>") # title
        button=crvapi.init_meta["help_btn_template"].replace('.html','.html#term-what-is-{}-in-the-{}-output'.format(scb,prog.nick)) if "help_btn_template" in crvapi.init_meta else ""
        judgement='&emsp;Judging from the obtained scores, the substructure is likely to be correctly determined.'
        score_here = self.score_adj*100  if self.GetProg(supported=True).nick=='prasa' else self.score
        rat = 4.0 if self.GetProg(supported=True).nick=='prasa' else 2.0
        if score_here<self.GetParam('threshold_weak')/rat:
          judgement='&emsp;Judging just from the obtained scores, there is a large chance that substructure has not been correctly determined.'
        elif score_here<self.GetParam('threshold_weak'):
          judgement='&emsp;From the obtained scores, it is unclear whether the substructure is correct or not.'
        elif score_here<self.GetParam('threshold_stop'):
          judgement='&emsp;Judging from the obtained scores, there is a reasonable chance that the substructure has been correctly determined.'
        self.rv_report.Text(judgement)
        self.rv_report.Text('&emsp;Maximum {0}{3}: {1:.2f} occurring at Trial {2}.'.format(scb, self.score_cut, self.trial, button))
        if occ:
          searched='.'
          if prog.nick=='shelxd' and prog.GetKey('FIND'):
            searched = ' (searched for {0} peaks'.format(prog.GetKey('FIND'))
            searched+= ' of which {0} modelled as disulphides)'.format(prog.GetKey('DSUL')) if prog.GetKey('DSUL') else ')'
          self.rv_report.Text('&emsp;The substructure has {1} atoms with occupancy of at least 25%, {0} in total{2}'.format(
            len(occ), len([o for o in occ if float(o)>=0.25]), searched))
        if not self.parent_process or not self.parent_process.ccp4i2:
          self.rv_report.Text("<BR><i><small>Output:</small></i>")
          self.rv_report.DataFile(self.out.Get('model').GetFileName(), "xyz", "<small>Substructure</small>", flush=True)
          if self.inp.Get(filetype='ins'):
            self.rv_report.Text("<BR><i><small>Input:</small></i>")
            self.rv_report.DataFile(self.inp.Get(filetype='ins').GetFileName(), "text", "<small>SHELXD input parameters</small>")


  def UpdateInfo(self,line,popen=None,empty=False,prog=None):
    # in shelx all comes in one line; in crunch not and we only have (optimized) fom outputted
    score, ccw, cc, fom_unopt, score_cut = None, None, None, None, None
    if prog is None:
      prog=self.GetProg(supported=True)
    trial=prog.GetStatGrep('try', from_str=line)
    if prog.nick=='prasa':
      patfom=0.
    else:
      patfom=prog.GetStatGrep('patfom', from_str=line)
    # early stop requested by user
    if (hasattr(self,'ccp4i2job') and self.ccp4i2job.testForInterrupt()) or os.path.isfile(self.stop_file):
      self.stop_message='Stopping early on user request!'
      if os.path.isfile(self.stop_file):
        os.remove(self.stop_file)
    stat_str='Best stats improved'
    # shelx part
    if prog.nick=='shelxd':
      ccw=prog.GetStatGrep('ccweak', from_str=line)
      cc=prog.GetStatGrep('cc', from_str=line)
      if ccw and cc:
        if ccw=='NaN' or cc=='NaN':
          self.stop_message='{} failed to obtain any solutions from this data.'.format(prog.name)
          self.PrintActualLogGraph(finished=True)
          common.Error(self.stop_message, nosuccess=True)
        score=ccw+cc
      if trial and score and score>self.score:
        stats_info=['corr.coef. {0}'.format(cc), 'corr.coef. weak {0}'.format(ccw), 'CFOM {0:.1f}'.format(cc+ccw)]
        self.Info('{0} in trial {1}: {2}'.format(stat_str,trial,', '.join(stats_info)))
    # prasa part
    if prog.nick=='prasa':
      cc=prog.GetStatGrep('cc', from_str=line)
      score=cc
      ccw=cc_range=prog.GetStatGrep('cc_range', from_str=line)
      if cc_range:  score_cut=cc_range
      #if trial and score and score>self.score and (not self.score_cut or score_cut>=self.score_cut):
      #  stats_info=['corr.coef. {0}'.format(cc),]
      #  if score_cut: stats_info.append( 'corr.coef.range {0}'.format(score_cut) )
      #  self.Info('{0} in trial {1}: {2}'.format(stat_str,trial,', '.join(stats_info)))
    # crunch part
    if prog.nick=='crunch2':
      fom_unopt=prog.GetStatGrep('fom_unopt', from_str=line)
      if trial is not None:
        self.trial_temp=trial
        self.fom_unopt=-1000.
      elif self.trial_temp and patfom is not None:
        trial=self.trial_temp
        score=patfom
        self.Info('FOM after optimization in trial {0}: {1}'.format(trial,patfom))
      elif self.trial_temp and fom_unopt is not None and fom_unopt>self.fom_unopt:
        self.fom_unopt=fom_unopt
        if len(self.fom_unopt_all)<=self.trial_temp:
          self.fom_unopt_all.append(fom_unopt)
        self.fom_unopt_all[self.trial_temp] = fom_unopt
      fom_unopt=self.fom_unopt
    if 'cutoff' in prog.stat and (not self.GetParam('high_res_cutoff_radius') or not self.GetParam('high_res_cutoff')):
      cutoff=prog.GetStatGrep('cutoff', from_str=line)
      if cutoff:
        self.PrintActualLogGraph(cutoff=cutoff)
    # save the best score and trial; save any stats to the table of stats; check thresholds; update loggraph
    prasa_check_thres = 10.
    score_now = score_cut  if score_cut  else score
    is_jobnum = hasattr(prog,'jobnum') and prog.jobnum>0
    if score_now and trial:
      if not self.stop_message:
        if self.GetParam('threshold_stop') and score_now>=self.GetParam('threshold_stop') and \
           (prog.nick not in ('shelxd','prasa') or trial>=45 or is_jobnum):
          self.stop_message='Specified score threshold ({0}) reached.'.format(self.GetParam('threshold_stop'))
          self.guess=2
        elif self.GetParam('threshold_stop') and prog.nick=='crunch2':
          self.CheckExtraStopConditions(score,trial,prog,popen)
        elif prog.nick=='prasa' and ( \
             (trial==50 and self.check_solution[0]==-1 and self.score_cut>prasa_check_thres and not self.refsol)):
          self.check_solution=(self.trial,self.score_cut)
      if self.score is None:  self.score=-1000
      score_sum = score_cut + score if score_cut else score
      if score_sum>=self.score and (prog.nick!='prasa' or score_cut is None or score_cut>=0.8*self.score_cut):
        if prog.nick=='prasa' and (self.prog is prog or score_cut is None or score_cut>=0.97*self.score_cut):
          stats_info=['corr.coef. {0}'.format(cc),]
          if score_cut: stats_info.append( 'corr.coef.range {0}'.format(score_cut) )
          self.Info('{0} in trial {1}: {2}'.format(stat_str,trial,', '.join(stats_info)))
          if score_now>prasa_check_thres and (trial>45 or is_jobnum) and not self.stop_message and not self.refsol and \
             (score_sum>1.01*self.score or score_sum>self.score+1.0 or self.solcheck_skipped):
            self.check_solution=(trial,score_now)
        self.score=score_sum
        prog.trial=trial
      if score_now>=self.score_cut and score_sum>=0.85*self.score:
        self.score_cut=score_now
        self.prog=prog
        self.trial=prog.trial
#      if score>self.score and (not self.score_cut or score_cut>self.score_cut):
#        self.score=score
#        self.score_cut=score_cut
#        self.trial=trial
#        self.prog=prog
#        if prog.nick=='prasa' and score_now>prasa_check_thres and (trial>20 or is_jobnum) and not self.stop_message:
#          self.check_solution=(trial,score_now)
      if score<self.min_score:
        self.min_score=score
      self.stats_table.append( (trial,ccw,cc,patfom,fom_unopt) )
      self.PrintActualLogGraph()
    # let's stop early; with crunch2, we need to wait till the trial is finished
    if self.stop_message and not self.stop_sent and \
       (prog.nick!='crunch2' or prog.GetStatGrep('cyc_in_try', self.trial, from_str=line)):
      self.Info(' Stopping: '+self.stop_message)
      prog.Stop(popen)
      self.stop_sent = True

  def CheckExtraStopConditions(self,score,trial,prog,popen):
    # relative score threshold
    if score and score>0.3 and score>2.25*self.min_score:
      self.stop_message='Relative score threshold reached.'
    # running bp3 in check mode to preliminary stop on bp3 Luzzati threshold
    elif score and score>0.25 and self.trial_temp%3==1:
      fixsubstr=self.AddProcess('fixsubstrpdb',propagate_out=False)
      mod = fixsubstr.inp.AddCopy(prog.out.Get('model'))
      #mod.SetFile( self.out.Get('model').GetFileName().replace('1',str(self.trial_temp)), 'pdb' )
      mod.SetFile( prog.GetTrialPdb(self.trial_temp), 'pdb' )
      fixsubstr.Run()
      phas=self.AddProcess('phas',propagate_out=False)
      phas.inp.Add(fixsubstr.out.Get('model'))
      phas.SetParam('target', self.GetParam('target'))
      phas.AddProg('bp3').SetKey('check')
      try:
        phas.Run()
      except (phas.ProgramRunError, common.CrankError) as e:
        common.Warning('Could not test solution by quick phasing: {0}'.format(e))
      else:
        if phas.GetProg('bp3').GetStat('luzzati')>0.7:
          self.stop_message='BP3 Luzzati parameter threshold reached.'
      self.processes.remove(phas)
      self.processes.remove(fixsubstr)

  def CheckExtraStopConditionsPrasa(self,score,trial,prog,popen):
    self.Info('Fast testing of solution from trial {0} with CC of {1} follows.'.format(trial,score))
    trial_sol=prog.out.GetCopy('model')
    trial_sol.SetFile( prog.GetTrialPdb(trial) )
    if not os.path.isfile( trial_sol.GetFileName('pdb') ):
      self.Info("Skipping check of trial solution {0}.".format(trial) )
      self.solcheck_skipped=True
      return
    fixsubstr=self.AddProcess('fixsubstrpdb')
#    fixsubstr.inp.Set(self.GetProg(supported=True).out.Get('model'))
    fixsubstr.inp.Set(trial_sol)
    fixsubstr.Run()
#    if not fixsubstr.inp.Get('model',has_atomtypes=True):
#      fixsubstr.inp.Get('model').SetAtomTypes(atomtypes,atomtype1)
#    if fixsubstr.inp.Get('model').GetType()=='unknown':
#      fixsubstr.inp.Get('model').SetType('substr')
#    fixsubstr.Run()
#    if score<20.0:
#      return
#    afm=self.AddProcess('atomsfrommap')
#    afm.inp.Add(prog.out.Get('mapcoef',filetype='map'))
#    afm.SetParam('rms_threshold', 5.0)
#    if self.GetParam('num_atoms'):
#      afm.SetParam('max_new_atoms', self.inp.GetParam('num_atoms')+2)
#    elif self.inp.Get(has_num_atoms=True):
#      afm.SetParam('max_new_atoms', self.inp.Get(has_num_atoms=True).exp_num_atoms+2)
#    afm.Run()
    # this should be done somewhere in afm and also the mtz from shelxc should have x/dname!
#    afm.out.Get('model').xname=self.inp.Get('model',typ='substr').xname
    handdet=self.AddProcess('handdet',propagate_out=False)
    handdet.SetParam('target', self.GetParam('target'))
#    handdet.inp.Add(afm.out.Get('model'))
    handdet.inp.Add(fixsubstr.out.Get('model'))
    handdet.TreatInOutPar()
    handdet.phas.SetParam('catch_output',False),  handdet.phas2.SetParam('catch_output',False)
    handdet.phas.SetParam('cycles', 3),  handdet.phas2.SetParam('cycles', 0)
    handdet.phas.SetParam('occ_cut',0.3),  handdet.phas2.SetParam('occ_cut',0.3)
    handdet.dmf[0].SetParam('dmcyc', 0),  handdet.dmf[1].SetParam('dmcyc', 0)
    handdet.Run()
    score_adj=score*(handdet.phas.fom+handdet.phas2.fom)/2*max(0.017,abs(max(-0.01,handdet.cld[1])-max(-0.01,handdet.cld[0])))#*abs(handdet.score[0]-handdet.score[1])
    self.Info('Combined score of the tested trial is {1}'.format(trial,score_adj))
#    if (handdet.cld[0]*score_adj>0.2 or handdet.cld[1]*score_adj>0.2):
    if score_adj>self.score_adj: #and not self.stop_message:
      self.score_adj = score_adj
      #self.trial_adj = trial
      #self.prog_adj = prog
    if (score_adj*100.)>(self.GetParam('threshold_stop')-8.):
      self.stop_message='CLD*FOM*CC parameter threshold reached.'
      self.guess=2
      prog.trial_stop=trial
      #self.score_adj = score_adj
      #self.trial_adj = trial
#    self.out.Set(afm.out.Get('model'))
    self.out.Set(trial_sol)
    self.processes.remove(handdet)
#    self.processes.remove(afm)
    self.processes.remove(fixsubstr)
    self.solcheck_skipped=False
