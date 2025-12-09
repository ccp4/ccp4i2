#!/usr/bin/python
import os,sys,shutil,subprocess,copy,multiprocessing
from process import process, crvapi
from program import program
import common,inout
import fileinput
par=common.parameter

class mbref(process):
  name="iterative model building and refinement"
  short_name="model tracing with refinement"
  supported_procs=["ref", "mb", "dmfull"]
  supported_progs=["arpwarp"]
  supported_params={}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True, share=True )
  supported_params['bigcyc'] = par( desc='Number of cycles of model building', typ=int )
  supported_params['minbigcyc'] = par( desc='Minimal number of model building cycles (only used for with_dm)', typ=int )
  supported_params['low_res_par'] = par( desc='Use low resolution parameters in refinement and model building.', typ=bool, share=True )
  supported_params['arp_sad'] = par( desc='Use the original ARP/wARP SAD procedure', typ=bool )
  supported_params['catch_output'] = par( desc='Catch program output continously', typ=bool )
  supported_params['with_dm'] = par( desc='Include dens.modif. iterations in building', typ=bool )
  supported_params['start_shelxe'] = par( desc='Use SHELXE in initial model building cycles then use Buccaneer', typ=bool, share=True )
  supported_params['handdet'] = par( desc='Run with both hands (not set by default unless prev.steps failed to determine hand)', typ=bool )

  def RunPreprocess(self, rundir=None, **kwargs):
    process.RunPreprocess(self,rundir,kwargs)
    # get the number of monomers for arpwarp if that may be needed
    if 'arpwarp' in self.programs and not self.inp.Get(has_monomers=True) and self.inp.Get('sequence'):
      matthews=self.AddProcess('matthews', propagate_out=False)
      matthews.Run()
      self.processes.remove(matthews)
    # run shelxc to get the ins file for shelxe
    if self.GetProcess('mb') and self.mb.GetProg('shelxe'):
      shelxc=self.AddProg('shelxc')
      shelxc.Run()
      self.inp.Add(shelxc.out.Get('datafile',filetype='ins'))
    if self.stop_file and self.GetCrankParent() and self.GetCrankParent().rundir:
      self.stop_file = os.path.join(self.GetCrankParent().rundir,self.stop_file)

  def Init(self):
    self.mb_res_all, self.ref_res_all = [], []
    self.stop_message=None
    self.cyc_fin,self.R_thres=None,0.42
    self.result_str = ""
    self.report_R = 1.0
    self.report_Rfree = 1.0
    self.hydrog_no_set = None
    self.queue,self.queue2=None,None
    self.num_proc=-1
    self.stop_file='stop_file'

  def TreatInOutPar(self, set_all_par=False):
    # defaults if no program and no subprocess is specified
    if not self.GetProg(supported=True) and not self.GetProcess(supported=True):
      # maybe we could use arpwarp by default if a) available, b) high res (<2.5) ?
      self.AddProcess('ref')
      self.AddProcess('mb')
    if self.GetProcess(supported=True):
      # assign the child processes to specific class attributes for convenience
      self.ph, self.mb, self.dm = self.GetOrAddProcess('ref'), self.GetOrAddProcess('mb'), None
      # check and set some defaults
      crank=self.GetCrankParent()
      if self.GetParam('with_dm') is None and crank and not crank.GetProcess('phdmmb'):
        self.SetParam('with_dm')
      if self.GetParam('with_dm'):
        if not self.GetParam('bigcyc'):
          self.SetParam('bigcyc',15)
        if not self.GetParam('minbigcyc'):
          self.SetParam('minbigcyc',min(5,self.GetParam('bigcyc')))
        if self.GetParam('bigcyc')<self.GetParam('minbigcyc'):
          self.SetParam('bigcyc', self.GetParam('minbigcyc'))
        self.dm = self.GetOrAddProcess('dmfull')
        if not self.dm.GetParam('dmcyc'):
          self.dm.SetParam('dmcyc',7)
        if not self.dm.inp.Get('mapcoef',custom='mlhl') and self.dm.inp.Get('mapcoef',typ='best'):
          self.dm.inp.AddCopy( self.dm.inp.Get('mapcoef',typ='best') ).custom.append('mlhl')
          self.ph.inp.AddCopy( self.dm.inp.Get('mapcoef',typ='best') ).custom.append('mlhl')
        self.dm.AddProcess('phcomb'),  self.dm.TreatInOutPar()
      if not self.GetParam('bigcyc'):
        self.SetParam('bigcyc',3)
        if not self.inp.Get('sequence'):
          self.SetParam('bigcyc',2)
        if crank and not crank.GetProcess('phdmmb'):
          self.SetParam('bigcyc',5)
      self.ph.SetRunDir()
      self.ph.TreatInOutPar()
      self.mb.SetRunDir()
      self.mb.TreatInOutPar()
      if not self.ph.IsInputtedParam('target') and crank and crank.GetProcess('phdmmb'):
        self.ph.SetParam('target','MLHL')
      # always merge for MLHL (to prevent the originally input HL's be rewritten by the new HL's)
      if self.ph.GetParam('target')=='MLHL':
        self.ph.GetProg().always_merge=True
      #if os.environ['CCP4_LIB'] and os.path.isfile(os.path.join(os.environ['CCP4_LIB'],'data','reference_structures','reference-1ajr.pdb')):
      #  self.mb.GetProg().SetKey('pdbin-ref',os.path.join(os.environ['CCP4_LIB'],'data','reference_structures','reference-1ajr.pdb'))
      #  self.mb.GetProg().SetKey('mtzin-ref',os.path.join(os.environ['CCP4_LIB'],'data','reference_structures','reference-1ajr.mtz'))
      if not self.ph.IsParam('cycles'):
        self.ph.SetParam('cycles',10)
      # start with no hydrogens, add them when R is low (and not at low res - as treated in ref)
      if self.ph.GetProg('refmac') and not self.ph.GetProg('refmac').IsKey('make'):
        self.ph.GetProg('refmac').SetKey('make',('hydrogens','no'))
        self.hydrog_no_set=True
      if not self.mb.GetProg('shelxe') and not self.IsInputtedParam('start_shelxe') and (not crank or not crank.GetProcess('phdmmb')) and \
         self.mb.resol and self.mb.resol<2.15 and program.from_name('shelxe',None).CheckBinary(silent=True):
        self.SetParam('start_shelxe',True)
      if not self.mb.GetProg('shelxe') and self.GetParam('start_shelxe'):
        self.mb.AddProg('shelxe')
      if not self.ph.inp.Get('model', typ='substr', has_atomtypes=True) and \
        self.ph.GetParam('target',capitalized=True) in ('SAD','SIRAS'):
          common.Error('No substructure model or atomtype(s) inputted for {0}.'.format(self.name))
    else:
      # the arpwarp specific part shall follow
      arp=self.GetProg('arpwarp')
      arp_auto=self.GetProg('arpwarp_auto')
      if arp and not arp_auto and not arp.inp.Get('datafile',typ='par'):
        arp_auto=self.AddProg('arpwarp_auto',ind=0)
        arp_auto.inp=copy.copy(arp.inp)
        for key in arp.key_list:
          if arp.IsKey(key):
            arp_auto.SetArg(key,arp.GetKey(key))
        arp_auto.SetArg('parfile','arp.par')
      # catch output continously (eg to stop early if certain stat is reached) - shelxd only as of now.
      if self.GetParam('catch_output') is False:
        self.Info('Output catching from {0} disabled, continous report wont be available.')
      else:
        arp.interact_output=True
    process.TreatInOutPar(self,set_all_par)

  def save_reset_Win(self):
    if os.name=='nt':  # windows cannot spawn instances with non-picklable attributes
      if hasattr(self,'ccp4i2job'):  i2job,self.ccp4i2job=self.ccp4i2job,None
      else:  i2job=None
      logfh,parent=self.logfilehandle,self.parent_process
      self.logfilehandle,self.parent_process=None,None
      return i2job,logfh,parent
    else:
      return None,None,None

  def restore_Win(self,data):
    i2job,logfh,parent=data[0],data[1],data[2]
    if os.name=='nt':  # reattaching the atrributes for windows
      if hasattr(self,'ccp4i2job'):  self.ccp4i2job=i2job
      self.logfilehandle,self.parent_process=logfh,parent

  def RunBody(self,*args,**kwargs):
    if not self.GetParam('handdet'):
      self.RunMBRef()
    else:
      windata=self.save_reset_Win() # windows cannot spawn instances with non-picklable attributes
      manager = multiprocessing.Manager()
      queue, queue2 = manager.Queue(), manager.Queue()
      mbref_hand,self.mbr_hand=[],[]
      for i in range(2):
        mbref_hand.append( multiprocessing.Process(target=self.RunMBRefHands, args=(i,queue,queue2)) )
        mbref_hand[i].daemon=True
        mbref_hand[i].start()
        self.mbr_hand.append(None)
      self.restore_Win(windata)
      num_finished, stop_other = 0, None
      while num_finished<len(mbref_hand):
        j,mbr_now,logdict = queue.get()
        if j<0:
          num_finished+=1
          if num_finished<len(mbref_hand):
            common.Warning(mbr_now)
          else:
            common.Error(mbr_now)
        else:
          self.ph, self.mb, self.dm = mbr_now.ph, mbr_now.mb, mbr_now.dm
          self.ref_res_all, self.mb_res_all = mbr_now.ref_res_all, mbr_now.mb_res_all
          for o in mbr_now.out.GetAll(stored_order=True):  self.out.Add(o)
          self.PrintActualLogGraph(hand=j+1,**logdict)
          # stop the other hand when this hand succeeds
          if mbr_now.cyc_fin and stop_other is None:
            stop_other = int(not j)
          if stop_other is j:
            queue2.put((j,'stop'))
          else:
            queue2.put((j,'continue'))
        if 'update' in logdict and logdict['update']=='end':
          num_finished+=1
          self.mbr_hand[j]=mbr_now
          self.mbr_hand[j].report_R,self.mbr_hand[j].report_Rfree=self.report_R,self.report_Rfree
          if stop_other is j:
            self.stop_message = mbr_now.stop_message
            self.PrintActualLogGraph(update='end')
            self.stop_message = None
      chosen_hand,chosen_mbr = 0,self.mbr_hand[0]
      if self.mbr_hand[0] is None or (self.mbr_hand[1] and self.mbr_hand[1].report_R<self.mbr_hand[0].report_R):
        chosen_hand,chosen_mbr = 1,self.mbr_hand[1]
      self.report_R,self.report_Rfree = chosen_mbr.report_R,chosen_mbr.report_Rfree
      for o in chosen_mbr.out.GetAll(stored_order=True):  self.out.Add(o)
      if self.mbr_hand[not chosen_hand]:
        self.out2=inout.input_output(is_output=True,parent=self)
        self.out2=self.mbr_hand[not chosen_hand].out
      self.ph, self.mb, self.dm = chosen_mbr.ph, chosen_mbr.mb, chosen_mbr.dm
      self.ref_res_all, self.mb_res_all = chosen_mbr.ref_res_all, chosen_mbr.mb_res_all
      self.stop_message,self.result_str = chosen_mbr.stop_message, chosen_mbr.result_str
      self.spacegroup_change=False
      if chosen_mbr.out.Get('mapcoef').GetSpacegroupNumber(self)!=self.inp.Get('fsigf').GetSpacegroupNumber(self):
        self.spacegroup_change=chosen_mbr.out.Get('mapcoef').GetSpacegroup(self)
        self.result_str+='  Spacegroup changed from {0} to {1}.'.format(self.inp.Get('fsigf').GetSpacegroup(self),self.spacegroup_change)
      self.PrintActualLogGraph(update='end')  # needed here since the messages are only now in the main thread

  def RunMBRefHands(self,num_proc=-1,queue=None,queue2=None):
    self.num_proc=num_proc
    self.queue,self.queue2=queue,queue2
    if self.logfilehandle:
      self.logfilehandle.close()
    try:
      if num_proc>0:
        self.SetRunDir(self.rundir+'_hand'+str(num_proc+1), reset_subtree=True)
        self.ReplaceInpOut('inp','inp2')
      self.TreatInOutPar()
      self.RunPreprocess()
      self.RunMBRef()
    except Exception as e:
      error='Parallel building process {0} failed with error: {1}'.format(num_proc, str(e))
      queue.put((-1,error,'err'))

  def RunMBRef(self):
    if self.processes:
      # if shelxe should be used then start with it
      if self.mb.GetProg('shelxe') and not self.mb.GetProg(supported=True).nick=='shelxe':
        self.mb.AddProgCopy( self.mb.programs.pop(self.mb.programs.index(self.mb.GetProg('shelxe'))), ind=0, deeper_copy=True )
      for self.mb.cyc in range(1, self.GetParam('bigcyc')+1):
        self.out.ClearAll(propagate=False)
        # build model
        if self.ph.out.mapcoef:
          self.mb.inp.Set(self.ph.out.mapcoef)
          if self.dm and (not self.cyc_fin or self.ph.R>self.R_thres-0.05) and self.mb.GetProg(supported=True).nick=='buccaneer':
            self.dm.GetProcess('dm').inp.Add(self.ph.out.Get('mapcoef',typ='combined'))
            ###!!!!!!
            if self.cyc_fin:  self.dm.SetParam('dmcyc',4)
            self.dm.Run()
            cphcomb=self.AddProg('cphasecombine', propagate_out=False)
            cphcomb.inp.Set([self.dm.out.Get('mapcoef',typ='combined'),]+[self.ph.out.Get('mapcoef',typ='combined')])
            if self.cyc_fin: # the order is actually reverted
              cphcomb.SetKey('weight-hl-1',0.45,keep_previous=False), cphcomb.SetKey('weight-hl-2',0.55,keep_previous=False)
            else:
              cphcomb.SetKey('weight-hl-1',0.55,keep_previous=False), cphcomb.SetKey('weight-hl-2',0.45,keep_previous=False)
            cphcomb.Run()
            self.mb.inp.Set(cphcomb.out.mapcoef)
            self.programs.remove(cphcomb)
        # skipping correlation mode if number of cycles is large enough
        if self.mb.cyc%3==0 and not self.cyc_fin and self.GetParam('bigcyc')-self.mb.cyc>3 and self.mb.GetProg('buccaneer'):
        #!!!!!!!
        #if self.mb.cyc%3==0 and self.GetParam('bigcyc')-self.mb.cyc>3 and self.mb.GetProg('buccaneer'):
          self.mb.GetProg('buccaneer').SetKey('fast',False,keep_previous=False)
          self.mb.inp.Clear('model')
          #!!!!!!!
          #self.ph.inp.Delete( self.ph.inp.Get('mapcoef',typ='best',custom='mlhl') )
        else:
          self.mb.GetProg('buccaneer').SetKey('fast',True,keep_previous=False)
        if not self.mb.GetProg():
          common.Error('No program specified for {0}'.format(self.mb.name))
        #self.mb.GetProg().runname=self.mb.GetProg().name+'_cyc'+str(self.mb.cyc)
        self.mb.GetProg(supported=True).SetRunDir(os.path.join(self.mb.rundir,self.mb.GetProg(supported=True).nick+'_cyc'+str(self.mb.cyc)))
        self.mb.Run()
        try:
          self.UpdateResults('mb')
        except common.Unsuccessful:
          if self.mb.GetProg('buccaneer'):
            self.mb.GetProg('buccaneer').SetKey('model-index',3,keep_previous=False)
          elif self.mb.GetProg('shelxe'):
            self.mb_parallel[i].GetProg().SetArg('G',0.6, keep_previous=False)
          self.mb_res_all.pop()
          self.mb.Run()
          try:
            self.UpdateResults('mb')
          except common.Unsuccessful:
            if self.mb.GetProg('buccaneer'):
              self.mb.GetProg('buccaneer').SetKey('model-index',5,keep_previous=False)
            elif self.mb.GetProg('shelxe'):
              self.mb_parallel[i].GetProg().SetArg('G',0.55, keep_previous=False)
            self.mb_res_all.pop()
            self.mb.Run()
            self.UpdateResults('mb')
        # run refinement
        if self.ph.out.Get('model',typ='substr',filetype='pdb'):
          self.ph.inp.Set(self.ph.out.Get('model',typ='substr',filetype='pdb'))
        self.ph.inp.Add(self.mb.out.Get('model'))
        self.ph.GetProg(supported=True).SetRunDir(os.path.join(self.ph.rundir,self.ph.GetProg(supported=True).nick+'_cyc'+str(self.mb.cyc)))
        self.ph.Run()
        self.UpdateResults('ph')
        # adding hydrogens if R is good (and not at low res)
        if self.hydrog_no_set and self.ref_res_all and self.ref_res_all[-1][0]<0.4:
          self.ph.GetProg('refmac').SetKey('make',False,keep_previous=False)
        # switch to buccaneer from shelxe 
        if self.mb.GetProg('buccaneer') and not self.mb.GetProg(supported=True).nick=='buccaneer':
          if self.mb.cyc>=min(3,self.GetParam('bigcyc')/2) or self.ph.R<self.R_thres+0.02:
            self.mb.AddProgCopy( self.mb.programs.pop(self.mb.programs.index(self.mb.GetProg('buccaneer'))), ind=0, deeper_copy=True )
            self.mb.GetProg('buccaneer').inp.Add( self.out.Get('model',typ='partial') )
            self.Info('Switching from ShelxE building to Buccaneer building.')
        # separate substructure/partial if needed
        if self.ph.out.Get('model',typ='partial+substr') and not self.ph.out.Get('model',typ='substr'):
          separate_models=self.ph.GetOrAddProcess('sepsubstrprot')
          separate_models.inp.Set(self.ph.out.Get('model',typ='partial+substr'))
          separate_models.Run()
        if (hasattr(self,'ccp4i2job') and self.ccp4i2job.testForInterrupt()) or os.path.isfile(self.stop_file):
          self.stop_message='Stopping early on user request!'
        # stopping if built
        if self.cyc_fin and self.mb.cyc>=self.cyc_fin:
          self.stop_message='Majority of the model should be successfully built.'
        if self.stop_message:
          self.Info(self.stop_message)
          break
        # set correlation mode for Buccaneer etc
        if self.mb.GetProg('buccaneer'):
          #!!!!!!!!
          #if not self.ph.inp.Get('mapcoef',typ='best',custom='mlhl'):
          #  self.ph.inp.AddCopy( self.dm.inp.Get('mapcoef',typ='best',custom='mlhl') )
          self.mb.GetProg('buccaneer').inp.Add( self.out.Get('model',typ='partial') )
          self.mb.GetProg('buccaneer').SetKey('fast',True,keep_previous=False)
          if not self.mb.GetProg('buccaneer').IsKey('cycles'):
            self.mb.GetProg('buccaneer').SetKey('cycles',2)
      self.UpdateResults('end')
    else:
      # arpwarp (or possibly other mb+ref programs in the future)
      arp_auto=self.GetProg('arpwarp_auto')
      arp=self.GetProg('arpwarp')
      if arp_auto:
        if not self.inp.Get('mapcoef',typ=('combined','best'),col=('ph','fom')) and self.inp.Get('mapcoef',typ=('combined','best'),col='hla'):
          fom2hl=self.AddProg('chltofom',propagate_out=False)
          fom2hl.Run()
          self.programs.remove(fom2hl)
        arp_auto.Run()
        if arp:
          arp.inp.Add(arp_auto.out.Get('datafile',typ='par'))
          arp.fsf=arp_auto.fsf
          arp.SetRunDir(arp_auto.rundir)
          for line in fileinput.input(arp.inp.Get(typ='par').GetFileName(), inplace=True):
            if line.startswith('set JOB_ID'):
              line='set JOB_ID = arp\n'
            sys.stdout.write(line)
      if arp:
        if self.GetParam('target') == 'SAD' and not self.GetParam('arp_sad'):
          parf=arp.inp.Get(typ='par').GetFileName()
          obj,N=self.inp.GetTargetObjects(self.GetParam('target'), modtyp='pdb', fsftyp='mtz')
          if obj is None:
            common.Error('Input data/model for ARP/wARP SAD could not be retrived.')
          with open(parf,'a') as f:
            f.write('set phaseref = SADH\n')
            f.write("set F_1 = '{0}'\n".format(obj['f+'].GetLabel('f')))
            f.write("set SIGF_1 = '{0}'\n".format(obj['f+'].GetLabel('sigf')))
            f.write("set F_2 = '{0}'\n".format(obj['f-'].GetLabel('f')))
            f.write("set SIGF_2 = '{0}'\n".format(obj['f-'].GetLabel('sigf')))
            for at in obj['mod'].GetAtomTypes(getlist=True):
              fp,fpp,dn,att=obj['mod'].Getfpfpp(at,obj['f+'].dname)
              if fp is not None and fpp is not None:
                f.write("set ANOM = ' FORM {0} {1} {2}'\n".format(obj['mod'].GetAtomType(),fp,fpp))
            f.write("set HEAVYIN = '{0}'\n".format(obj['mod'].GetFileName('pdb')))
            ### this should be disabled and the various heavy atom exceptions may be removed from dummy/refmac
            #f.write("set heavyin = '{0}'\n".format(obj['mod'].GetFileName('pdb')))
          dummy_env = os.path.join( os.getenv('CRANK2'), 'bin', 'dummy' ) if os.getenv('CRANK2') else None
          if not dummy_env or not os.path.isdir(dummy_env):
            dummy_env = os.path.join( os.getenv('CCP4'), 'share', 'ccp4i', 'crank', 'bin', 'dummy' )
          inp_parf=os.path.join(os.path.dirname(parf),'input.par')
          shutil.copy(parf, inp_parf)
          shutil.copy(obj['mod'].GetFileName('pdb'), obj['mod'].GetFileName('pdb')+'.ref')
          arp.env['PATH'] = "{0}{2}{1}".format(dummy_env, os.getenv('PATH'), os.pathsep)
          arp.env['CCP4I_TOP'] = dummy_env
          arp.env['CBIN'] = dummy_env
        arp.Run()
        self.UpdateResults('end')
          #with open('arp_warp.log','w') as g:
          #  subprocess.call( [os.path.join(os.getenv('warpbin'),'warp_tracing.sh'),inp_parf,'1'], stdout=g )


  def UpdateInfo(self,line):
    # this is arpwarp specific; just pushes results from arpwarp log through UpdateResults()
    prog=self.GetProg(supported=True)
    if not prog or prog.nick!='arpwarp':
      return
    self.UpdateResults('arpwarp',line)

  def UpdateResults(self, update, line=None):
    if update=='arpwarp':
      mbref_stat, update = self.GetProg(supported=True).GetStatGrep, ''
      if mbref_stat('build_cycle',from_str=line):
        self.cyc=mbref_stat('build_cycle',from_str=line)
        self.act_mb_res = []
      elif mbref_stat('rfact', from_str=line):
        self.ref_res_all.append( (mbref_stat('rfact',from_str=line), 0.0, mbref_stat('rfree',from_str=line)) )
        self.report_R = self.ref_res_all[-1][0]
        update='ph'
      elif mbref_stat('res_built',from_str=line):
        self.act_mb_res.append( (mbref_stat('res_built',from_str=line), mbref_stat('frag_built',from_str=line)) )
      elif mbref_stat('round_cycle',from_str=line):
        res_frac = self.act_mb_res[mbref_stat('round_cycle',from_str=line)-1]
        self.mb_res_all.append( (self.cyc, res_frac[0], res_frac[1], None,None) )
        update='mb'
    elif update=='mb':
      mb_stat = self.mb.GetProg(supported=True).GetStat
      if self.mb.GetProg(supported=True).nick=='shelxe':
        self.mb_res_all.append( (self.mb.cyc, mb_stat('res_built'), len(list(filter(None,mb_stat('res_built_per_frag')))), 0, 0) )
      else:
        self.mb_res_all.append( (self.mb.cyc, mb_stat('res_built'), mb_stat('frag_built'), 
                              mb_stat('compl_chain'), mb_stat('compl_res')) )
      if not self.mb_res_all[-1][1]:
        self.result_str = "No residues could be traced, map too noisy. Structure solution was unsuccessful."
        common.Error(self.result_str, nosuccess=True)
    elif update=='ph':
      ph_stat = self.ph.GetProg(supported=True).GetStat
      self.ref_res_all.append( (ph_stat('rfact')[-1], ph_stat('fom')[-1], (ph_stat('rfree',accept_none=True) or [None])[-1]) )
      self.report_R, self.report_Rfree = self.ref_res_all[-1][0], self.ref_res_all[-1][2]
      if self.cyc_fin is None and self.report_R<self.R_thres and self.GetParam('minbigcyc'):
        self.cyc_fin=self.GetParam('minbigcyc')+self.mb.cyc+min(self.mb.cyc,self.GetParam('minbigcyc'))-2
    elif update=='end' and self.processes:
      # disabled for arpwarp - R factor not a good criterion there
      if self.ref_res_all[-1][0]<0.4:
        self.result_str = "Majority of the model should be successfully built!"
      elif self.ref_res_all[-1][0]<0.5 or (self.GetCrankParent() and self.GetCrankParent().GetProcess('phdmmb') and self.GetCrankParent().GetProcess('phdmmb').best_cc>0.25):
        self.result_str = "Partial model was probably correctly built."
      else:
        self.result_str = "The structure solution was unsuccessful."
      #self.result_str = "Final model contains {0} residues and provides R-factor of {1}.".format(
      #                    self.mb_res_all[-1][1], self.ref_res_all[-1][0])
    if update:
      self.PrintActualLogGraph(update=update)

  def PrintActualLogGraph(self, update=None, paral_num=None, hand=-1):
    if self.queue and self.num_proc>=0:
      # in case of both hands, passing to/from the main mbref process
      queue,queue2,parent,lfh,phlfh,mblfh,dmlfh = self.queue,self.queue2,self.parent_process,self.logfilehandle,self.ph.logfilehandle,self.mb.logfilehandle,self.GetProcess('dmfull').logfilehandle if self.GetProcess('dmfull') else None
      # the queue and parent process are removed and reattached after putting to the queue
      self.queue,self.queue2,self.parent_process,self.logfilehandle,self.mb.logfilehandle,self.ph.logfilehandle = None,None,None,None,None,None
      if dmlfh: self.GetProcess('dmfull').logfilehandle = None
      if hasattr(self,'ccp4i2job'):  i2job,self.ccp4i2job=self.ccp4i2job,None
      queue.put((self.num_proc,self,{'update':update,'paral_num':paral_num}))
      self.queue,self.queue2,self.parent_process,self.logfilehandle,self.ph.logfilehandle,self.mb.logfilehandle = queue,queue2,parent,lfh,phlfh,mblfh
      if dmlfh: self.GetProcess('dmfull').logfilehandle = dmlfh
      if hasattr(self,'ccp4i2job'):  self.ccp4i2job=i2job
      numproc,check_stop=self.queue2.get()
      if check_stop=='stop' and self.num_proc==numproc and not self.stop_message:
        self.SetParam('minbigcyc',0)
        self.stop_message = "Hand {} building stopped since the other hand succeeded.".format(numproc+1)
        self.Info(self.stop_message)
      return
    if self.opened_loggraph:
      self.GetLogGraphHandle().seek(0,0)
      self.LGInfo(self.GetCCP4Header())
      self.LGInfo(self.GetLogGraphPrograms())
      if self.ref_res_all:
        self.LGInfo('\n $TABLE : R-factor and FOM per cycle:')
        self.LGInfo('$GRAPHS    : R-factor and FOM per cycle :A:1,2,3:\n$$')
        self.LGInfo('Cycle  R-factor   FOM  $$ $$')
        for cyc,(R,fom,Rfree) in enumerate(self.ref_res_all):
          self.LGInfo('{0:3}  {1:8.3f}  {2:8.3f}'.format(cyc+1,R,fom))
        self.LGInfo('$$\n')
      if self.mb_res_all:
        self.LGInfo('\n $TABLE : Built per building cycle:')
        self.LGInfo('$GRAPHS    : Residues and fragments per cycle :A:1,2,3:\n$$')
        self.LGInfo('Cycle Residues Fragments  $$ $$')
        for cyc,res,fr,c1,c2 in self.mb_res_all:
          self.LGInfo('{0:3}  {1:8}  {2:8}'.format(cyc,res,fr))
        self.LGInfo('$$\n')
      if update=='end':
        self.LGInfo('\n$TEXT:Result: $$ Final result $$')
        self.LGInfo("{0}\n$$\n".format(self.result_str))
    if self.rv_report is not None:
      if not hasattr(self,'rv_plotmb'):
        handstr = '' if hand<0 else ' - hand 1'
        self.rv_plotmb = self.rv_report.Plot( 'Residues and fragments vs build cycle'+handstr, "Cycle", \
              "Number of residues and fragments", block="Model building statistics", legendloc='nw', ymin=0., intx=True )
        self.rv_res = self.rv_plotmb.PlotLine(["x"],["residues"])
        self.rv_fr = self.rv_plotmb.PlotLine(["x"],["fragments"])
        self.rv_fr.CustomTick( "x", 0, "0" )
        if hand>0:
          handstr = ' - hand 2'
          self.rv_plotmb2 = self.rv_plotmb.parent.Plot( 'Residues and fragments vs build cycle'+handstr, "Cycle", \
              "Number of residues and fragments", legendloc='nw', ymin=0., intx=True )
          self.rv_res2 = self.rv_plotmb2.PlotLine(["x"],["residues"])
          self.rv_fr2 = self.rv_plotmb2.PlotLine(["x"],["fragments"])
          self.rv_fr2.CustomTick( "x", 0, "0" )
      if self.mb_res_all and update=='mb':
        if hand<=1:
          rv_res, rv_fr = self.rv_res, self.rv_fr
        else:
          rv_res, rv_fr = self.rv_res2, self.rv_fr2
        if len(self.ref_res_all)==10:
          rv_fr.Reset(reset_data=False)
        cyc,res,fr,c1,c2 = self.mb_res_all[-1]
        cst = len(self.mb_res_all) if len(self.mb_res_all)<10 else False
        rv_res.Data( cyc, res, int )
        rv_fr.Data( cyc, fr, int, custom_x_tick=cst, flush=True )
      if self.ref_res_all and update=='ph':
        # this displays the actual pdb+mtz updated cycle by cycle
        # the files are continuously rewritten which might cause issues - we'll need to test!
        out_mtz = None
        out_pdb = self.out.Get('model',typ=('partial','partial+substr')).GetFileName()
        if out_pdb and os.path.isfile(out_pdb): # skip for arpwarp where the pdb only exists at the end of the job
          out_mtz = crvapi.SetMtzFWT(self,self.out.Get('mapcoef',typ='combined'))
          #out_map = self.out.Get('mapcoef',typ='weighted',filetype='map').GetFileName('map')
          #out_dmap = self.out.Get('mapcoef',typ='diff',filetype='map',conv_opts=['outfilename','diff.map']).GetFileName('map')
          rv_dir=os.path.join(self.rundir,'rvapi_display')
          if not os.path.isdir(rv_dir):
            os.mkdir(rv_dir)
          shutil.copy(out_pdb,rv_dir),shutil.copy(out_mtz,rv_dir)#,shutil.copy(out_dmap,rv_dir),shutil.copy(out_map,rv_dir)
          out_pdb = os.path.join(rv_dir,os.path.basename(out_pdb)) #shutils.copy only returns the file from py3.3
          out_mtz = os.path.join(rv_dir,os.path.basename(out_mtz))
          #out_map = os.path.join(rv_dir,os.path.basename(out_map))
          #out_dmap = os.path.join(rv_dir,os.path.basename(out_dmap))
        if not hasattr(self,'rv_plotref'):
          handstr = '' if hand<0 else ' - hand 1'
          self.rv_plotref = self.rv_plotmb.parent.parent.Plot( 'R-factor and FOM vs refinement cycle'+handstr, \
                  "Cycle", "R-factor and FOM", block="Refinement statistics", legendloc='nw', ymin=0., intx=True )
          self.rv_fom = self.rv_plotref.PlotLine(["x"],["FOM"])
          self.rv_R = self.rv_plotref.PlotLine(["x"],["R"])
          self.rv_R.CustomTick( "x", 0, "0" )
          self.rv_Rfree = self.rv_plotref.PlotLine(["x"],["Rfree"], color='darkorchid')
          if hand>0:
            handstr = ' - hand 2'
            self.rv_plotref2 = self.rv_plotref.parent.Plot( 'R-factor and FOM vs refinement cycle'+handstr, \
                  "Cycle", "R-factor and FOM", legendloc='nw', ymin=0., intx=True )
            self.rv_fom2 = self.rv_plotref2.PlotLine(["x"],["FOM"])
            self.rv_R2 = self.rv_plotref2.PlotLine(["x"],["R"])
            self.rv_R2.CustomTick( "x", 0, "0" )
            self.rv_Rfree2 = self.rv_plotref2.PlotLine(["x"],["Rfree"], color='darkorchid')
          if (not self.parent_process or not self.parent_process.ccp4i2) and self.out.Get('mapcoef',typ='combined'):
            self.rv_report.Text("<BR><i><small>Output (updating after each refinement cycle):</small></i>")
            #self.rv_files = self.rv_report.DataFile(out_pdb, "xyz", "<small>Built model and map</small>")
            self.rv_files = self.rv_report.DataFile("", "", "<small>Built model and map - hand 1</small>")
            if hand>0:
              self.rv_files2 = self.rv_report.DataFile("", "", "<small>Built model and map - hand 2</small>")
        handstr='' if hand<=1 else '2'
        if (not self.parent_process or not self.parent_process.ccp4i2) and not hasattr(self,'rv_files_created'+handstr) and out_mtz:
          setattr(self,'rv_files_created'+handstr, getattr(self,'rv_files'+handstr).DataFile(out_pdb, "xyz"))
          getattr(self,'rv_files'+handstr).DataFile(out_mtz, "hkl:map")
          #getattr(self,'rv_files'+handstr).DataFile(out_map, "hkl:ccp4_map")
          #getattr(self,'rv_files'+handstr).DataFile(out_dmap, "hkl:ccp4_dmap")
        if hand<=1:
          rv_fom, rv_R, rv_Rfree = self.rv_fom, self.rv_R, self.rv_Rfree
        else:
          rv_fom, rv_R, rv_Rfree = self.rv_fom2, self.rv_R2, self.rv_Rfree2
        if len(self.ref_res_all)==10:
          rv_R.Reset(reset_data=False)
        R,fom,Rfree = self.ref_res_all[-1]
        cst = len(self.ref_res_all) if len(self.ref_res_all)<10 else False
        if fom:
          rv_fom.Data( len(self.ref_res_all), fom, int )
        if Rfree:
          rv_Rfree.Data( len(self.ref_res_all), Rfree, int, custom_x_tick=cst )
        rv_R.Data( len(self.ref_res_all), R, int, custom_x_tick=cst )
        # flushing only after both hands outputs are created
        if hasattr(self,'rv_files_created') and (hasattr(self,'rv_files_created2') or hand<0):
          crvapi.Flush()
      if update=='end':
        self.report_R, self.report_Rfree = self.ref_res_all[-1][0], self.ref_res_all[-1][2]
        if self.stop_message:
          self.rv_report.Text("<i><small>Stop:</small></i>")
          self.rv_report.Text('&emsp;'+self.stop_message)
        if self.result_str:
          self.rv_report.Text("<i><small>Result:</small></i>")
          self.rv_report.Text('&emsp;'+self.result_str, flush=True)

