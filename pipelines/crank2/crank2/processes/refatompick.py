#!/usr/bin/python
import os,sys
from process import process
import common
par=common.parameter


class refatompick(process):
  name="iterative substructure improvement and phasing"
  short_name="substructure improvement"
  supported_procs=["ref", "phas", "atomsfrommap"]
  supported_params = {}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True, share=True )
  supported_params['refcyc'] = par( desc='Number of refinement cycles in each iteration', typ=int )
  supported_params['num_iter'] = par( desc='Number of refinement+picking iterations', typ=int )
  supported_params['rms_threshold'] = par( desc='RMS threshold for atom picking from map (negative disables picking)', typ=float, share=True )
  supported_params['occ_cut'] = par( desc='Occupancy threshold for atom deletion (relative to largest occ; negative to disable)', typ=float )
  supported_params['skip_init_ref'] = par( desc='Skip the initial refinement', typ=bool )
  supported_params['skip_last_pick'] = par( desc='Skip atom picking in the last iteration', typ=bool )
  supported_params['skip_last_occ_cut'] = par( desc='Skip atom deletion in the last iteration', typ=bool )
  supported_params['skip_init_occ_cut'] = par( desc='Skip atom deletion in the initial model', typ=bool )
  supported_params['max_new_atoms'] = common.parameter( desc='Maximal number of atoms added in one iteration', typ=(int,bool), share=True )
  supported_params['bfactor'] = par( desc='Set B factor of the new atoms to this value', typ=(float,bool), share=True )
  supported_params['res_cut'] = par( desc='Use (automatic if True) resolution cutoff in refinement', typ=(float,bool), share=True )
  supported_params['res_cut_cyc'] = par( desc='Number of refinement cycles with the resolution cutoff (default all)', typ=int, share=True )

  def Init(self):
    self.lastcut=None
    self.noimpr=None

  def TreatInOutPar(self, set_all_par=False):
    # input
    self.ph = self.GetProcess('phas')
    if not self.ph:
      self.ph = self.GetProcess('ref')
    if not self.ph:
      if self.inp.Get('model',typ=('partial','partial+substr'),filetype='pdb'):
        self.ph = self.AddProcess('ref')
      else:
        self.ph = self.AddProcess('phas')
    if not self.ph.GetProg():
      self.ph.AddProg('refmac')
    self.ph.TreatInOutPar()
    self.afm = self.GetOrAddProcess('atomsfrommap')
    self.afm.TreatInOutPar()
    # parameters
    if not self.IsParam('num_iter'):
      self.SetParam('num_iter', 3)
    if self.GetParam('rms_threshold') is None:
      self.SetParam('rms_threshold', 4.75)
    if self.GetParam('skip_last_pick') is None:
      self.SetParam('skip_last_pick')
    if self.GetParam('occ_cut') is None:
      if self.ph.nick=='phas':
        self.SetParam('occ_cut',0.3)
        if self.inp.Get('model',has_atomtypes=True) and self.inp.Get('model',has_atomtypes=True).GetAtomType() in ('S','BR'):
          self.SetParam('occ_cut',0.25)
        elif self.inp.Get('model',has_atomtypes=True) and self.inp.Get('model',has_atomtypes=True).GetAtomType() in ('I',):
          self.SetParam('occ_cut',0.2)
      else:
        self.SetParam('occ_cut',0.15)
    if self.GetVirtPar('refcyc') is None and self.ph.nick=='ref':
      self.SetVirtPar('refcyc', 10)
      if self.ph.GetParam('low_res_par'):
        self.SetVirtPar('refcyc', 20)
    #if self.GetVirtPar('refcyc') is None and self.ph.nick=='phas':
    #  self.SetVirtPar('refcyc', 25)
    if self.IsVirtPar('refcyc') and not self.ph.IsVirtPar('cycles'):
      self.ph.SetParam('cycles', self.GetVirtPar('refcyc'))
#    if not self.IsParam('res_cut') and self.GetParam('target')=='SAD' and self.ph.nick=='phas' and \
#       self.inp.Get('model',has_atomtypes=True) and self.inp.Get('model',has_atomtypes=True).GetAtomType() in ('S','CA','CR','MN','FE','CO','NI','CU','ZN','SE','BR'):
       #self.inp.Get('model',has_atomtypes=True) and self.inp.Get('model',has_atomtypes=True).GetAtomType() in ('S','CA'):
#      self.SetParam('res_cut')
    process.TreatInOutPar(self,set_all_par)

  def GetOccB(self,pdb,occ=1):
    if not pdb:
      return None,None
    pdbcur=self.GetOrAddProg('pdbcur')
    pdbcur.inp.Set(pdb)
    B_or_occ='occup' if occ else 'B'
    with open(pdb.GetFileName('pdb')) as f:
      occ = pdbcur.GetStatGrep(B_or_occ,from_str=f.read())
    return occ,pdbcur

  def CutOcc(self,occ,pdbcur,init=False):
    if not occ:
      return None
    maxocc = max(occ) if occ else 1.
    pdbcur.SetKey('cutocc', self.GetParam('occ_cut')*maxocc, keep_previous=False)
    pdbcur.Run()
    message=pdbcur.GetStat('atoms_deleted_occ',accept_none=True)
    if message:
      self.Info(message)
      if self.rv_report is not None:
        self.rv_report.Text('&emsp;'+message if not init else message)
    self.inp.Add(pdbcur.out.Get('model'))
    return message

  def RunBody(self,*args,**kwargs):
    if self.GetParam('res_cut') and not self.ph.GetProg(supported=True).IsKey('reso'):
      if self.GetParam('res_cut') is True:
        substrdet=self.AddProcess('substrdet', propagate_out=False)
        cutoff = substrdet.EstimateCutOff()
        self.processes.remove(substrdet)
      else:
        cutoff = float( self.GetParam('res_cut') )
      if cutoff:
        self.ph.GetProg(supported=True).SetKey('reso', cutoff)
    model_type = 'substr'  if self.ph.nick=='phas'  else 'partial+substr'
    if self.ph.nick=='phas' and not self.GetParam('skip_init_occ_cut'):
      occ,pdbcur=self.GetOccB( self.inp.Get('model',typ=model_type,filetype='pdb') )
      self.CutOcc( occ,pdbcur, init=True )
    refprob=0
    for i in range(self.GetParam('num_iter')):
      self.LGInfo('\n$TEXT: Iteration {0}: $$ Atom picking result $$'.format(i+1))
      self.lastcut=None
      if self.rv_report is not None:
        self.rv_report.Text("<i><small>Iteration {0}.</small></i>".format(i+1))
        self.afm.rv_report = self.rv_report
      #if self.GetParam('res_cut') and i==self.GetParam('num_iter')-1:
      #  self.ph.GetProg().SetKey('reso',False,keep_previous=False)
      if i>0 or not self.GetParam('skip_init_ref'):
        self.ph.GetProg().runname=self.ph.GetProg().name+'_cyc'+str(i)
        self.ph.GetProg().outfilename['pdb']=self.ph.GetProg().name+'_cyc'+str(i)+'.pdb'
        if self.GetParam('res_cut_cyc')==i:
          self.ph.GetProg().SetKey('reso',False,keep_previous=False)
        # try bootstrapping from no substr. atoms inputted
        if self.ph.nick=='ref' and not self.ph.inp.Get('model',('partial+substr','substr'),filetype='pdb') and self.ph.inp.Get('model','partial',filetype='pdb'):
          ph_backup=self.ph.GetProg().BackupAnyPars()
          self.ph.GetProg().SetKey('anom','mapo')
          self.ph.GetProg().SetKey('ncyc',0)
          self.ph.SetParam('no_substr',True)
          self.ph.Run()
          self.ph.GetProg().RestoreAnyPars(ph_backup)
          self.ph.SetParam('no_substr',False)
          self.afm.inp.AddCopy(self.out.Get('model',typ='partial'))
        else:
          self.ph.Run()
          # to cover the obscure but possible case of substr only with ref
          if model_type=='partial+substr' and not self.ph.out.Get('model',typ='partial+substr') and self.ph.out.Get('model',typ='substr'):
            model_type='substr'
          if self.GetParam('occ_cut')>=0.0 and (i<self.GetParam('num_iter')-1 or not self.GetParam('skip_last_occ_cut')):
            occ,pdbcur=self.GetOccB( self.out.Get('model',typ=model_type,filetype='pdb') )
            max_occ=max(occ) if occ else 1.0
            B,pdbcur=self.GetOccB( self.out.Get('model',typ=model_type,filetype='pdb'), occ=0 )
            aver_B=sum(B)/len(B) if B else 500.0
            if ((occ and occ[0]<=0.7*max_occ) or (B and B[0]>=1.5*aver_B) or (len([o for o in occ if o<0.15*max_occ])>0.7*len(occ))): #or len(self.ph.GetProg().GetStat('rfact'))>=50:
              cyc=self.ph.GetParam('cycles')
              problemcyc = 0  if refprob  else 1
              self.ph.SetParam('cycles',problemcyc)
              self.ph.Run()
              self.ph.SetParam('cycles',cyc)
              occ,pdbcur=self.GetOccB( self.out.Get('model',typ=model_type,filetype='pdb') )
              self.CutOcc( occ,pdbcur )
              if refprob:
                self.Info('Substructure improvement not likely to improve the substructure. Stopping.')
                self.noimpr=True
                break
              refprob+=1
            self.lastcut=self.CutOcc( occ,pdbcur )
          model_out=self.out.Get('model',typ=model_type)
          if model_out:
            self.afm.inp.AddCopy(model_out)
        self.afm.inp.Set(self.ph.out.GetAll('mapcoef'))
      if self.GetParam('rms_threshold')>0.0 and (i<self.GetParam('num_iter')-1 or not self.GetParam('skip_last_pick')):
        self.afm.GetProg().runname=self.afm.GetProg().name+'_cyc'+str(i)
        self.afm.Run()
        if (not self.afm.GetProg().out.model or self.afm.GetProg().added_atoms==0) and i>=self.GetParam('num_iter')-2:  # run only short last cycle if no atoms added (useful eg for 5b82)
          self.ph.SetParam('cycles',1)
      if self.out.Get('model'):
        self.ph.inp.AddCopy(self.out.Get('model'))
      else:
        common.Error('No substructure inputted and no substructure atoms found from anomalous maps.',nosuccess=True)
      self.LGInfo('\n$$')
    # quick phasing with all reflections if there was resolution cutoff
    if self.GetParam('res_cut') and (not self.GetParam('res_cut_cyc') or self.GetParam('res_cut_cyc')==self.GetParam('num_iter')):
      self.ph.GetProg().SetKey('reso',False,keep_previous=False)
      self.ph.GetProg().runname=self.ph.GetProg().name+'_phas'
      self.ph.GetProg().outfilename['pdb']=self.ph.GetProg().name+'.pdb'
      self.ph.SetParam('cycles',3)
      self.ph.Run()
      if not self.GetParam('skip_last_occ_cut'):
        occ,pdbcur=self.GetOccB( self.inp.Get('model',typ=model_type,filetype='pdb') )
        self.lastcut=self.CutOcc( occ,pdbcur, init=True )

  def RunPostprocess(self,restore=True,*args,**kwargs):
    if hasattr(self.ph,'report_fom'): 
      self.report_fom=self.ph.report_fom
      if not self.report_fom:
        self.result_str = "Wrong or no substructure: it's refinement resulted in FOM of 0. Structure solution unsuccessful."
        common.Error(self.result_str, nosuccess=True)
    if self.rv_report is not None:
      self.ph.rv_report=self.rv_report
      if self.ph.nick=='phas':
        self.ph.PrintActualLogGraph(finished=True)
      #self.rv_report.Text('&emsp;'+info, flush=True)
      if (self.lastcut or self.noimpr) and self.GetCrankParent() and self.GetCrankParent().prep.GetProcess('handdet'):
        self.GetCrankParent().prep.GetProcess('handdet').SetParam('ignore_passed_phas')
    process.RunPostprocess(self,restore,*args,**kwargs)
