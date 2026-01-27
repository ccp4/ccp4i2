#!/usr/bin/python
import os,sys
from ..program import program
from .. import common
import math

class prasa(program):
  name="PRASA"
  binary="prasa"
  #labelout_prefix="PRAS_"
  modif='-'
  always_merge=True
  stat={}
  stat['cc'] = common.stats(name='correlation coef. (all reflections mean)', 
    regexp=r"nd of trial \d+, final CC is (\S+)")
  stat['cc_range'] = common.stats(name='average correlation coef. in specified res. cutoff range', 
    regexp=r"nd of trial \d+, final CC is \S+\s+,\s+CCrange is (\S+)")
  stat['try'] = common.stats(name='try number', regexp=r"nd of trial (\d+), final CC is \S+")
  stat['cutoff'] = common.stats(name='used high resolution cutoff', regexp=r"\srescut:\s+(\S+)")
  references = ( "Skubak P (2018) Substructure determination using phase-retrieval techniques. Acta Cryst D74.", )

  def Init(self):
    self.SetArg('stdin')
    self.outfilename = { 'pdb': self.nick+'.pdb', \
                         'mtz': self.nick+'.mtz', \
                       }
    self.out_mapc = None

  def Interact_output(self, line, popen, empty):
    self.process.UpdateInfo(line,popen,prog=self)
    # this is done here in order to make sure that the same solution is not checked multiple times
    # polling to see whether prasa finished yet: if so then don't check so that we can proceed quickly (assuming checks are primarily to save time by preliminary stopping and don't need to be consistent)
    #if empty and self.process.check_solution[0]:
    if self.process.check_solution[0]>0 and popen.poll() is None:
      #self.process.CheckExtraStopConditionsPrasa(self.process.score,self.process.trial,self,popen)
      self.process.CheckExtraStopConditionsPrasa(self.process.check_solution[1],self.process.check_solution[0],self,popen)
      self.process.check_solution=(0,0)


  def Stop(self,popen=None):
    with open('stop_prasa', 'w') as f:
      f.write('stop_prasa')


  def TreatInput(self):
    self.ea = None
    if self.process.nick!='dm':
      self.ea = self.inp.Get('fsigf',typ='fa',filetype='mtz', inp_cont=self.inp.Get('fsigf',typ='average',filetype='mtz',try_convert=False))
      if not self.ea:
        self.ea = self.inp.Get('fsigf',typ='delta-anom',filetype='mtz')
    if not self.ea:
      self.ea = self.inp.Get('fsigf',typ='average',filetype='mtz')
    if self.ea is None:
      common.Error('FA/EA values not inputted to {0}'.format(self.name))
    self.SetKey('mtzin', self.ea.GetFileName('mtz'))
    if self.ea.GetLabel('e'):
      self.SetKey('colin-fo', self.ea.GetFullLabel('e','sige'))
    elif self.ea.GetLabel('f'):
      self.SetKey('colin-fo', self.ea.GetFullLabel('f','sigf'))
    if self.process.nick=='dm' and self.inp.Get('mapcoef',col='f',typ=('combined','best','weighted')):
      self.SetKey('colin-fc', self.inp.Get('mapcoef',col='ph',typ=('combined','best','weighted')).GetFullLabel('f','ph'))
    # atom type
    self.subs_inp = self.inp.Get('model',typ='substr',has_atomtypes=True,is_native=False)
    if self.subs_inp:
      self.SetKey('atom', self.subs_inp.GetAtomType())

  def TreatParams(self):
    if self.process.nick=='dm':
      if not self.GetKey('ntrials'):
        self.SetKey('ntrials', 1)
      if not self.GetKey('ncycles'):
        self.SetKey('ncycles', 8)
      if not self.GetKey('delta'):
        self.SetKey('delta', 2.0)
      if not self.GetKey('shannon'):
        self.SetKey('shannon', 1.5)
#      if not self.GetKey('chargeflip'):
#        self.SetKey('chargeflip', 2)
      if not self.GetKey('beta'):
        self.SetKey('beta', 0.35)
#      if not self.GetKey('recirestrdiff'):
#        self.SetKey('recirestrdiff', 0.00001)
#      if not self.GetKey('histmatch'):
#        self.SetKey('histmatch', 1)
    else:
      # expected number of atoms from substr. object
      #if self.inp.Get(typ='substr',has_num_atoms=True) and not self.IsKey('natoms'):
      #  self.SetKey('natoms', self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms)
      if not self.IsKey('ncycles') and ((self.process.GetParam('num_atoms') and self.process.GetParam('num_atoms')>=20) or
         (self.inp.Get(typ='substr',has_num_atoms=True) and self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms>=20)):
        if not self.process.IsInputtedParam('num_trials') and self.process.GetParam('num_trials') and not self.GetKey('pdbin'):
          self.SetKey('ncycles', 750)
          if not self.IsKey('ntrials'):
            self.SetKey('ntrials', int(self.process.GetParam('num_trials')*0.5), keep_previous=False)
          else:
            self.SetKey('ntrials', int(self.GetKey('ntrials')*0.5), keep_previous=False)
        else:
          self.SetKey('ncycles', 250)
      if self.process.GetVirtPar('num_trials') and not self.GetKey('ntrials') and not self.GetKey('pdbin'):
        self.SetKey('ntrials', self.process.GetVirtPar('num_trials'))
      if self.process.GetVirtPar('high_res_cutoff') and not self.GetKey('rescut'):
        self.SetKey('rescut', self.process.GetVirtPar('high_res_cutoff'))
      if self.process.GetParam('num_atoms') and not self.GetKey('natoms'):
        self.SetKey('natoms', self.process.GetParam('num_atoms'))
      if not self.GetKey('minpeaks'): # or not not self.GetKey('natoms'):
        num_at = self.process.GetParam('num_atoms')
        if not num_at and self.inp.Get(typ='substr',has_num_atoms=True):
          num_at = self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms
        if num_at:
          if num_at>13 and not self.GetKey('minpeaks'):
            if num_at<=40:
              self.SetKey('minpeaks', int(0.2*num_at))  #eg ssec
            else:
              self.SetKey('minpeaks', 8+int(0.1*(num_at-40)))  # eg 8dop
          if not self.GetKey('natoms'):
            self.SetKey('natoms', int(3*num_at))
      if self.process.IsParam('min_dist_symm_atoms') and not self.IsKey('specialpos'):
        #self.SetKey('specialpos', min(0.6,int(bool(not(self.process.GetParam('min_dist_symm_atoms')>0)))))
        self.SetKey('specialpos', int(bool(not(self.process.GetParam('min_dist_symm_atoms')>0))))
      if self.process.GetParam('num_threads') and not self.GetKey('numthreads'):
        self.SetKey('numthreads', self.process.GetParam('num_threads'))
    program.TreatParams(self)

  def DefineOutput(self):
    #self.out.AddNew( 'mapcoef', self.nick+'.map', filetype='map', typ='anomalous', xname=self.ea.GetCrystalName(), dname=self.ea.GetDataName() )
    if not self.GetKey('generaterefdist'):
      if self.subs_inp:
        subs_out=self.out.AddCopy(self.subs_inp)
        self.out.AddFileToChild( subs_out, self.outfilename['pdb'], filetype='pdb' )
      else:
        self.out.AddNew( 'model', self.outfilename['pdb'], filetype='pdb', typ='substr', xname=self.ea.GetCrystalName() )
    if self.process.nick=='dm':
      self.out_mapc = self.out.AddNew( 'mapcoef', self.outfilename['mtz'] )
      self.out_mapc.SetType('densmod')
      self.out_mapc.SetLabel( ['f','ph'], ['pras.F_phi.F','pras.F_phi.phi'] )  # todo: adjustable output col labels
      self.out_mapc.SetCrystalName(self.ea.GetCrystalName())
      self.out_mapc.SetDataName(self.ea.GetDataName())

  def TreatOutput(self):
    #self.SetKey( 'mapout', self.out.mapcoef[-1].GetFileName('map') )
    if not self.GetKey('generaterefdist'):
      self.SetKey( 'pdbout', self.out.model[-1].GetFileName('pdb') )
    if self.process.nick=='dm':
      self.SetKey( 'mtzout', self.out_mapc.GetFileName('mtz') )

  def GetTrialPdb(self,trial):
    # returns input trial's pdb name (incl. path) 
    # must be called after DefineOutput(), otherwise there is an error.
    return self.out.Get('model').GetFileName()+"_"+str(trial)
