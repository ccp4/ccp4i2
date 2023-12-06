#!/usr/bin/python
import os,sys
from process import process
import common
par=common.parameter

class dm(process):
  name="real space density modification"
  short_name="real space DM"
  supported_progs=["parrot", "solomon", "shelxe", "prasa"]
  supported_params = {}
  supported_params['solvent_content'] = par(desc='(Expected) solvent fraction of crystal', typ=float)
  supported_params['ncs_det'] = par( desc='Determine NCS from heavy atoms or partial model (Parrot only)', typ=bool )
  supported_params['ncs_det_ha'] = par( desc='Determine NCS from heavy atoms (Parrot only)', typ=bool )
  supported_params['ncs_det_mr'] = par( desc='Determine NCS from partial model (Parrot only)', typ=bool )
  supported_params['solventmask_radius'] = par( desc='Use the specified solvent mask radius (Parrot only)', typ=float )
  supported_params['solvent_perturb'] = par(desc='Adjust the solvent by the specified change', typ=(float,bool))


  def TreatInOutPar(self, set_all_par=False):
    if not self.GetProg(supported=True):
      if self.parent_process and self.parent_process.parent_process and \
         self.parent_process.parent_process.nick == 'handdet':
        self.AddProg('solomon')
      else:
        self.AddProg('parrot')
    prog=self.GetProg(supported=True)
    # changing the default ncs-mask-filter-radius for ha here!
    if prog.nick=='parrot' and not prog.IsKeyOrArg('ncs-mask-filter-radius') and \
       (self.GetParam('ncs_det') is not False and self.inp.Get('model',typ='substr',filetype='pdb')) and \
       (self.GetParam('ncs_det_ha') or not self.GetParam('ncs_det_mr')):
      prog.SetKey('ncs-mask-filter-radius', 9.0)
    # always using anisotropy for parrot as of now
    if prog.nick=='parrot':
      if not prog.IsKeyOrArg('anisotropy-correction'):
        prog.SetKey('anisotropy-correction',True)
    process.TreatInOutPar(self,set_all_par)

  ### disabled as map<->mtz is considered conversion and performed automatically on demand now!
  #def RunPreprocess(self, rundir=None, **kwargs):
    # let's transform to/from map/mtz if needed (at the moment, only for solomon)
  #  process.RunPreprocess(self,rundir,kwargs)
  #  if self.GetProg('solomon') and not self.inp.Get('mapcoef',filetype='map') and \
  #     self.inp.Get('mapcoef',filetype='mtz',col=('f','ph')):
  #      if not self.GetProg('fft'):
  #        self.AddProg('fft', ind=self.programs.index(self.GetProg('solomon')))
  #      if not self.GetProg('invfft'):
  #        self.AddProg('invfft', ind=self.programs.index(self.GetProg('solomon'))+1)

  def RunPreprocess(self, rundir=None, **kwargs):
    process.RunPreprocess(self,rundir,kwargs)
    # get solvent content estimation if not supplied
    if not self.GetParam('solvent_content') and not self.inp.Get(has_solvent_content=True) and \
       (self.inp.Get('sequence') or self.inp.Get(has_residues_mon=True)):
      matthews=self.AddProcess('matthews', propagate_out=False)
      matthews.Run()
    # ncs treatment (parrot only)
    if self.GetParam('ncs_det') is not False and self.GetProg('parrot'):
      if not self.IsInputtedParam('ncs_det',test_parents=True):
        mon_obj = self.inp.Get(has_monomers_asym=True)
        if mon_obj and mon_obj.monomers_asym==1 and (not mon_obj.seq_monomers or mon_obj.seq_monomers==1):
          self.SetParam('ncs_det',False)
          self.Info('One monomer in asymmetric unit estimated: NCS detection by Parrot disabled.', stdout=False)
      if self.GetParam('ncs_det'):
        self.GetProg('parrot').SetKey('ncs-operator',False,keep_previous=False)
    # if ncs_det is disabled then ha/mr suboptions can be ignored (and thus no warnings reported)
    if self.GetParam('ncs_det') is False:
      self.GetParam('ncs_det_ha')
      self.GetParam('ncs_det_mr')
