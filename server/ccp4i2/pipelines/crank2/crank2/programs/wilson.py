#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class wilson(program):
  name="WILSON"
  binary="wilson"
  stat={}
  stat['wilson_scale'] = common.stats(name='Wilson scale estimate', regexp=r"Least squares straight line gives:\s+B\s+=\s+\S+\s+SCALE\s+=\s+(\S+)")
  stat['wilson_B'] = common.stats(name='Wilson B estimate', regexp=r"Least squares straight line gives:\s+B\s+=\s+(\S+)\s+SCALE\s+=\s+\S+")

  def TreatInput(self):
    obj = self.inp.Get('fsigf', filetype='mtz', typ='average', col='f')
    if not obj:
      obj = self.inp.Get('fsigf', filetype='mtz', typ=('plus','minus'), col='f')
    if not obj:
      common.Error('No mtz with F\'s inputted to {0}'.format(self.name))
    self.SetArg('hklin', obj.GetFileName('mtz'))
    self.AddToKey('labin', ('FP='+obj.GetLabel('f'),'SIGFP='+obj.GetLabel('sigf')) )
    if self.inp.Get(has_residues_mon=True) and self.inp.Get(has_monomers_asym=True):
      self.SetKey('NRES',self.inp.Get(has_residues_mon=True).residues_mon*self.inp.Get(has_monomers_asym=True).monomers_asym)

