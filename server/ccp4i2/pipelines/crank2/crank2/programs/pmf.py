#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class pmf(program):
  name="PMF"
  binary="pmf"
  ccp4_parsing=True

  def TreatInput(self):
    fa = self.inp.Get('fsigf',typ='fa',col='f')
    ea = self.inp.Get('fsigf',typ='fa',filetype='drear',inp_cont=fa)
    if ea is None:
      common.Error('FA/EA values not inputted to {0}'.format(self.name))
    self.SetArg('HKLIN', ea.GetFileName('drear'))
    # set cell and spacegroup
    if not self.IsKey('CELL') or not self.IsKey('SYMM'):
      self.SetKey('CELL', ' '.join((str(c) for c in fa.GetCell(self))))
      self.SetKey('SYMM', ea.GetSpacegroupNumber(self))
    # expected number of atoms from substr. object
    if self.inp.Get(typ='substr',has_num_atoms=True) and not self.GetKey('NATO'):
      self.SetKey('NATO', self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms)

  def DefineOutput(self):
    # not sure what patterson attributes should patterson have nor whether it should exist at all...
    self.out.AddNew('model', self.name+'.xyz', filetype='xyz', typ='patterson')

  def TreatOutput(self):
    self.SetArg('xyzout', self.out.Get('model').GetFileName('xyz'))
