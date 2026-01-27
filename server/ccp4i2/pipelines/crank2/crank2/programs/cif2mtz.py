#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class cif2mtz(program):
  name="CIF2MTZ"
  binary="cif2mtz"

  def TreatInput(self):
    self.cif = self.inp.Get('fsigf',filetype='cif')
    if self.cif:
      self.SetArg('HKLIN',self.cif.GetFileName('cif'))
    else:
      common.Error('No cif file inputted to {0}'.format(self.name))
    if self.cif.GetCrystalName()!='Undefined':
      self.SetKey('name',('crystal',self.cif.GetCrystalName()))

  def DefineOutput(self):
    for o in set(self.inp.GetAll(filetype='cif',filename=self.cif.GetFileName('cif'))):
      out=self.out.Add(o)
      self.out.AddFileToChild(out, self.nick+'.mtz', 'mtz', no_rewrite=True)
      if out.GetLabel('i') and out.GetLabel('sigi'):
        out.SetLabel(('i','sigi'))
      if out.GetLabel('f') and out.GetLabel('sigf'):
        out.SetLabel(('f','sigf'))

  def TreatOutput(self):
    self.SetArg('hklout', self.out.Get().GetFileName('mtz'))
    for o in self.out.GetAll():
      if o.GetLabel('i') and o.GetLabel('sigi'):
        if o.GetType()=='plus':
          self.AddToKey('labout', ('I(+)='+o.GetLabel('i'), 'SIGI(+)='+o.GetLabel('sigi')))
        elif o.GetType()=='minus':
          self.AddToKey('labout', ('I(-)='+o.GetLabel('i'), 'SIGI(-)='+o.GetLabel('sigi')))
        elif o.GetType()=='average':
          self.AddToKey('labout', ('I='+o.GetLabel('i'), 'SIGI='+o.GetLabel('sigi')))
      if o.GetLabel('f') and o.GetLabel('sigf'):
        if o.GetType()=='plus':
          self.AddToKey('labout', ('F(+)='+o.GetLabel('f'), 'SIGF(+)='+o.GetLabel('sigf')))
        elif o.GetType()=='minus':
          self.AddToKey('labout', ('F(-)='+o.GetLabel('f'), 'SIGF(-)='+o.GetLabel('sigf')))
        elif o.GetType()=='average':
          self.AddToKey('labout', ('F='+o.GetLabel('f'), 'SIGF='+o.GetLabel('sigf')))
