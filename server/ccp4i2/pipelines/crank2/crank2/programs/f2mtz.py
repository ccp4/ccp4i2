#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class f2mtz(program):
  name="F2MTZ"
  binary="f2mtz"
  ccp4_parsing=True

  # warning:  shelx outputs assumed as of now!  the formats should be coupled with programs later.
  def TreatInput(self):
    inp_format = self.process.GetParam('inputformat')
    self.inp_cont = None
    if inp_format == 'hkl':
      self.inp_cont = self.inp.Get(filetype=inp_format,typ='fa',try_convert=False)
      if self.inp_cont:
        self.SetKey('format','\'(3I4,2F8.2,I4)\'')
      else:
        self.inp_cont = self.inp.Get(filetype=inp_format,typ='average',try_convert=False)
        if self.inp_cont:
          self.SetKey('format','\'(3I4,2F8.2)\'')
    elif inp_format == 'phs':
      self.inp_cont = self.inp.Get(filetype=inp_format,try_convert=False)
      if self.inp_cont:
        self.SetKey('format','\'(3I4,2F9.2,F8.1,F9.2)\'')
    if self.inp_cont:
      self.SetArg('hklin', self.inp_cont.GetFileName(inp_format))
    else:
      common.Error('No suitable file inputted to {0}'.format(self.name))
    symm, cell = self.inp_cont.GetSpacegroup(self,accept_none=True), self.inp_cont.GetCell(self,accept_none=True)
    if not symm or not cell:
      inp = self.inp.Get('fsigf',xname=self.inp_cont.GetCrystalName(),try_convert=False)
      if not inp:
        inp = self.inp.Get('fsigf',try_convert=False)
      if inp:
        if not symm:
          symm = inp.GetSpacegroup(self)
        if not cell:
          cell = inp.GetCell(self)
    if symm and cell:
      self.SetKey('symm', symm.replace(" ",""))
      self.SetKey('cell', ','.join(str(c) for c in cell))
    else:
      common.Error('{0} could not determine cell/symmetry'.format(self.name))
    if self.inp_cont.GetCrystalName()!='Undefined':
      self.SetKey('name',('crystal',self.inp_cont.GetCrystalName()))


  def DefineOutput(self):
    self.outmtz=self.out.Add(self.inp_cont)
    self.out.AddFileToChild(self.outmtz, self.nick+'.mtz', 'mtz', no_rewrite=True)
    if self.inp_cont.GetType()=='fa':
      self.outmtz.SetLabel(('f','sigf','alpha'))
    elif self.inp_cont.GetType()=='average':
      self.outmtz.SetLabel(('i','sigi'))
    elif self.process.GetParam('inputformat')=='phs':
      self.outmtz.SetLabel(('f','ph','fom'))

  def TreatOutput(self):
    self.SetArg('hklout', self.out.Get().GetFileName('mtz'))
    if self.inp_cont.GetType()=='fa':
      self.AddToKey('labout', ('H','K','L',self.outmtz.GetLabel('f'), self.outmtz.GetLabel('sigf'), self.outmtz.GetLabel('alpha')))
      self.AddToKey('ctype', ('H','H','H','F','Q','P'))
    elif self.inp_cont.GetType()=='average':
      self.AddToKey('labout', ('H','K','L',self.outmtz.GetLabel('i'), self.outmtz.GetLabel('sigi')))
      self.AddToKey('ctype', ('H','H','H','J','Q'))
    elif self.process.GetParam('inputformat')=='phs':
      self.AddToKey('labout', ('H','K','L',self.outmtz.GetLabel('f'), self.outmtz.GetLabel('fom'), self.outmtz.GetLabel('ph'), 'SIGF'))
      self.AddToKey('ctype', ('H','H','H','F','W','P','Q'))
