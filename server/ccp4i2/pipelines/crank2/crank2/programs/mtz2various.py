#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class mtz2various(program):
  name="MTZ2VARIOUS"
  binary="mtz2various"
  ccp4_parsing=True
  # mapping between the filetype ('outputformat' parameter) and mtz2various 'output' key value
  out_key = { 'hkl': 'shelx', 'sca': 'scal', 'phs': "user '(3I5,F9.2,F8.4,F8.1)'" }

  def Init(self):
    self.outfilename = {}
    self.opts = self.process.opts  if self.process and hasattr(self.process,'opts')  else {}
    for out in self.out_key: 
      self.outfilename[out] = self.name+'.'+out

  def TreatInput(self):
    if self.GetKey('output').startswith('sca'):
      self.inpcont=self.inp.GetAll(filetype='mtz',col='i')
      if 'amplitudes' in self.opts:
        common.Error('{} asked to convert amplitudes but {} only supports intensities'.format(self.name,self.GetKey('output')))
    else:
      self.inpcont=self.inp.GetAll(filetype='mtz')
    inten = True  if not 'amplitudes' in self.opts and self.GetKey('output')!=self.out_key['phs']  else False
    if self.inpcont:
      self.SetArg('hklin',self.inpcont[0].GetFileName('mtz'))
      for cont in self.inpcont:
        # intensities get priority over f's unless specified otherwise
        if cont.GetType()=='plus' and cont.GetLabel('i') and inten:
          self.AddToKey('labin', ('I(+)='+cont.GetLabel('i'), 'SIGI(+)='+cont.GetLabel('sigi')))
        elif cont.GetType()=='plus' and cont.GetLabel('f'):
          self.AddToKey('labin', ('F(+)='+cont.GetLabel('f'), 'SIGF(+)='+cont.GetLabel('sigf')))
        elif cont.GetType()=='minus' and cont.GetLabel('i') and inten:
          self.AddToKey('labin', ('I(-)='+cont.GetLabel('i'), 'SIGI(-)='+cont.GetLabel('sigi')))
        elif cont.GetType()=='minus' and cont.GetLabel('f'):
          self.AddToKey('labin', ('F(-)='+cont.GetLabel('f'), 'SIGF(-)='+cont.GetLabel('sigf')))
        elif cont.GetType()=='average' and cont.GetLabel('i') and inten:
          self.AddToKey('labin', ('I='+cont.GetLabel('i'), 'SIGI='+cont.GetLabel('sigi')))
        elif cont.GetType()=='average' and cont.GetLabel('f'):
          if self.GetKey('output')==self.out_key['phs']:
            self.AddToKey('labin', 'FP='+cont.GetLabel('f'), ind=0)
          else:
            self.AddToKey('labin', ('FP='+cont.GetLabel('f'), 'SIGFP='+cont.GetLabel('sigf')))
        elif cont.GetType() in ('best','combined') and cont.GetLabel('ph') and cont.GetLabel('fom'):
          self.AddToKey('labin', ('FOM='+cont.GetLabel('fom'), 'PHIB='+cont.GetLabel('ph')))
    else:
      common.Error('No suitable mtz file inputted to {0}'.format(self.name))

  def TreatParams(self):
    if not self.process.IsParam('outputformat') and not self.IsKey('output'):
      common.Error('Output file type not specified for {0}'.format(self.name))
    elif not self.IsKey('output'):
      self.SetKey('output', self.out_key[self.process.GetParam('outputformat')])
    else:
      if self.GetKey('output') in out_key.itervalues:
        of = [of for of,o in out_key if o==self.GetKey('output')][0]
        self.process.SetParam('outputformat', of)
      else:
        common.Warning('Output "{0}" not predefined for {1}'.format(self.GetKey('output'),self.name))
    program.TreatParams(self)

  def DefineOutput(self):
    outtyp = self.process.GetParam('outputformat')
    for cont in self.inpcont:
      if outtyp!='phs' or cont.nick=='mapcoef':
        outcont=self.out.Add(cont)
        self.out.AddFileToChild(outcont, self.outfilename[outtyp], outtyp, no_rewrite=True)

  def TreatOutput(self):
    self.SetArg('hklout', self.out.Get().GetFileName(self.process.GetParam('outputformat')))
