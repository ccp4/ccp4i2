#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class coordconv(program):
  name="COORDCONV"
  binary="coordconv"
  never_merge=True
  ccp4_parsing=True
  # mapping between the filetype ('inputformat' parameter) and coordconv 'input' key value
  in_key = { 'pdb': 'pdb', 'frac': 'frac', 'hat': 'shelx-s', 'res': 'shelx-s' }
  # mapping between the filetype ('outputformat' parameter) and coordconv 'output' key value
  out_key = { 'pdb': 'pdb', 'frac': 'frac' }

  def Init(self):
    self.outfilename = { 'pdb': self.name+'.pdb', 'frac': self.name+'.frac' }

  def TreatInput(self):
    in_ft = next((inf for inf,i in self.in_key.items() if i==self.GetKey('input')), None)
    if not in_ft:
      common.Error('Wrong input format "{1}" specified for {0}: only one of {2} supported'.format(self.name,self.GetKey('input'),self.in_key.values()))
    self.model = self.inp.Get('model',filetype=in_ft)
    if self.model:
      self.SetArg('xyzin', self.model.GetFileName(in_ft))
    else:
      common.Error('No such input file in format {1} inputted to {0}'.format(self.name,in_ft))

  def TreatParams(self):
    if self.process and 'outputformat' in self.process.supported_params:
      if not self.process.IsParam('outputformat') and not self.IsKey('output'):
        common.Error('Output file type not specified for {0}'.format(self.name))
      elif self.process.GetParam('outputformat') not in self.out_key:
        common.Error('Output format {1} not supported by {0}'.format(self.name, self.process.GetParam('outputformat')))
      elif not self.IsKey('output'):
        self.SetKey('output', self.out_key[self.process.GetParam('outputformat')])
      else:
        if self.GetKey('output') in self.out_key.itervalues:
          of = next(of for of,o in self.out_key.items() if o==self.GetKey('output'))
          self.process.SetParam('outputformat', of)
        else:
          common.Warning('Output "{0}" not predefined for {1}'.format(self.GetKey('output'),self.name))
    if self.process and 'inputformat' in self.process.supported_params:
      if not self.process.IsParam('inputformat') and not self.IsKey('input'):
        common.Error('Input file type not specified for {0}'.format(self.name))
      elif not self.IsKey('input'):
        self.SetKey('input', self.in_key[self.process.GetParam('inputformat')])
      else:
        if self.GetKey('input') in self.in_key.itervalues:
          inf = next(inf for inf,i in self.in_key.items() if i==self.GetKey('input'))
          self.process.SetParam('inputformat', inf)
        else:
          common.Warning('Input "{0}" not predefined for {1}'.format(self.GetKey('input'),self.name))
    program.TreatParams(self)

  def DefineOutput(self):
    self.out_ft = next((of for of,o in self.out_key.items() if o==self.GetKey('output')), None)
    if not self.out_ft:
      common.Error('Wrong output format "{1}" specified for {0}: only one of {2} supported'.format(self.name,self.GetKey('output'),self.out_key.values()))
    outcont=self.out.Add(self.model)
    self.out.AddFileToChild(outcont, self.outfilename[self.out_ft], self.out_ft, no_rewrite=True)

  def TreatOutput(self):
    self.SetArg('xyzout', self.out.Get().GetFileName(self.out_ft))
