#!/usr/bin/python
import os,sys
from program import program
import common

class pdbset(program):
  name="PDBSET"
  binary="pdbset"
  never_merge=True
  ccp4_parsing=True
  # pdbset cannot output these chain ids (complaining about longer than 1 character or not working properly...)
  bad_chain_ids=['&','\'','"']


  def Init(self):
    self.outfilename = { 'pdb': self.name+'.pdb' }

  def TreatInput(self):
    self.mdl = self.inp.Get('model',filetype='pdb')
    if not self.mdl:
      common.Error('No model supplied for {0}.'.format(self.name))
    self.SetArg('xyzin', self.mdl.GetFileName('pdb'))

  def DefineOutput(self):
    newmdl=self.out.AddCopy(self.mdl)
    self.out.SetFileToChild(newmdl, self.outfilename['pdb'], 'pdb')

  def TreatOutput(self):
    self.SetArg('xyzout', self.out.Get('model').GetFileName('pdb'))
