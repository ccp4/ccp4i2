#!/usr/bin/python
import os,sys,copy
from program import program
import common

class pdbmerge(program):
  name="PDB_MERGE"
  binary="pdb_merge"
  never_merge=True
  ccp4_parsing=True

  def TreatInput(self):
    models = self.inp.GetAll('model',filetype='pdb')
    self.models = [m  for i,m in enumerate(models)  if models.index(m)==i]
    if len(self.models)<2:
      common.Error('{0} models supplied for {1}, at least 2 needed for merging.'.format(
        len(self.models),self.name))
    self.SetArg('xyzin1', self.models[0].GetFileName('pdb'))
    self.SetArg('xyzin2', self.models[1].GetFileName('pdb'))

  def DefineOutput(self):
    outmdl=self.out.AddCopy(self.models[0])
    self.out.SetFileToChild(outmdl, self.name+'.pdb', 'pdb')
    typ=self.models[0].GetType()
    typ2=self.models[1].GetType()
    if ('substr' in typ or 'substr' in typ2) and ('partial' in typ or 'partial' in typ2):
      typ='partial+substr'
    if typ=='unknown':
      typ=typ2
    outmdl.SetType(typ)
    if 'substr' in typ and not 'substr' in typ2:
      atomtypes=copy.copy(self.models[1].atomtypes)
      atomtypes.update(self.models[0].atomtypes)
    else:
      atomtypes=copy.copy(self.models[0].atomtypes)
      atomtypes.update(self.models[1].atomtypes)
    atomtype1=self.models[0].atomtype1
    if not atomtype1:
      atomtype1=self.models[1].atomtype1
    outmdl.SetAtomTypes(atomtypes,atomtype1)

  def TreatOutput(self):
    self.SetArg('xyzout', self.out.Get('model').GetFileName('pdb'))
