#!/usr/bin/python
import os,sys
from program import program
import common

class pdbcur(program):
  name="PDBCUR"
  binary="pdbcur"
  never_merge=True
  ccp4_parsing=True
  # pdbcur cannot select these chain ids (due to syntax restrictions)
  bad_chain_ids=['[',']','*','&','-']
  stat={}
  stat['atoms_deleted'] = common.stats(name='(list of) numbers of atoms deleted', multiple=True,
            regexp=r"(\S+)\s+atoms were deleted.")
  stat['atoms_deleted_left'] = common.stats(name='(list of) numbers of atoms deleted and left', multiple=True,
            regexp=r"(\S+)\s+atoms were deleted,\s+(\S+)\s+left.")
  stat['atoms_deleted_occ'] = common.stats(name='message about atoms deleted due to low occupancy',
            regexp=r"Atoms with occupancy less than or equal to\s+\S+\s+have been removed.")
  stat['chain_ids'] = common.stats(name='list of all chain IDs in the input pdb file', multiple=True,
            regexp=r"Chain \"(\S+)\" has \d+ residues, of which", convert=False)
  # this is regexp for the PDB, not log
  stat['occup'] = common.stats(name='atom occupancies', regexp=r'(?:(?:ATOM  )|(?:HETATM)).{{50}}(\d+.?\d*)', multiple=True)
  stat['B'] = common.stats(name='atomic B factors', regexp=r'(?:(?:ATOM  )|(?:HETATM)).{{54}}\s?(\d+.?\d*)', multiple=True)
  stat['spgr'] = common.stats(name='spacegroup', regexp=r'(?:CRYST1).{{49}}(\S+\s*\S*\s*\S*\s*\S*)')


  def Init(self):
    self.outfilename = { 'pdb': self.name+'.pdb' }

  def TreatInput(self):
    self.mdl = self.inp.Get('model',filetype='pdb')
    if not self.mdl:
      common.Error('No model supplied for {0}.'.format(self.name))
    self.SetArg('xyzin', self.mdl.GetFileName('pdb'))

  def DefineOutput(self):
    # the typ and atomtypes could change in principle; for the moment just this simplistic implementation
    newmdl=self.out.AddCopy(self.mdl)
    self.out.SetFileToChild(newmdl, self.outfilename['pdb'], 'pdb')

  def TreatOutput(self):
    self.SetArg('xyzout', self.out.Get('model').GetFileName('pdb'))
