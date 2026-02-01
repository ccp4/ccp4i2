#!/usr/bin/python
import os,sys,re,string
from ..process import process
from .. import common


class sepsubstrprot(process):
  name="separation of substructure and protein into 2 PDB files"
  short_name="PDB models separation"
  # assumes substr+prot pdb and substructure atom types on input 
  # output is a substructure pdb and a protein pdb (with empty intersection between them)
  # at the moment, just performed in a simplistic way, based on substructure atomtypes
  # this may be not indended eg for S+methionine/cysteine - perhaps (some previous) substr pdb 
  # should be inputted too and distance criterion from previous substr atoms added?

  def RunBody(self,*args,**kwargs):
    if not self.GetProg('pdbcur'):
      self.pdbcur=self.AddProg('pdbcur')
    else:
      self.pdbcur=self.GetProg('pdbcur')
    if not self.inp.Get('model',typ='partial+substr'):
      common.Error('No protein+substructure model inputted - separation cannot be done')
    if not self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True):
      common.Error('Substructure atom types not specified - separation of substructure and protein cannot be done')
    # quick hack to fix the TER card issues by removing them - not needed anymore thanks to Eugene's new pdbcur functions!
    #pdb = self.inp.Get('model',typ='partial+substr')
    #file_fixed = pdb.GetFileName('pdb')+'_fixed.pdb'
    #with open(file_fixed,'w') as g:
    #  with open(pdb.GetFileName('pdb')) as f:
    #    for line in f:
    #      if not line.startswith('TER   '):
    #        g.write(line)
    #pdb_fixed=self.pdbcur.inp.AddCopy(pdb)
    #self.pdbcur.inp.SetFileToChild(pdb_fixed,file_fixed,filetype='pdb')
    # if substr. is provided in inp. then also filter according to the chains in it
    substr_chains_filter=''
    if self.inp.Get('model',typ='substr'):
      self.pdbcur.ClearAnyParams()
      self.pdbcur.runname=self.pdbcur.name+'_checkinpsub'
      self.pdbcur.inp.Set(self.inp.Get('model',typ='substr'))
      self.pdbcur.SetKey('summ')
      self.pdbcur.Run()
      substr_chains=self.pdbcur.GetStat('chain_ids',accept_none=True)
      self.pdbcur.inp.Set(self.inp.Get('model',typ='partial+substr'))
      if set(self.pdbcur.bad_chain_ids).intersection(substr_chains):
        # rename substr. chains if ids that pdbcur does not accept are present
        self.pdbcur.runname=self.pdbcur.name+'_checkinpmerge'
        self.pdbcur.Run()
        all_chains=self.pdbcur.GetStat('chain_ids',accept_none=True)
        pdbset=self.GetOrAddProg('pdbset')
        pdbset.ClearAnyParams()
        pdbset.inp.Set(self.inp.Get('model',typ='partial+substr'))
        for bad_id in set(self.pdbcur.bad_chain_ids).intersection(substr_chains):
          new_id = next(nid for nid in string.printable if nid not in all_chains and nid not in self.pdbcur.bad_chain_ids and nid not in self.pdbset.bad_chain_ids)
          pdbset.SetKey('chain',(bad_id,new_id))
          all_chains.append(new_id)
          substr_chains.append(new_id), substr_chains.remove(bad_id)
        pdbset.Run()
        self.pdbcur.inp.Set(pdbset.out.Get('model'))
        substr_chains_filter='//'+','.join(substr_chains)+'//' if substr_chains else ''
    # if the prot+substr input model does not contain ha types then include them from substr
    if not self.pdbcur.inp.Get('model',has_atomtypes=True):
      substr=self.inp.Get('model',typ='substr',has_atomtypes=True)
      self.pdbcur.inp.Get('model').SetAtomTypes( substr.GetAtomTypes(), atomtype1=substr.GetAtomType() )
    # first get heavy atoms PDB
    self.pdbcur.ClearAnyParams()
    self.pdbcur.runname=self.pdbcur.name+'_heavy'
    at_str=','.join(self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True).GetAtomTypes())
    self.pdbcur.SetKey('delter')
    # at the moment, an exception is needed for TA as its name is TA1 in the only residue existing
    # we could only check the atom type in the future but for the moment this is safer
    #if at_str=='TA':
    #  self.pdbcur.SetKey('lvatom', '"{1}{0}1[{0}]:*"'.format(at_str,substr_chains_filter))
    #else:
    self.pdbcur.SetKey('lvatom', '"{1}{0}[{0}]:*"'.format(at_str,substr_chains_filter))
      #self.pdbcur.SetKey('lvatom', '"(!MSE)/{0}[{0}]:*"'.format(at_str))
    self.pdbcur.outfilename['pdb']='heavy.pdb'
    self.pdbcur.Run()
    self.pdbcur.out.Get('model').SetType('substr')
    # manually remove all cards before CRYST1 (otherwise wrong cards eg non-existing bonds may remain 
    # and cause trouble later - pdbcur should remove them but does not)
    with open(self.pdbcur.out.Get('model').GetFileName('pdb')) as f:
      pdb=f.read()
      re_res=re.search('CRYST1.+(\n.*)+',pdb)
      if not re_res:
        common.Error('CRYST1 card missing in PDB file {0}'.format(self.pdbcur.out.Get('model').GetFileName('pdb')))
    with open(self.pdbcur.out.Get('model').GetFileName('pdb'),'w') as g:
      g.write(re_res.group(0))
    fixsubstr=self.AddProcess('fixsubstrpdb')
    fixsubstr.CheckMergeChains(self.pdbcur.out.Get('model'))
    self.processes.remove(fixsubstr)
    # now get partial model PDB
    self.pdbcur.ClearAnyParams()
    self.pdbcur.runname=self.pdbcur.name+'_part'
    for at in self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True).GetAtomTypes():
      self.pdbcur.SetKey('rmatom', '"{1}{0}[{0}]:*"'.format(at,substr_chains_filter))
      if at=='TA':
        self.pdbcur.SetKey('rmatom', '"{1}{0}1[{0}]:*"'.format(at,substr_chains_filter))
    self.pdbcur.outfilename['pdb']='part.pdb'
    self.pdbcur.Run()
    self.pdbcur.out.Get('model').SetType('partial')
