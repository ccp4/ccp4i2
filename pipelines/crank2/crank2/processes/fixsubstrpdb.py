#!/usr/bin/python
import os,sys
from process import process
import common


class fixsubstrpdb(process):
  name="fixing substructure PDB file (due to programs treatment differences)"
  short_name="fixing sustructure PDB"
  no_reporting=True

  def RunBody(self,*args,**kwargs):
    pdbcur=self.GetOrAddProg('pdbcur')
    # set the atom type to the "main" atomtype
    if self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True):
      atomtype=self.inp.Get('model',typ=('substr','partial+substr'),has_atomtypes=True).GetAtomType()
      # we assume that the fixing is called from the same process that created the PDB
      if self.parent_process and (self.parent_process.GetProg('shelxd') or self.parent_process.GetProg('shelxe')):
        def_res='HAT'
        pdbcur.SetKey('renelement', "({1})/S[S] '{0}'".format(atomtype,def_res))
        pdbcur.SetKey('renatom', "({1})/S[{0}] '{0}'".format(atomtype,def_res))
      if self.parent_process and self.parent_process.GetProg('peakmax'):
        def_res='HOH'
        pdbcur.SetKey('renelement', "({1})/O[O] '{0}'".format(atomtype,def_res))
        pdbcur.SetKey('renatom', "({1})/O[{0}] '{0}'".format(atomtype,def_res))
      if self.parent_process and self.parent_process.GetProg('crunch2'):
        def_res='AA'
        pdbcur.SetKey('renelement', "({1})/S:A '{0}'".format(atomtype,def_res))
        pdbcur.SetKey('renatom', "({1})/S[{0}]:A '{0}'".format(atomtype,def_res))
      # these exceptions are due to refmac not accepting I/TA residue types...
      # S will be accepted soon:  the next 3 lines can be commented out then
      #if atomtype=='S':
      #  pdbcur.SetKey('renresidue', "[{0}]:* 'SO2'".format(atomtype))
      #  pdbcur.SetKey('renatom', "[{0}]:* '{0}'".format(atomtype))
      if atomtype=='I':
        pdbcur.SetKey('renresidue', "[{0}]:* 'IOD'".format(atomtype))
        pdbcur.SetKey('renatom', "[{0}]:* '{0}'".format(atomtype))
      elif atomtype=='AS':
        pdbcur.SetKey('renresidue', "[{0}]:* 'ARS'".format(atomtype))
        pdbcur.SetKey('renatom', "[{0}]:* '{0}'".format(atomtype))
      elif atomtype=='Y':
        pdbcur.SetKey('renresidue', "[{0}]:* 'YT3'".format(atomtype))
        pdbcur.SetKey('renatom', "[{0}]:* '{0}'".format(atomtype))
      #elif atomtype=='TA':
        #pdbcur.SetKey('renresidue', "[{0}]:* 'TBR'".format(atomtype))
        #pdbcur.SetKey('renatom', "(TBR)/{0}[{0}]:* 'TA1'".format(atomtype))
      elif atomtype=='U':
        pdbcur.SetKey('renresidue', "[{0}]:* 'U1'".format(atomtype))
        pdbcur.SetKey('renatom', "[{0}]:* '{0}'".format(atomtype))
      else:
        pdbcur.SetKey('renresidue', "[{0}]:* '{0}'".format(atomtype))
        pdbcur.SetKey('renatom', "[{0}]:* '{0}'".format(atomtype))
      # another refmac expception: early 5.8 refmac versions do not accept no chain ID
      pdbcur.SetKey('renchain', "/*// W")
      pdbcur.outfilename['pdb']='heavy.pdb'
      pdbcur.Run()
      if pdbcur.out.Get('model').GetType() not in ('substr','partial+substr'):
        pdbcur.out.Get('model').SetType('substr')
      self.CheckMergeChains(pdbcur.out.Get('model'))
    else:
      common.Warning('Substructure atom type not defined. This may be fatal for many non-shelx steps.')
    self.programs.remove(pdbcur)

  def CheckMergeChains(self, pdbobj=None):
    # merge substr. chains if there are too many (perventing possible crashes as happened eg for serca)
    # this is not needed if we do not have limits on the number of chains (eg mmcif)
    pdbcur=self.AddProg('pdbcur')
    pdbcur.runname += '_checkchains'
    pdbcur.SetKey('summarise')
    if pdbobj:
      pdbcur.inp.Set(pdbobj)
    pdbcur.Run()
    chain_ids = pdbcur.GetStat('chain_ids')
    if len(chain_ids)>8:
      pdbset=self.GetOrAddProg('pdbset')
      pdbset.outfilename['pdb']='heavy_chains_merged.pdb'
      pdbset.runname += '_mergechains'
      pdbset.SetKey('chain W')
      pdbset.SetKey('renumber 1')
      pdbset.inp.Set(pdbcur.out.Get('model'))
      pdbset.Run()
      if pdbset.out.Get('model').GetType() not in ('substr','partial+substr'):
        pdbset.out.Get('model').SetType('substr')
      self.programs.remove(pdbset)
    self.programs.remove(pdbcur)
