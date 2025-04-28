#!/usr/bin/python

from .. import common
from ..program import program


class matthews_coef(program):
  name="MATTHEWS_COEF"
  binary="matthews_coef"
  ccp4_parsing=True
  stat={}
  stat['solvent_content_all'] = common.stats(name='solvent content', xpath="result", attrib='percent_solvent', convert=True, multiple=True)
  stat['monomers_asym_all'] = common.stats(name='number of monomers in asymetric unit', xpath="result", attrib='nmol_in_asu', convert=True, multiple=True)
  stat['probab_matth_all'] = common.stats(name='Matthew coef. probability', xpath="result", attrib='prob_matth', convert=True, multiple=True)
  #stat['residues_mon'] = common.stats(name='number of residues in a monomer', xpath="parameters/residues")


  def TreatInput(self):
    # cell and spacegroup
    if not self.IsKey('cell'):
      self.cell=self.inp.Get('fsigf').GetCell(self)
      self.SetKey('cell', ' '.join((str(c) for c in self.cell)))
    if not self.IsKey('symm'):
      self.spgr=self.inp.Get('fsigf').GetSpacegroup(self)
      self.SetKey('symm', self.spgr.replace(' ',''))
    if not self.IsKey('molweight') and not self.IsKey('nres'):
      if self.inp.Get(has_residues_mon=True):
        self.SetKey('NRES',self.inp.Get(has_residues_mon=True).residues_mon)
      else:
        common.Error('Neither mol.weight nor number of residues inputted for {0}'.format(self.name))
    # number of monomers
    if self.inp.Get(has_monomers_asym=True) and not self.GetKey('NMOL'):
      self.SetKey('NMOL', self.inp.Get(has_monomers_asym=True).monomers_asym)
    else:
      self.SetKey('auto')

  def DefineOutput(self):
    self.out.AddNew('datafile',typ='xml',filename=self.name+'.xml',filetype='xml')

  def TreatOutput(self):
    self.SetKey('XMLOUTPUT')
    self.SetArg('XMLFILE', self.out.Get('datafile').GetFileName('xml'))
