#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class gcx(program):
  name="GCX"
  binary="gcx"
  ccp4_parsing=True
  stat={}
  stat['solvent_content'] = common.stats(name='solvent content', xpath="parameters/solvent_content")
  stat['monomers_asym'] = common.stats(name='number of monomers in asymetric unit', xpath="parameters/nmonomers")
  stat['wilson_B'] = common.stats(name='Wilson\'s overall B factor', xpath="parameters/bfactor")
  stat['residues_mon'] = common.stats(name='number of residues in a monomer', xpath="parameters/residues")
  stat['monomers_asym_all'] = common.stats(name='number of monomers in asymetric unit', regexp=r' Monomers  Matthews-Coeff  Solvent-content    Probability  \$\$\s\$\$'+r'(?:\s+(\d+)\s+\S+\s+\S+\s+\S+\s)?'*99+r'\$\$')
  stat['solvent_content_all'] = common.stats(name='solvent content', regexp=r' Monomers  Matthews-Coeff  Solvent-content    Probability  \$\$\s\$\$'+r'(?:\s+\d+\s+\S+\s+(\S+)\s+\S+\s)?'*99+r'\$\$')
  stat['probab_matth_all'] = common.stats(name='Matthew coef. probability', regexp=r' Monomers  Matthews-Coeff  Solvent-content    Probability  \$\$\s\$\$'+r'(?:\s+\d+\s+\S+\s+\S+\s+(\S+)\s)?'*99+r'\$\$')


  def TreatInput(self):
    # input sequence
    seq=self.inp.Get('sequence')
    if seq and seq.GetFileName():
      # setting an absolute path here for backwards compatibility - there used to be a bug in gcx
      self.SetArg('seqin',seq.GetFileName(),relpath=False)
    elif self.inp.Get(has_residues_mon=True):
      self.SetKey('NRES',self.inp.Get(has_residues_mon=True).residues_mon)
    else:
      common.Error('Neither sequence nor number of residues inputted for {0}'.format(self.name))
    # data
    f, fmin, fpl = self.inp.Get('fsigf',filetype='mtz',typ='average'), self.inp.Get('fsigf',filetype='mtz',typ='minus'), self.inp.Get('fsigf',filetype='mtz',typ='plus')
    f_use = f if f else fmin
    #if not f:
     # f=self.inp.Get('fsigf',filetype='mtz',typ='plus')
      #fmin=self.inp.Get('fsigf',filetype='mtz',typ='minus')
    if not f and (not fmin or not fpl):
      common.Error('No X-ray data inputted for {0}'.format(self.name))
    self.model=self.inp.Get('model')
    if self.model: 
      if self.model.xname!=self.model.default_unknown:
        self.SetKey('XTAL', self.model.xname)
      if self.model.dname!=self.model.default_unknown:
        self.SetKey('DNAME', self.model.dname)
    if not self.IsKey('XTAL'):
      self.SetKey('XTAL', 'dummy')
    if not self.IsKey('DNAME'):
      self.SetKey('DNAME', 'dummy')
    self.SetKey('MTZIn',f_use.GetFileName('mtz'))
    if f_use.GetLabel('sigf'):
      if fmin and fpl and fmin.GetLabel('sigf') and fpl.GetLabel('sigf'):
        self.SetKey('COLU', ('F+='+fpl.GetLabel('f'),'SF+='+fpl.GetLabel('sigf')) )
        self.SetKey('COLU', ('F-='+fmin.GetLabel('f'),'SF-='+fmin.GetLabel('sigf')) )
      else:
        self.SetKey('COLU', ('F='+f.GetLabel('f'),'SF='+f.GetLabel('sigf')) )
    else:
      if fmin and fpl and fmin.GetLabel('sigi') and fpl.GetLabel('sigi'):
        self.SetKey('COLU', ('I+='+fpl.GetLabel('i'),'SI+='+fpl.GetLabel('sigi')) )
        self.SetKey('COLU', ('I-='+fmin.GetLabel('i'),'SI-='+fmin.GetLabel('sigi')) )
      else:
        self.SetKey('COLU', ('I='+f.GetLabel('i'),'SI='+f.GetLabel('sigi')) )
    # number of monomers
    if self.inp.Get(has_monomers_asym=True) and not self.GetKey('NMON'):
      self.SetKey('NMON', self.inp.Get(has_monomers_asym=True).monomers_asym)

  def DefineOutput(self):
    self.outputbase='gcx'
    self.out.AddNew('datafile',typ='xml',filename=self.outputbase+'.xml',filetype='xml')

  def TreatOutput(self):
    self.SetKey('outp {0}'.format(self.outputbase))
