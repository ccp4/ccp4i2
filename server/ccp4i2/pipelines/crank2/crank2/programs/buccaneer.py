#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class buccaneer(program):
  name="Buccaneer"
  binary="cbuccaneer"
  labelout_prefix="BUCC_"
  modif='-'
  always_merge=True
  stat={}
  stat['res_built'] = common.stats(name='number of residues built', regexp=r"(\d+)\s+residues were built in\s+\d+\s+fragments,")
  stat['frag_built'] = common.stats(name='number of fragments built', regexp=r"\d+\s+residues were built in\s+(\d+)\s+fragments,")
  stat['compl_res'] = common.stats(name='completness by residues built', regexp=r"Completeness by residues built:\s+(\S+)\%")
  stat['compl_chain'] = common.stats(name='completness of chains', regexp=r"Completeness of chains \(number\):\s+(\S+)\%")
  references = ( "Cowtan K. (2006) The Buccaneer software for automated model building. 1. "+
                 "Tracing protein chains. Acta Cryst. D62, 1002-1011.", )

  def Init(self):
    self.SetArg('stdin')

  def TreatInput(self):
    # F,sigF - prefer native
#    self.fsf = self.inp.Get('fsigf',typ=('average','average-derived'),col='f',is_native=True)
    self.fsf = self.inp.Get('fsigf',typ=('average'),col='f',is_native=True)
    if not self.fsf:
      self.fsf = self.inp.Get('fsigf',typ=('average-derived'),col='f',is_native=True)
    if not self.fsf:
      self.fsf = self.inp.Get('fsigf',typ=('average'),col='f')
    if not self.fsf:
      self.fsf = self.inp.Get('fsigf',typ=('average-derived'),col='f')
    if self.fsf:
      self.SetKey('mtzin-wrk', self.fsf.GetFileName('mtz'))
      self.SetKey('colin-wrk-fo', self.fsf.GetFullLabel('f','sigf'))
    else:
      common.Error('No data inputted to {0}'.format(self.name))
    # prefer combined phases over "best" phases
    for mapc in self.inp.GetAll('mapcoef',typ='combined',filetype='mtz',inp_cont=self.fsf) + \
                self.inp.GetAll('mapcoef',typ='best',filetype='mtz',inp_cont=self.fsf):
      if mapc.GetLabel('hla'):
        self.SetKey('colin-hl', mapc.GetFullLabel('hla','hlb','hlc','hld'))
        break
      if mapc.GetLabel('ph') and mapc.GetLabel('fom'):
        self.SetKey('colin-phifom', mapc.GetFullLabel('ph','fom'))
        break
    else:
      common.Error('No phase distribution inputted for {0}'.format(self.name))
    # free assignment
    free=self.inp.Get('exclude',typ='freeR',col='free')
    if free and (not self.IsKey('colin-free') or self.GetKey('colin-free') is True):
      self.SetKey( 'colin-free', free.GetFullLabel('free'), keep_previous=False )
    # input map (optional - phase distr. and Fo used if not inputted)
    mapc = self.inp.Get('mapcoef', col=('f','ph'), custom='build')
    if not mapc and self.process.GetVirtPar('from_weighted_map'):
      mapc = self.inp.Get('mapcoef',typ='weighted',col=('f','ph'))
    if mapc:
      self.SetKey('colin-wrk-fc', mapc.GetFullLabel('f','ph'))
    # sequence
    if self.inp.Get('sequence') and self.inp.Get('sequence').GetFileName():
      self.SetKey('seqin-wrk', self.inp.Get('sequence').GetFileName())
    else:
      if not self.GetKey('find'):
        self.SetKey('find'),self.SetKey('grow'),self.SetKey('join'),self.SetKey('link')
        self.SetKey('filter'),self.SetKey('prune'),self.SetKey('rebuild'),self.SetKey('tidy') # should test whether >=1.6.3?
      else:
        self.SetKey('sequence',False),self.SetKey('correct',False),self.SetKey('ncsbuild',False)
    if self.inp.Get('model',typ=('partial','partial+substr'),filetype='pdb'):
      if self.process.GetVirtPar('from_mr_model'):
        self.SetKey('pdbin-mr', self.inp.Get('model',typ=('partial','partial+substr')).GetFileName('pdb'))
      else:
        self.SetKey('pdbin-wrk', self.inp.Get('model',typ=('partial','partial+substr')).GetFileName('pdb'))
        # this is not needed anymore (set by buc. automatically now) but still keeping for case older versions were used
        if not self.IsKey('correlation-mode') and not self.IsKey('no-correlation-mode'):
          self.SetKey('correlation-mode', True)

  def DefineOutput(self):
    out_model = self.out.AddNew( 'model', self.nick+'.pdb', typ='partial', xname=self.fsf.xname )
    if 'native' in self.fsf.custom:
      out_model.custom.append('native')

  def TreatOutput(self):
    self.SetKey( 'pdbout-wrk', self.out.model[-1].GetFileName('pdb') )
