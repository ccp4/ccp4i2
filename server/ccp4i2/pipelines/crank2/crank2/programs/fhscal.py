#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class fhscal(program):
  name="FHSCAL"
  binary="fhscal"
  ccp4_parsing=True
  stat={}
  #stat['mol_weight'] = common.stats(name='total molecular weight in Da', regexp=r"Total\s+Molecular\s+Weight\s*(\d+)\s+\(Da\)")

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    self.aver_data, self.anom_data_list, self.native = {}, [], None
    for f in self.inp.GetAll(filetype='mtz',typ='average',col='sigf'):
      if f.GetCrystalName()=='native' or 'native' in f.custom:
        if not self.native:
          self.native = f
      elif (f.GetCrystalName(),f.GetDataName()) not in self.aver_data:
        self.aver_data[(f.GetCrystalName(),f.GetDataName())] = f
    if len(self.aver_data)<1:
      common.Error('{0} requires at least two datasets for scaling.'.format(self.name))
    if not self.native:
      common.Error('{0} requires native dataset.'.format(self.name))
    self.SetArg( 'hklin', self.native.GetFileName() )
    self.SetKey( 'LABI', ( 'FP={0}'.format(self.native.GetLabel('f')), 'SIGFP={0}'.format(self.native.GetLabel('sigf'))) )
    for i,(cd,f) in enumerate(self.aver_data.items()):
      self.AddToKey( 'LABI', ( 'FPH{0}={1}'.format(i+1,f.GetLabel('f')), 'SIGFPH{0}={1}'.format(i+1,f.GetLabel('sigf'))) )
      fpl  = self.inp.Get(filetype='mtz',typ='plus', col='sigf',xname=f.GetCrystalName(),dname=f.GetDataName())
      fmin = self.inp.Get(filetype='mtz',typ='minus',col='sigf',xname=f.GetCrystalName(),dname=f.GetDataName())
      if fpl and fmin:
        self.anom_data_list.extend((fpl,fmin))
        self.AddToKey( 'LABI', ( 'FPH{0}(+)={1}'.format(i+1,fpl.GetLabel('f')),  'SIGFPH{0}(+)={1}'.format(i+1,fpl.GetLabel('sigf'))) )
        self.AddToKey( 'LABI', ( 'FPH{0}(-)={1}'.format(i+1,fmin.GetLabel('f')), 'SIGFPH{0}(-)={1}'.format(i+1,fmin.GetLabel('sigf'))) )


  def DefineOutput(self):
    for f in self.aver_data.values() + self.anom_data_list:
      self.inp.AddFileToChild( f, self.outfilename['mtz'] )

  def TreatOutput(self):
    self.SetArg('hklout', self.outfilename['mtz'])
