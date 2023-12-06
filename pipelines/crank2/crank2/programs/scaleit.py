#!/usr/bin/python
import os,sys
from program import program
import common

class scaleit(program):
  name="SCALEIT"
  binary="scaleit"
  ccp4_parsing=True

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }
    self.revert_if_one = False

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
    # reverting native and derivative in case we only have one, so that native is scaled and anom. untouched!
    if len(self.aver_data)==1 and self.revert_if_one:
      self.native, self.aver_data = self.aver_data[list(self.aver_data.keys())[0]], { (self.native.GetCrystalName(),self.native.GetDataName()): self.native }
    self.SetArg( 'hklin', self.native.GetFileName() )
    self.SetKey( 'LABI', ( 'FP={0}'.format(self.native.GetLabel('f')), 'SIGFP={0}'.format(self.native.GetLabel('sigf'))) )
    for i,(cd,f) in enumerate(self.aver_data.items()):
      self.AddToKey( 'LABI', ( 'FPH{0}={1}'.format(i+1,f.GetLabel('f')), 'SIGFPH{0}={1}'.format(i+1,f.GetLabel('sigf'))) )
      if f.GetLabel('sigi'):
        self.AddToKey( 'LABI', ( 'IMEAN{0}={1}'.format(i+1,f.GetLabel('i')), 'SIGIMEAN{0}={1}'.format(i+1,f.GetLabel('sigi'))) )
      # anomalous data deterioriate after scaling by scaleit! (eg fibro cannot be solved)
      # there might be issues though if we scale averages and +/- are not scaled - a better solution should be found
      #fpl  = self.inp.Get(filetype='mtz',typ='plus', col='sigf',xname=f.GetCrystalName(),dname=f.GetDataName())
      #fmin = self.inp.Get(filetype='mtz',typ='minus',col='sigf',xname=f.GetCrystalName(),dname=f.GetDataName())
      #if fpl and fmin:
      #  self.anom_data_list.extend((fpl,fmin))
      #  self.AddToKey( 'LABI', ( 'FPH{0}(+)={1}'.format(i+1,fpl.GetLabel('f')),  'SIGFPH{0}(+)={1}'.format(i+1,fpl.GetLabel('sigf'))) )
      #  self.AddToKey( 'LABI', ( 'FPH{0}(-)={1}'.format(i+1,fmin.GetLabel('f')), 'SIGFPH{0}(-)={1}'.format(i+1,fmin.GetLabel('sigf'))) )
      #  if fpl.GetLabel('sigi') and fmin.GetLabel('sigi'):
      #    self.AddToKey( 'LABI', ( 'I{0}(+)={1}'.format(i+1,fpl.GetLabel('i')),  'SIGI{0}(+)={1}'.format(i+1,fpl.GetLabel('sigi'))) )
      #    self.AddToKey( 'LABI', ( 'I{0}(-)={1}'.format(i+1,fmin.GetLabel('i')), 'SIGI{0}(-)={1}'.format(i+1,fmin.GetLabel('sigi'))) )


  def DefineOutput(self):
    if not self.GetKey('anal'):
      for f in list(self.aver_data.values()) + self.anom_data_list:
        self.inp.AddFileToChild( f, self.outfilename['mtz'] )

  def TreatOutput(self):
    if not self.GetKey('anal'):
      self.SetArg('hklout', self.outfilename['mtz'])
