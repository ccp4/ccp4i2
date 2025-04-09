#!/usr/bin/python

from .. import common
from ..program import program


class invfft(program):
  name="invfft"
  binary="cinvfft"
  modif='-'
  always_merge=True

  def Init(self):
    self.SetArg('stdin')
    self.outfilename = { 'mtz': self.nick+'.mtz', }

  def TreatInput(self):
    # use the first map found in inp.mapcoef
    self.mapcin=self.inp.Get('mapcoef',filetype='map')
    if self.mapcin:
      self.SetKey('mapin', self.mapcin.GetFileName('map'))
    else:
      common.Error('No map inputted to {0}'.format(self.name))
    # fo is needed too
    mapc=self.inp.Get('fsigf',filetype='mtz',col='f')
    if mapc:
      self.SetKey('mtzin',mapc.GetFileName('mtz'))
    else:
      common.Error('No observed data inputted to {0}'.format(self.name))

  def DefineOutput(self):
    out_mapc=self.out.Add(self.mapcin)
    self.out.AddFileToChild(out_mapc,self.outfilename['mtz'],'mtz')
    out_mapc.SetLabel( ['f','ph'] )

  def TreatOutput(self):
    self.SetKey( 'mtzout', self.out.mapcoef[-1].GetFileName('mtz') )
    self.SetKey( 'colout', self.out.mapcoef[-1].GetFullLabel('f','ph') )
