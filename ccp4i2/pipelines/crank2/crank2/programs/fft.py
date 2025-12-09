#!/usr/bin/python
import os,sys
from program import program
import common

class fft(program):
  name="cfft"
  binary="cfft"
  modif='-'
  always_merge=True

  def Init(self):
    self.SetArg('stdin')
    self.outfilename = { 'map': self.nick+'.map' }

  def TreatInput(self):
    # phase information
    for mapc in self.inp.GetAll('mapcoef'):
      if mapc.GetFile('mtz') and (mapc.GetLabel('hla') or mapc.GetLabel('ph')):
        self.SetKey('mtzin',mapc.GetFileName('mtz'))
        if mapc.GetLabel('f') and mapc.GetLabel('ph'):
          self.SetKey('colin-fc',mapc.GetFullLabel('f','ph'))
        else:
          self.SetKey('colin-hl',mapc.GetFullLabel('hla','hlb','hlc','hld'))
        # saving the object so that we can easily access it by the output
        self.mapcin = mapc
        break
    else:
      common.Error('No phase information inputted for {0}'.format(self.name))
    # Fo
    self.fsf=self.inp.Get('fsigf',typ=('average','delta-anom'),filetype='mtz')
    if self.fsf:
      self.SetKey('colin-fo', self.fsf.GetFullLabel('f','sigf'))
    elif not self.GetKey('colin-fc'):
      common.Error('No data inputted for {0}'.format(self.name))

  def DefineOutput(self):
    out_mapc=self.out.Add(self.mapcin)
    self.out.AddFileToChild(out_mapc,self.outfilename['map'],'map')

  def TreatOutput(self):
    self.SetKey( 'mapout', self.out.mapcoef[-1].GetFileName('map') )
