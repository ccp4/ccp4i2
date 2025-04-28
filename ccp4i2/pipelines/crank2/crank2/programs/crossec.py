#!/usr/bin/python

from .. import common
from ..program import program


class crossec(program):
  name="CROSSEC"
  binary="crossec"
  ccp4_parsing=True
  never_merge=True
  stat={}
  stat['anom_scat_coef'] = common.stats(name='anomalous scattering coefficients (f\', f\'\')', 
                                        regexp=r"{0}\s+\d+\.\d+\s+(\S+)\s+(\S+)", multiple=True)


  def TreatInput(self):
    self.dnames, self.wavels = [], []
    # only doing the estimation for the "main" atom type (could be easily extended)
    if not self.GetKey('ATOM'):
      self.model = self.inp.Get('model',has_atomtypes=True)
      if self.model:
        self.SetKey('ATOM', self.model.GetAtomType())
      else:
        common.Error('No model with atom type inputted for {0}'.format(self.name))
    if not self.GetKey('NWAV') and not self.GetKey('CWAV'):
      self.dnames = list(set( [f.GetDataName() for f in self.inp.GetAll('fsigf',xname=self.model.GetCrystalName())] ))
      wavels = [(dn,self.inp.Get('fsigf',xname=self.model.GetCrystalName(),dname=dn).GetWavelength(self,accept_none=True)) for dn in self.dnames]
      self.wavels = sorted( [(dn,wl) for dn,wl in wavels if wl], key=lambda w: w[1] )
      if self.wavels:
        self.SetKey( 'NWAV', [len(self.wavels)]+list(list(zip(*self.wavels))[1]) )
      else:
        common.Error('No or zero wavelength inputted for {0}'.format(self.name))

