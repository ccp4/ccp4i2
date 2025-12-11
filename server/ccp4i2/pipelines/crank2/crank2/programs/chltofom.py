#!/usr/bin/python
import os,sys
from program import program
import common

class chltofom(program):
  name="chltofom"
  binary="chltofom"
  modif='-'
  always_merge=True

  def Init(self):
    self.SetArg('stdin')

  def TreatInput(self):
    # we just take the first f + first hl/phfom object found
    # fo is optional
    #fsf=self.inp.Get('fsigf',filetype='mtz',col='f',typ='average')
    #if fsf:
    #  self.SetKey('colin-fo', fsf.GetFullLabel('f','sigf'))
    # HL -> PH+FOM
    self.hl=self.inp.Get('mapcoef',filetype='mtz',col=('hla','hlb'))#,inp_cont=fsf)
    self.phfom=self.inp.Get('mapcoef',filetype='mtz',col=('ph','fom'))#,inp_cont=fsf)
    if self.hl:
      self.SetKey('mtzin', self.hl.GetFileName('mtz'))
      self.SetKey('colin-hl', self.hl.GetFullLabel('hla','hlb','hlc','hld'))
    # PH+FOM -> HL
    elif self.phfom:
      self.SetKey('mtzin', self.phfom.GetFileName('mtz'))
      self.SetKey('colin-phifom', self.phfom.GetFullLabel('ph','fom'))
    if not self.hl and not self.phfom:
      common.Error('Neither HL coefficients nor phase+fom estimates in mtz format inputted to {0}'.format(self.name))

  def DefineOutput(self):
    outfilename=self.nick+'.mtz'
    if self.hl:
      self.inp.AddFileToChild(self.hl, outfilename, 'mtz')
      self.hl.SetLabel( ['ph','fom'] )
    else:
      self.inp.AddFileToChild(self.phfom, outfilename, 'mtz')
      self.phfom.SetLabel( ['hla','hlb','hlc','hld'] )

  def TreatOutput(self):
    if self.hl:
      self.SetKey( 'mtzout', self.hl.GetFileName('mtz') )
      self.SetKey( 'colout', self.hl.GetFullLabel('ph','fom') )
    else:
      self.SetKey( 'mtzout', self.phfom.GetFileName('mtz') )
      self.SetKey( 'colout', self.phfom.GetFullLabel('hla','hlb','hlc','hld') )
