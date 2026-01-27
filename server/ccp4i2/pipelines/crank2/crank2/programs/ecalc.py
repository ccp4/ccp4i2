#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class ecalc(program):
  name="ECALC"
  binary="ecalc"
  never_merge=True
  ccp4_parsing=True

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    self.inp_cont = self.inp.Get('fsigf',typ='fa',col='f',filetype='mtz')
    if not self.inp_cont:
      self.inp_cont = self.inp.Get('fsigf',typ='delta-anom',col='f',filetype='mtz')
    if not self.inp_cont:
      self.inp_cont =  self.inp.Get('fsigf',typ='average',col='f',filetype='mtz')
    if not self.inp_cont:
      common.Error("No suitable file inputted for E-values generation by {0}".format(self.name))
    self.SetArg('hklin', self.inp_cont.GetFileName('mtz'))
    if self.inp_cont.GetType()=='fa' or self.inp_cont.GetType()=='delta-anom':
      self.AddToKey('labin', ('DPH='+self.inp_cont.GetLabel('f'), 'SIGDPH='+self.inp_cont.GetLabel('sigf')))
    # we'll need to add eaver or adjust the data model to support this
    elif self.inp_cont.GetType()=='average':
      self.AddToKey('labin', ('FP='+self.inp_cont.GetLabel('f'), 'SIGFP='+self.inp_cont.GetLabel('sigf')))
#      self.AddToKey('labin', 'FP='+self.inp_cont.GetLabel('f'), 'SIGFP='+self.inp_cont.GetLabel('sigf'))

  def DefineOutput(self):
    self.outmtz=self.out.AddCopy(self.inp_cont)
    self.out.AddFileToChild(self.outmtz, self.outfilename['mtz'], 'mtz')
    self.outmtz.SetLabel(('e','sige'))

  def TreatOutput(self):
    self.SetArg('hklout', self.out.Get().GetFileName('mtz'))
    self.AddToKey( 'labout', ('E='+self.outmtz.GetLabel('e'),'SIGE='+self.outmtz.GetLabel('sige')) )
