#!/usr/bin/python
import os,sys
from program import program
import common

class solomon(program):
  name="Solomon"
  binary="solomon"
  labelout_prefix="SOLO_"
  stat={}
  stat['contrast_inv'] = common.stats(name='inverse map contrast', regexp=r"ratio:\s+(\S+)")
  references = ( "Abrahams JP and Leslie AGW (1996) Methods used in the " +
                 "structure determination of bovine mitochondrial F1 ATPase. "+
                 "Acta Cryst. D52, 30-42.", )

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    # input map
    inp_mtz=self.inp.Get('fsigf',typ='average',filetype='mtz')
    if inp_mtz: inp_mtz.GetFileName()  # register mtz has been asked 
    self.mapc=self.inp.Get('mapcoef',filetype='map',inp_cont=inp_mtz)
    if self.mapc:
      self.SetArg('mapin',self.mapc.GetFileName('map'))
    else:
      common.Error('No map inputted for {0}'.format(self.name))
    # solvent content - if supplied with an input container
    solv_obj = self.inp.Get(has_solvent_content=True)
    if solv_obj and not self.IsKey('slvfrc'):
      self.SetKey('slvfrc', solv_obj.solvent_content)
    if self.process.IsNonFalseVirtPar('solvent_perturb') and self.GetKey('slvfrc'):
      self.SetKey('slvfrc', self.GetKey('slvfrc') + self.process.GetVirtPar('solvent_perturb'), keep_previous=False )
    if not self.GetKey('slvfrc'):
      common.Error('No solvent content inputted for {0}'.format(self.nick))

  def TreatParams(self):
    if self.process.GetVirtPar('solvent_content'):
      self.SetKey('slvfrc', self.process.GetVirtPar('solvent_content'))
    program.TreatParams(self)

  def DefineOutput(self):
    out_mapc = self.out.AddNew( 'mapcoef', self.nick+'.map', 'map' )
    out_mapc.SetType('densmod')
    out_mapc.SetCrystalName(self.mapc.GetCrystalName())
    out_mapc.SetDataName(self.mapc.GetDataName())

  def TreatOutput(self):
    self.SetArg('mapout',self.out.mapcoef[-1].GetFileName('map'))
    self.SetKey('mapout')
