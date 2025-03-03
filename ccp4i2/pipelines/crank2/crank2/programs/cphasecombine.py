#!/usr/bin/python
import os,sys
from program import program
import common

class cphasecombine(program):
  name="cphasecombine"
  binary="cphasecombine"
  modif='-'
  always_merge=True

  def Init(self):
    self.SetArg('stdin')

  def TreatInput(self):
    self.hls=self.inp.GetAll('mapcoef',filetype='mtz',col=('hla','hlb'))#,inp_cont=fsf)
    if len(self.hls)<2:
      common.Error('{0} requires two sets of HL coefficients; {1} inputted.'.format(self.name,len(self.hls)))
    self.SetKey('mtzin', self.hls[0].GetFileName('mtz'))
    self.SetKey('colin-hl-1', self.hls[0].GetFullLabel('hla','hlb','hlc','hld'))
    self.SetKey('colin-hl-2', self.hls[1].GetFullLabel('hla','hlb','hlc','hld'))

  def DefineOutput(self):
    outfilename=self.nick+'.mtz'
    self.hl_comb=self.out.AddNew('mapcoef', outfilename, typ='combined', filetype='mtz', xname=self.hls[0].xname)
    self.hl_comb.SetLabel( ['hla','hlb','hlc','hld'], bad_lbl_obj=(self.hls[0],self.hls[1]) )

  def TreatOutput(self):
    self.SetKey( 'mtzout', self.hl_comb.GetFileName('mtz') )
    self.SetKey( 'colout', self.hl_comb.GetFullLabel('hla','hlb','hlc','hld') )
