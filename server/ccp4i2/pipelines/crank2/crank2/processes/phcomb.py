#!/usr/bin/python
import os,sys
from process import process
import common
par=common.parameter


class phcomb(process):
  name="reciprocal space phase combination"
  short_name="phase combination"
  supported_progs=["multicomb","refmac"]
  supported_params = {}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True )
  supported_params['cycles'] = par( desc='Number of substructure refinement cycles', typ=int )
  supported_params['beta'] = par( desc='Bias parameter for the map combination', typ=float )
  supported_params['comb_with_typ'] = par( desc='Combine phases with the specified mapcoef type', typ=str )

  def TreatInOutPar(self, set_all_par=False):
    # consult with Raj, adjust accordingly if needed
    if not self.GetProg(supported=True):
      if self.GetVirtPar('target') == 'MLHL':
        self.AddProg('multicomb')
      else:
        self.AddProg('refmac')
    if not self.GetVirtPar('comb_with_typ'):
      self.SetVirtPar('comb_with_typ', 'densmod')
    # disabling hydrogens for refmac by default
    if self.GetProg(supported=True).nick=='refmac' and not self.GetProg(supported=True).IsKey('make'):
      self.GetProg(supported=True).SetKey('make',('hydrogens','no'))
    process.TreatInOutPar(self,set_all_par)

