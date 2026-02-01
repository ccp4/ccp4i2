#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class seqwt(program):
  name="SEQWT"
  binary="seqwt"
  ccp4_parsing=True
  stat={}
  stat['mol_weight'] = common.stats(name='total molecular weight in Da', regexp=r"Total\s+Molecular\s+Weight\s*(\d+)\s+\(Da\)")

  def TreatInput(self):
    # input sequence - there might be problems with fas (?)
    seq=self.inp.Get('sequence',filetype=('seq','pir','fas'))
    if seq and seq.GetFileName():
      self.SetArg('sequence',seq.GetFileName())
    else:
      common.Error('Sequence in neither pir nor seq inputted to {0}'.format(self.name))

