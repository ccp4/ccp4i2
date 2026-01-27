#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class mapro(program):
  name="MAPRO"
  binary="mapro"
  modif='-'
  stat={}
  stat['cld'] = common.stats(name='correlation of map with its local deviation', regexp=r"Correlation of the map with its local deviation is\s+(\S+)")
  stat['skewness'] = common.stats(name='skewness of map', regexp=r"Skewness of the map is\s+(\S+)")


  def TreatInput(self):
    # input map
    mapc=self.inp.Get('mapcoef',filetype='map')
    if mapc:
      self.SetArg('mapin',mapc.GetFileName('map'))
    else:
      common.Error('No map inputted for {0}'.format(self.name))
    model=self.inp.Get('model',filetype='pdb')
    if model:
      self.SetArg('pdbin',model.GetFileName('pdb'))
