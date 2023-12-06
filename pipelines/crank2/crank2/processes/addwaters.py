#!/usr/bin/python
import os,sys
from process import process
import common

class addwaters(process):
  name="adding waters into electron density map"
  short_name="adding waters"
  supported_progs=["coot"]


  def TreatInOutPar(self, set_all_par=False):
    self.coot=self.GetOrAddProg('coot')
    self.model=self.inp.Get('model',typ=('partial','partial+substr'))
    if not self.model:
      common.Error('Model required for {}'.format(self.name))
    mapc=self.inp.Get('mapcoef',typ=('weighted','combined'),col=('f','ph'))
    if not mapc:
      common.Error('Map coefitients required for {}'.format(self.name))
    with open(self.coot.inpfilename['py'],'w') as f:
      f.write('modl=read_pdb(r\'{}\')\n'.format(self.model.GetFileName()))
      f.write('map=make_and_draw_map(r\'{}\', \'{}\', \'{}\', \'{}\', 0, 0)\n'.format(mapc.GetFileName(),mapc.GetLabel('f'),mapc.GetLabel('ph'),mapc.GetLabel('ph')))
      f.write('set_ligand_water_to_protein_distance_limits(2.2, 3.3)\n')
      f.write('execute_find_waters_real(map, modl, 0, 1.3)\n')
      f.write('write_pdb_file(modl, \'{}\')\n'.format(self.coot.outfilename['pdb']))
      f.write('coot_real_exit(0)')
    process.TreatInOutPar(self,set_all_par)

  def RunBody(self,*args,**kwargs):
    mtzobj = self.inp.GetAll(filetype='mtz')
    self.coot.Run()
    self.outmodel = self.out.AddCopy(self.model)
    self.out.SetFileToChild(self.outmodel, self.coot.outfilename['pdb'], 'pdb')
