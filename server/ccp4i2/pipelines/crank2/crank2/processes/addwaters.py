#!/usr/bin/python
import os,sys,string
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
    self.mapc=self.inp.Get('mapcoef',typ=('weighted','combined'),col=('f','ph'))
    if not self.mapc:
      common.Error('Map coefitients required for {}'.format(self.name))
    process.TreatInOutPar(self,set_all_par)

  def RunPreprocess(self,*args,**kwargs):
    process.RunPreprocess(self,*args,**kwargs)
    # coot seems to add waters as chain Z also if another Z already exists...
    pdbcur=self.AddProg('pdbcur')
    pdbcur.SetKey('summarise')
    pdbcur.inp.Set(self.model)
    pdbcur.Run()
    chain_ids = pdbcur.GetStat('chain_ids')
    if 'Z' in chain_ids:
      for new_id in string.printable:
        if new_id not in chain_ids and new_id not in pdbcur.bad_chain_ids:
          pdbcur.SetKey('renchain', ('\"Z\"', "\""+new_id+"\"") )
          break
      if pdbcur.IsKey('renchain'):
        pdbcur.Run()
        self.model=self.inp.Add(pdbcur.out.Get('model'))
    # create the coot script
    if not os.path.isdir('coot'):
      os.mkdir('coot')
    with open(os.path.join('coot',self.coot.inpfilename['py']),'w') as f:
      f.write('modl=read_pdb(r\'{}\')\n'.format(os.path.join(self.model.GetFileName())))
      f.write('map=make_and_draw_map(r\'{}\', \'{}\', \'{}\', \'{}\', 0, 0)\n'.format(self.mapc.GetFileName(),self.mapc.GetLabel('f'),self.mapc.GetLabel('ph'),self.mapc.GetLabel('ph')))
      f.write('set_ligand_water_to_protein_distance_limits(2.2, 3.3)\n')
      f.write('execute_find_waters_real(map, modl, 0, 1.3)\n')
      f.write('write_pdb_file(modl, \'{}\')\n'.format(self.coot.outfilename['pdb']))
      f.write('coot_real_exit(0)')


  def RunBody(self,*args,**kwargs):
    self.coot.Run()
