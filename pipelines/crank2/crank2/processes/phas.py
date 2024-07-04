#!/usr/bin/python
import os,sys
from process import process,crvapi
from program import program
import common
par=common.parameter


class phas(process):
  name="phasing and substructure refinement"
  short_name="phasing"
  supported_progs=["bp3","refmac"]
  supported_params = {}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True )
  supported_params['cycles'] = par( desc='Number of substructure refinement cycles', typ=int )
  supported_params['catch_output'] = par( desc='Catch program output continously', typ=bool )
  supported_params['occ_cut'] = par( desc='Occupancy threshold for atom deletion (<0 to disable)', typ=float )

  def Init(self):
    self.tables, self.in_table = [], 0

  def TreatInOutPar(self, set_all_par=False):
    if not self.GetProg(supported=True):
      if self.GetVirtPar('target') in ('SAD',):
        self.AddProg('refmac')
      else:
        self.AddProg('bp3')
    if self.GetParam('catch_output') is False:
      self.GetProg(supported=True).interact_output=False
      #self.Info('Continous log output disabled.')
    elif self.GetProg('refmac') or self.GetParam('catch_output') is True:
      self.GetProg(supported=True).interact_output=True
    # disabling hydrogens for refmac by default
    if self.GetProg(supported=True).nick=='refmac' and not self.GetProg(supported=True).IsKey('make'):
      self.GetProg(supported=True).SetKey('make',('hydrogens','no'))
    if self.GetParam('occ_cut') is None and self.inp.Get('model',filetype='res',typ='substr') and \
       (not self.parent_process or self.parent_process.nick not in ('refatompick','handdet')):
      self.SetParam('occ_cut',0.3)
    process.TreatInOutPar(self,set_all_par)

  def UpdateInfo(self,line,popen=None):
    if self.in_table:
      self.tables[-1].append(line)
      if line.strip()=='$$':
        self.in_table += 1
      if self.in_table==4:
        self.in_table = 0
        self.tables[-1] = ''.join(self.tables[-1])+'\n'
        self.PrintActualLogGraph()
    else:
      fom_table_start=self.GetProg(supported=True).GetStatGrep('lg_fom', from_str=line)
      if fom_table_start:
        self.tables.append([line,])
        self.in_table = 1

  def PrintActualLogGraph(self, finished=False):
    if finished:
      result = "Mean figure of merit (all reflections) is {0}".format(self.GetProg().GetStat('fom')[-1])
      self.report_fom = self.GetProg().GetStat('fom')[-1]
    if self.opened_loggraph:
      self.GetLogGraphHandle().seek(0,0)
      self.LGInfo(self.GetCCP4Header())
      self.LGInfo(self.GetLogGraphPrograms())
      for table in self.tables:
        self.LGInfo(table)
      if finished:
        self.LGInfo('\n$TEXT:Result: $$ Final result $$')
        self.LGInfo("{0}\n$$\n".format(result))
    if self.rv_report is not None:
      if finished:
        self.rv_report.Text('<BR><i><small>Result:</small></i>')
        button=crvapi.init_meta["help_btn_template"].replace('.html','.html#term-what-is-fom') if "help_btn_template" in crvapi.init_meta else ""
        self.rv_report.Text('&emsp;'+result.replace("merit ","merit "+button+" "))
        conclusion='&emsp;The value of FOM indicates the phases may be weaker but could be still sufficient for structure solution.'
        if self.report_fom<0.15:
          conclusion='&emsp;Judging from the value of FOM, the phases are weak or wrong.'
        if self.report_fom>0.3:
          conclusion='&emsp;The value of FOM indicates a reasonable quality of the phases.'
        if self.report_fom>0.4:
          conclusion='&emsp;The value of FOM indicates a good quality of the phases.'
        self.rv_report.Text(conclusion)
        if not self.GetCrankParent() or not self.GetCrankParent().ccp4i2:
          self.rv_report.Text("<BR><i><small>Output:</small></i>")
          self.rv_files = self.rv_report.DataFile(self.out.Get('model').GetFileName(), "xyz", "<small>Refined substructure and experimental map</small>")
          out_mtz = crvapi.SetMtzFWT(self,self.out.Get('mapcoef',typ='best'))
          self.rv_files.DataFile(out_mtz, "hkl:map", flush=True)
          #self.rv_files.DataFile(self.out.Get('mapcoef',typ='best',filetype='map').GetFileName('map'), "hkl:ccp4_map", flush=True)

  def RunBody(self,*args,**kwargs):
    self.GetProg().Run()
    self.Info( "Mean phasing figure of merit (all reflections) is {0}".format(self.GetProg().GetStat('fom')[-1]) )
    self.PrintActualLogGraph(finished=True)

  def RunPostprocess(self,restore=True,*args,**kwargs):
    if self.GetParam('occ_cut') is not None and self.GetParam('occ_cut')>=0.0:
      pdbcur=self.GetOrAddProg('pdbcur')
      pdbcur.inp.Set(self.out.Get('model',typ='substr'))
      occ_thr = self.GetParam('occ_cut')
      if pdbcur.inp.Get('model',filetype='pdb',typ='substr'):
        with open(pdbcur.inp.Get('model',filetype='pdb',typ='substr').GetFileName('pdb')) as f:
          occ_thr = occ_thr*max(pdbcur.GetStatGrep('occup',from_str=f.read()))
      pdbcur.SetKey('cutocc', occ_thr, keep_previous=False)
      pdbcur.Run()
    # a workaround for bp3 replacing the residue type (for example IOD always becomes I)
    # this can (and is better to) be removed after this is not done in bp3 or not needed anymore
    if self.GetProg().nick=='bp3':
      fixsub=self.AddProcess('fixsubstrpdb',propagate_inp=False)
      hand1=self.out.Get('model',typ='substr',filetype='pdb')
      hand2=self.out.Get('model',typ='substr',filetype='pdb',custom='otherhand')
      if hand2:
        fixsub.inp.Set(hand2)
        fixsub.SetRunDir(os.path.join(self.rundir,'fixsub-hand2'))
        fixsub.Run()
      fixsub.inp.Set(hand1)
      fixsub.SetRunDir(os.path.join(self.rundir,'fixsub-hand1'))
      fixsub.Run()
      self.processes.remove(fixsub)
    if not self.report_fom:
      self.result_str = "Wrong or no substructure: it's refinement resulted in FOM of 0. Structure solution unsuccessful."
      common.Error(self.result_str, nosuccess=True)
    process.RunPostprocess(self,restore,*args,**kwargs)
