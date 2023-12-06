#!/usr/bin/python
import os,sys
from process import process
import common

class createfree(process):
  name="definition of a free set of reflections"
  short_name="definition of free set"
  supported_progs=["sftools"]
  supported_params = {}
  supported_params['fraction'] = common.parameter( desc='fraction of the free reflections', typ=float )
  supported_params['freetype'] = common.parameter( desc='type of the free reflections set', typ=str )


  def TreatInOutPar(self, set_all_par=False):
    self.sft=self.GetOrAddProg('sftools')
    if self.GetParam('fraction') is None:
      self.SetParam('fraction', 0.05)
    if self.GetParam('freetype') is None:
      self.SetParam('freetype', 'freeR')
    process.TreatInOutPar(self,set_all_par)

  def RunBody(self,*args,**kwargs):
    mtzobj = self.inp.GetAll(filetype='mtz')
    if not mtzobj:
      common.Error('Cannot create free set - no mtz file inputted.')
    allres = [o.GetResolution(self,accept_none=True) for o in mtzobj]
    highres = min([r if r else 100. for r in allres])
    mtzobj = [mtzobj[i] for i,r in enumerate(allres) if r==highres]
    allreslow = [o.GetResolution(self,accept_none=True,lowres=True) for o in mtzobj]
    mtzobj = mtzobj[ allreslow.index(max(allreslow))  if allreslow else 0]
    exclude_out=self.sft.out.AddNew('exclude', self.sft.nick+'_free.mtz', typ=self.GetParam('freetype'))
    exclude_out.SetLabel('free')
    backup=self.sft.BackupAnyPars()
    self.sft.SetKey( 'read', '"'+mtzobj.GetFileName('mtz')+'"' )
    # increase the fraction if too small - for bias est. only as of now
    if not self.IsInputtedParam('fraction') and self.GetParam('freetype'):
      self.sft.Run(clear_out=False)
      num_ref = self.sft.GetStat('num_refl', accept_none=True)
      if num_ref and int(num_ref)<15000:
        self.SetParam('fraction', 0.08)
      if num_ref and int(num_ref)<10000:
        self.SetParam('fraction', 0.1)
    self.sft.SetKey( 'calc', ('I', 'col', exclude_out.GetLabel('free'), '=', 'rfree({0})'.format(self.GetParam('fraction'))) )
    self.sft.SetKey( 'write', (exclude_out.GetFile('mtz').name, '\nY') )
    self.sft.SetKey( 'quit', '\nY' )
    self.sft.Run(clear_out=False)
    self.sft.RestoreAnyPars(backup)

  def RunPostprocess(self,restore=True,*args,**kwargs):
    if self.rv_report:
      self.rv_report.Text("<i><small>Result:</small></i>")
      self.rv_report.Text('{}% of reflections have been flagged as free for cross-validation.'.format(self.GetParam('fraction')*100), flush=True)
    process.RunPostprocess(self,restore,*args,**kwargs)


