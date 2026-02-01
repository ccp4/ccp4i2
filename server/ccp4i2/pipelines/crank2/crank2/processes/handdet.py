#!/usr/bin/python
import os,sys
from ..process import process, crvapi
from .. import common
par=common.parameter


class handdet(process):
  name="hand determination"
  short_name=name
  supported_progs=["mapro"]
  supported_procs=["phas","dmfull","refatompick"]
  supported_params = {}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True )
  supported_params['solvent_content'] = par( desc='(Expected) solvent fraction of crystal', typ=float, share=True )
  supported_params['ignore_passed_phas'] = par( desc='Ignore passed phasing (if any)', typ=bool )
  supported_params['threshold_discrim'] = par( desc='Consider the hand unresolved if the difference in score is smaller than this value', typ=float )
  # there may be previous phas/dmfull job passed to hand determination
  phas,phas2=None,None
  dmfull=None

  # enant. pairs: P41,P43 P4122,P4322 P41212,P43212 P31,P32  P3112,P3212 P3121,P3221  P61,P65
  enant_pairs = ( (76,78),  (91,95),     (92,96),  (144,145), (151,153),  (152,154), (169,170), \
                  (171,172), (178,179),  (180,181),  (212,213)    )
  #                P62,P64  P6122,P6522 P6222,P6422 P4332,P4132

  def TreatInOutPar(self, set_all_par=False):
    self.dmf = [None,None]
    self.dmf[0] = self.GetProcess('dmfull')
    if not self.dmf[0]:
      for i in range(2):
        self.dmf[i] = self.AddProcess('dmfull')
        self.dmf[i].AddProcess('dm').AddProg('solomon')
        self.dmf[i].AddProcess('phcomb').AddProg('multicomb')
        self.dmf[i].TreatInOutPar(set_all_par)
        if not self.dmf[i].GetProcess('phcomb').GetVirtPar('target'):
          self.dmf[i].GetProcess('phcomb').SetVirtPar('target', 'MLHL')
    else:
      if not self.dmf[0].GetVirtPar('target'):
        self.dmf[0].SetVirtPar('target', 'MLHL')
      self.dmf[0].TreatInOutPar(set_all_par)
      self.dmf[1] = self.AddProcessCopy(self.dmf[0], deeper_copy=True)
      self.dmf[1].TreatInOutPar(set_all_par)
    self.dmf[0].SetParam('threshold_stop', 1.), self.dmf[1].SetParam('threshold_stop', 1.)
    self.dmf[0].SetParam('solvent_perturb', 0.), self.dmf[1].SetParam('solvent_perturb', 0.)
    # if hand1 phasing has not been passed then create a new one
    if not self.phas:
      if not self.GetParam('ignore_passed_phas'):
        self.phas = self.GetProcess(('phas','refatompick'))
      if not self.phas:
        self.phas = self.AddProcess('phas')
      if self.GetVirtPar('target') and not self.phas.GetVirtPar('target'):
        self.phas.SetVirtPar('target', self.GetVirtPar('target'))
      if not self.phas.GetProg():
        self.phas.TreatInOutPar(set_all_par)
    if not self.phas2:
      if self.phas.nick=='phas':
        self.phas2 = self.AddProcessCopy(self.phas, deeper_copy=True)
      elif self.phas.GetProcess('phas'):
        self.phas2 = self.AddProcessCopy(self.phas.GetProcess('phas'), deeper_copy=True)
      if self.phas.out.Get('model',custom='otherhand') and self.phas.out.Get('mapcoef',custom='otherhand'):
        self.phas2.out.Set(self.phas2.out.GetAll('model',custom='otherhand',stored_order=True))
        self.phas2.out.Set(self.phas2.out.GetAll('mapcoef',custom='otherhand',stored_order=True))
      #else:
        #self.phas2 = self.AddProcess('phas')
        #if self.phas.GetParam('target'):
        #  self.phas2.SetParam('target', self.phas.GetParam('target'))
        #self.phas2.AddProg(self.phas.GetProg().nick)
    self.mapro = self.GetOrAddProg('mapro')
    if not self.IsParam('threshold_discrim'): #and self.GetParam('target')=='SAD':  # only for SAD now, add other soon
      self.SetParam('threshold_discrim', 12.)
    self.guess=None
    process.TreatInOutPar(self,set_all_par)

  def RunBody(self,*args,**kwargs):
    self.pdbcur = self.GetOrAddProg("pdbcur")
    self.cld, self.fomcontr, self.score = [], [], []
    changedmtz=False
    # now run the actual hand determination processes
    self.ProcessPhasHand1()
    if not self.phas2.out.Get('mapcoef',typ='best',custom='otherhand'):
      # next step may be skipped for bp3 as it can do the other hand model
      self.CreateHand2Model(firsthand=self.phas.out.Get('model',typ='substr'))
    changedmtz=self.ChangeHand2Spacegroup(self)
    self.ProcessPhasHand2(changedmtz)
    self.PrintActualLogGraph()
    self.RunDMHand(0,self.phas)
    self.RunDMHand(1,self.phas2,changedmtz)
    self.Evaluate()

  def ProcessPhasHand1(self):
    if not self.phas.out.Get('mapcoef',typ='best'):
      self.phas.Run()
      self.phas.fom = self.phas.GetProg(supported=True).GetStat('fom')[-1]
    self.GetCLD(1, self.phas)

  def ProcessPhasHand2(self,changedmtz):
    if changedmtz:
      self.phas2.inp.Set(changedmtz)
    if not self.phas2.out.Get('mapcoef',typ='best',custom='otherhand'):
      # skipped for bp3
      self.phas2.inp.Set(self.out.Get('model',custom='otherhand'))
      if self.phas2.GetProg().nick=='refmac' and self.phas.out.Get('datafile',typ='dluz'):
        self.phas2.inp.Set(self.phas.out.Get('datafile',typ='dluz'))
        if not self.phas2.IsParam('cycles'):
          self.phas2.SetParam('cycles',1)
      self.phas2.Run(rundir=os.path.join(self.rundir,'phas2'))
      self.phas2.out.model[-1].custom.append('otherhand')
      self.phas2.fom = self.phas2.GetProg().GetStat('fom')[-1]
    self.GetCLD(2, self.phas2)

  def GetCLD(self,hand,phas):
    if not self.mapro.IsArg('cld'):
      self.mapro.SetArg('cld')
      self.mapro.SetArg('skew')
    self.mapro.inp.SetCopy(phas.out.Get('mapcoef',typ='best'))
    self.mapro.inp.SetCopy(phas.out.Get('model',typ='substr'))
    self.mapro.inp.Clear('fsigf')
    self.mapro.runname='mapro_hand{0}'.format(hand)
    self.mapro.Run()
    self.cld.append(round(self.mapro.GetStat('cld'),4))

  def RunDMHand(self,hand,phas,changedmtz=None):
    self.dmf[hand].inp.Set(phas.out.Get('mapcoef',col='hla',typ='best'))
    if changedmtz:
      self.dmf[hand].inp.Set(changedmtz)
    if hand==1:
      self.dmf[hand].inp.Set(self.out.Get('model',custom='otherhand'))
    self.dmf[hand].Run(rundir=os.path.join(self.rundir,'dmfull_hand'+str(hand+1)))
    self.fomcontr.append( [b/a for a,b in zip(self.dmf[hand].dm.contr_inv_all,self.dmf[hand].ph.fom_all) if a>0.] )

  def CreateHand2Model(self,firsthand,spgr_num=None):
    self.inversion_operation = '-X,-Y,-Z'
    if spgr_num is None:
      spgr_num = self.phas.out.Get('mapcoef',typ='best').GetSpacegroupNumber(self)
    if spgr_num==80:
      self.inversion_operation = '-X+1/2,-Y,-Z'
    elif spgr_num==98:
      self.inversion_operation = '-X+1/2,-Y,-Z+1/4'
    elif spgr_num==210:
      self.inversion_operation = '-X+1/4,-Y+1/4,-Z+1/4'
    self.pdbcur.inp.Set(firsthand, propagate=False)
    self.pdbcur.SetKey('symop', self.inversion_operation)
    self.pdbcur.SetKey('symcommit')
    self.pdbcur.outfilename['pdb'] = os.path.splitext(firsthand.GetFile('pdb').name)[0]+'.oh.pdb'
    self.pdbcur.Run()
    self.pdbcur.out.model[-1].custom.append('otherhand')
    self.programs.remove(self.pdbcur)

  @staticmethod
  def ChangeHand2Spacegroup(self,spgr_num=None,spgr_num2=None,fsigfs=[]):
    outputmtz=''
    if not spgr_num:
      firsthandmtz=self.phas.out.Get('mapcoef',typ='best')
      spgr_num = firsthandmtz.GetSpacegroupNumber(self)
      #outputmtz=os.path.splitext(firsthandmtz.GetFileName())[0]+'.oh.mtz'
    if not spgr_num2:
      is_enant = [list(ep) for ep in self.enant_pairs if spgr_num in ep]
      if is_enant:
        is_enant[0].remove(spgr_num)
        spgr_num2 = is_enant[0][0]
    if spgr_num and spgr_num2:
      if not fsigfs:
        fsigfs = self.phas.inp.GetAll('fsigf',filetype='mtz')
      sft = self.AddProg('sftools',propagate_out=False)
      fsigf_lst,fsigf_lst2=[],[]
      for i,fsigf in enumerate(fsigfs):
        fname=fsigf.GetFileName('mtz')
        if fname and fname not in fsigf_lst:
          fsigf_lst.append(fname)
          fsigf_lst2=[fsf for fsf in fsigfs if fsf.GetFileName('mtz') and fsf.GetFileName('mtz')==fname]
          labels=[]
          for fsf in fsigf_lst2:
            labels.extend(fsf.GetAllLabels(labels_only=True))
          labels=list(set(labels))
          outputmtz=os.path.join(self.rundir,os.path.basename(os.path.splitext(fname)[0]+'.oh.mtz'))
          sft.runname='sftools_setspgr_'+str(i)
          sft.SetRunDir()
          sft.SetKey( 'read', ('"'+fsigf.GetFileName('mtz')+'"', 'col', '"'+'" "'.join(labels)+'"') )
          sft.SetKey('set spacegroup\n',spgr_num2)
          sft.SetKey('reduce','\n')
          sft.SetKey('merge','\n')
          sft.SetKey('write', (outputmtz,'\nY'))
          sft.SetKey('quit\nY')
          sft.Run(clear_out=False)
          sft.ClearAnyParams()
          for fsigf2 in fsigf_lst2:
            fsigf_new=sft.out.AddCopy(fsigf2)
            fsigf_new.SetFile(outputmtz,filetype='mtz')
            fsigf_new.spgr_num=spgr_num2
            fsigf_new.spgr=None
            fsigf_new.custom.append('otherhand')
      self.programs.remove(sft)
      return sft.out.fsigf
    else:
      return None

  def Evaluate(self):
    if len(self.cld)<2 or len(self.fomcontr)<2:
      common.Error('Hand cannot be evaluated - not all information was acquired.')
    # evaluate the score
    weight=0.95*max(1,len(self.fomcontr[0]))
    # it seems that cld is much less reliable for SIRAS?  needs to be better analyzed...
    if self.GetParam('target')!='SAD':
      weight*=0.5
    for hand in (0,1):
      otherhand=(hand+1)%2
      #self.score.append( weight*int(self.cld[hand]>self.cld[otherhand]) + \
      #                   len([1 for h,oh in zip(self.fomcontr[hand],self.fomcontr[otherhand]) if h>oh]) )
      self.score.append( weight*min(1.,abs(self.cld[hand]-self.cld[otherhand]+max(self.cld[otherhand]/3.,0.))/0.02)*int(self.cld[hand]>self.cld[otherhand]) + \
                         len([1 for h,oh in zip(self.dmf[hand].ph.fom_all,self.dmf[otherhand].ph.fom_all) if h>oh]) )
    # decide the hand
    self.spacegroup_change=False
    self.chosen_hand=0
    if not self.IsParam('threshold_discrim') or abs(self.score[1]-self.score[0])>=self.GetParam('threshold_discrim'):
      if self.score[1]>self.score[0]:
        self.chosen_hand=2
        self.chosen_phas=self.phas2
        self.other_phas=self.phas
        if self.chosen_phas.out.Get('mapcoef').GetSpacegroupNumber(self)!=self.other_phas.out.Get('mapcoef').GetSpacegroupNumber(self):
          self.spacegroup_change=self.chosen_phas.out.Get('mapcoef').GetSpacegroup(self)
      else:
        self.chosen_hand=1
        self.chosen_phas=self.phas
        self.other_phas=self.phas2
    else:
      # if the discrimination threshold is not passed then we keep the hands order as is, keep "chosen_hand" at 0 and tell the next steps to run both hands
      self.chosen_phas=self.phas
      self.other_phas=self.phas2
      both_hand_procs = ('dmfull','comb_phdmmb','mbref')
      if self.GetCrankParent() and not self.parent_process.nick=='substrdet':
        for pr in both_hand_procs:
          if self.GetCrankParent().prep.GetProcess(pr):
            if not self.GetCrankParent().prep.GetProcess(pr).IsParam('handdet'):
              self.GetCrankParent().prep.GetProcess(pr).SetParam('handdet')
    # set the output to chosen hand
    self.out.Set(self.chosen_phas.out.Get('model',typ='substr'))
    self.out.Set(self.chosen_phas.out.GetAll('mapcoef',stored_order=True))
    self.out.Set(self.chosen_phas.inp.GetAll('fsigf',stored_order=True))
    # Try standalone import first (original behavior), fall back to package import (ccp4i2 context)
    try:
      import inout
    except ImportError:
      from crank2 import inout
    self.out2 = inout.input_output(is_output=True,parent=self)
    self.out2.Set(self.other_phas.out.Get('model',typ='substr'),propagate=False)
    self.out2.Set(self.other_phas.out.GetAll('mapcoef',stored_order=True),propagate=False)
    self.out2.Set(self.other_phas.inp.GetAll('fsigf',stored_order=True),propagate=False)
    self.guess=1
    if abs(self.score[1]-self.score[0])<=8:
      self.guess=0
    elif abs(self.score[1]-self.score[0])>=13:
      self.guess=2

  def PrintActualLogGraph(self, result1=None, result2=None):
    dm = hasattr(self.dmf[1],'ph') and hasattr(self.dmf[1].ph,'fom_all') and \
         self.dmf[0].GetParam('dmcyc')>0 and self.dmf[1].GetParam('dmcyc')>0
    if self.opened_loggraph:
      self.GetLogGraphHandle().seek(0,0)
      self.LGInfo(self.GetCCP4Header())
      self.LGInfo(self.GetLogGraphPrograms())
      if self.cld[0] and self.cld[1]:
        self.LGInfo('\n$TEXT: Phasing CLD result: $$ Phasing CLD result $$')
        self.LGInfo('Phasing CLD for the first hand is {0} and for the other hand {1}\n$$'.format(round(self.cld[0]*100,2),round(self.cld[1]*100,2)))
      if dm:
        self.LGInfo('\n $TABLE : Fom and contrast of each hand per cycle:')
        self.LGInfo('$GRAPHS    : FOM for each hand per cycle :A:1,2,3:')
        self.LGInfo('           : Contrast for each hand per cycle :A:1,4,5:\n$$')
        self.LGInfo('Cycle FOM_Hand1 FOM_Hand2 Contrast_Hand1 Contrast_Hand2 $$ $$')
        for i,n in enumerate(zip(self.dmf[0].ph.fom_all,self.dmf[1].ph.fom_all,self.dmf[0].dm.contr_inv_all,self.dmf[1].dm.contr_inv_all)):
          self.LGInfo('{0:3}  {1:8.3f}  {2:8.3f}   {3:8.3f}   {4:8.3f}'.format(i+1,n[0],n[1],n[2],n[3]))
        self.LGInfo('$$\n')
      if result1:
        self.LGInfo('\n$TEXT:Result: $$ Final result $$')
        self.LGInfo(result1), self.LGInfo(result2+'\n$$')
    if self.rv_report is not None:
      if self.cld[0] and self.cld[1] and not dm:
        self.rv_report.Text("<BR><i><small>Phasing CLD result:</small></i>")
        button=crvapi.init_meta["help_btn_template"].replace('.html','.html#term-what-is-the-cld-value-reported-in-hand-determination') if "help_btn_template" in crvapi.init_meta else ""
        self.rv_report.Text('&emsp;Phasing CLD'+button+' for the first hand is {0} and for the other hand {1}'.format(round(self.cld[0]*100,2),round(self.cld[1]*100,2)), flush=True)
      if dm:
        if not hasattr(self,'rv_plot_fom'):
          cyc = [i+1 for i in range(len(self.dmf[0].ph.fom_all))]
          self.rv_plot_fom = self.rv_report.Plot( 'FOM for each hand vs cycle', "Cycle", "Figure of merit", \
                   block="FOM and contrast of hands vs DM cycle", legendloc='nw', intx=True)
          self.rv_plot_con = self.rv_plot_fom.parent.Plot( 'Contrast for each hand vs cycle', "Cycle", "Inverse contrast", intx=True )
          self.rv_fom1 = self.rv_plot_fom.PlotLine( ["x",cyc], ["hand 1",self.dmf[0].ph.fom_all] )
          self.rv_fom2 = self.rv_plot_fom.PlotLine( ["x",cyc], ["hand 2",self.dmf[1].ph.fom_all] )
          self.rv_con1 = self.rv_plot_con.PlotLine( ["x",cyc], ["hand 1",self.dmf[0].dm.contr_inv_all] )
          self.rv_con2 = self.rv_plot_con.PlotLine( ["x",cyc], ["hand 2",self.dmf[1].dm.contr_inv_all], flush=True )
      if result1:
        self.rv_report.Text("<i><small>Result:</small></i>")
        button=crvapi.init_meta["help_btn_template"].replace('.html','.html#term-what-is-fom') if "help_btn_template" in crvapi.init_meta else ""
        self.rv_report.Text('&emsp;'+result1.replace('FOM','FOM'+button))
        if self.chosen_hand:
          if abs(self.score[0]-self.score[1])<15:
            result2+=" There is a small chance that the chosen hand is not correct."
          else:
            result2+=" It is unlikely that the chosen hand would not be correct as there is a good discrimination between the hands."
        self.rv_report.Text('&emsp;'+result2)
        if not self.parent_process or not self.parent_process.ccp4i2:
          self.rv_report.Text("<BR><i><small>Output:</small></i>")
          hand1str='Chosen hand' if self.chosen_hand else 'Hand 1'
          self.rv_files1 = self.rv_report.DataFile(self.out.Get('model').GetFileName(), "xyz", "<small>{}: Substructure and experimental map</small>".format(hand1str))
          out_mtz1 = crvapi.SetMtzFWT(self,self.out.Get('mapcoef',typ=('best','combined'),filetype='mtz'))
          self.rv_files1.DataFile(out_mtz1, "hkl:map")
          #self.rv_files1.DataFile(self.out.Get('mapcoef',typ='best',filetype='map').GetFileName('map'), "hkl:ccp4_map")
          self.rv_files2 = self.rv_report.DataFile(self.other_phas.out.Get('model').GetFileName(), "xyz", "<small>Other hand: Substructure and experimental map</small>")
          out_mtz2 = crvapi.SetMtzFWT(self,self.other_phas.out.Get('mapcoef',typ=('best','combined'),filetype='mtz'))
          self.rv_files2.DataFile(out_mtz2, "hkl:map", flush=True)
          #self.rv_files2.DataFile(self.other_phas.out.Get('mapcoef',typ='best',filetype='map',conv_opts=['no_rewrite']).GetFileName('map'), "hkl:ccp4_map", flush=True)


  def RunPostprocess(self,restore=True,*args,**kwargs):
    # print log info
    self.Info('Phasing CLD score for hand1 is {0} and for hand 2 is {1}'.format(round(self.cld[0]*100,2),round(self.cld[1]*100,2)))
    result1,result2=None,None
    if self.dmf[0].GetParam('dmcyc')>0 and self.dmf[1].GetParam('dmcyc')>0:
      self.Info('FOM in DM:')
      self.Info('Cycle \tHand1 \tHand2')
      for i,pair in enumerate(zip(self.dmf[0].ph.fom_all,self.dmf[1].ph.fom_all)):
        self.Info('  {0} \t{1} \t{2}'.format(i,pair[0],pair[1]))
      result1='The combined DM FOM and phasing CLD score for hand 1 is {0} and for hand 2 is {1}'.format(round(self.score[0],2),round(self.score[1],2))
      if self.chosen_hand:
        result2='Hand {0} is chosen.'.format(self.chosen_hand)
      else:
        result2='Hand was not chosen, difference between the scores of hands is smaller than the threshold.'
      if self.spacegroup_change:
        result2+='  Spacegroup changed from {0} to {1}.'.format(self.other_phas.out.Get('mapcoef').GetSpacegroup(self),self.spacegroup_change)
      self.Info(result1)
      self.Info(result2)
    # print loggraph info
    self.PrintActualLogGraph(result1,result2)
    process.RunPostprocess(self,restore,*args,**kwargs)
