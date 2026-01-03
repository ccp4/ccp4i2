#!/usr/bin/python
import os,sys
from process import process, crvapi
import threading
import common
par=common.parameter

class phdmmb(process):
  name="phasing, density modification and model building"
  short_name="density mod. & poly-Ala tracing"
  supported_progs = ['shelxe']
  supported_params={}
  supported_params['catch_output'] = par( desc='Catch program output continously (True by default)', typ=bool )
  supported_params['bigcyc'] = par( desc='Number of model building cycles', typ=int )
  supported_params['dmcyc'] = par( desc='Number of density modification cycles', typ=int )
  supported_params['solvent_content'] = par( desc='(Expected) solvent fraction of crystal', typ=float, share=True )
  supported_params['handdet'] = par( desc='Run with both hands (True by default)', typ=bool )
  supported_params['substr_in_native'] = par( desc='Whether substr. atoms are present in native (if applicable)', typ=bool )
  supported_params['threshold_stop'] = par( desc='Stop when this CC threshold is reached', typ=(float,bool) )
  supported_params['threshold_hand_stop'] = par( desc='Stop the other hand when this CC threshold is reached', typ=(float,bool) )
  supported_params['threshold_failure'] = par( desc='Failure is reported if the final CC is lower than this threshold', typ=(float,bool) )
  supported_params['thorough_build'] = par( desc='Thorough building (False by default)', typ=bool )

  def Init(self):
    self.stats_table={}
    self.stats_table['1'], self.stats_table['2'] = {}, {}
    self.stat_nicks = ('cc','contrast','res_built','res_built_per_frag')
    self.stat_names = ("Correlation coefficient", "Contrast", "Residues", "Residues per fragment")
    self.stop_file='stop_file'

  def TreatInOutPar(self, set_all_par=False):
    shelxe=self.GetOrAddProg('shelxe')
    if self.GetParam('catch_output') is False:
      self.Info('Early stop and continuous reporting are disabled.')
    else:
      self.GetProg(supported=True).interact_output=True
    if self.GetVirtPar('handdet') is None:
      self.SetVirtPar('handdet', True)
    # rewriting the program defaults here!
    if self.IsTrueOrNoneParam('bigcyc'):
      self.SetVirtPar('bigcyc', 8)
    if self.IsTrueOrNoneParam('threshold_stop'):
      self.SetParam('threshold_stop', 40.)
    if self.IsTrueOrNoneParam('threshold_hand_stop'):
      self.SetParam('threshold_hand_stop', 30.)
    if self.IsTrueOrNoneParam('threshold_failure') and (not self.GetCrankParent() or not self.GetCrankParent().GetProcesses('comb_phdmmb')):
      self.SetParam('threshold_failure', 20.)
    if self.GetVirtPar('substr_in_native') is None:
      # we try to "guess" whether the heavy atoms are in the "native"...  -not here anymore
      at_obj = self.inp.Get(has_atomtypes=True)
      #if at_obj and not self.GetProg(supported=True).IsArg('h') and at_obj.GetAtomType() not in ('SE','I','BR'):
      #if at_obj and not self.GetProg(supported=True).IsArg('h') and at_obj.GetAtomType() in ('S','FE'):
      # substr_in_native trivially true if no native inputted;  needed for the current shelx (S,SE excluded)
      if not self.GetProg(supported=True).IsArg('h') and at_obj and at_obj.GetAtomType() not in ('S','SE') and \
         (not self.inp.Get(is_native=True) or self.inp.Get('model',typ=('substr','partial+substr'),is_native=True)):
        self.SetVirtPar('substr_in_native', True)

  def RunHand(self,hand_str,lock=None):
    if lock:
      lock.acquire()
    prog = getattr(self,'prog'+hand_str)
    prog.hand=hand_str
    prog.failed=False
    self.Info(" Running shelxe with hand {0} ".format(hand_str))
    try:
      prog.Run( rundir=os.path.join(self.rundir,'hand'+hand_str), lock=lock )
    except Exception as e:
      common.Warning('Hand {0} tracing failed with error message: {1}'.format(hand_str,e))
      if lock and lock.locked():  lock.release()
      prog.failed=str(e)
    else:
      if prog.GetStat('untraceable'):
        self.Info(" shelxe could not trace the map for hand {0}.".format(hand_str))
      elif self.GetParam('bigcyc')==0:
        self.Info(" shelxe finished density modification with hand {0}.".format(hand_str))
      else:
        self.Info(" shelxe finished tracing with hand {0}, best correlation coef. {1}.".format(hand_str,prog.GetStat('best_cc')))
    #self.stats_table.append( [prog.GetStat("cc"),prog.GetStat("contrast"),prog.GetStat("res_built")] )
    #self.PrintActualLogGraph()

  def RunBody(self,*args,**kwargs):
    if self.stop_file and self.GetCrankParent() and self.GetCrankParent().rundir:
      self.stop_file = os.path.join(self.GetCrankParent().rundir,self.stop_file)
    if self.GetProg('shelxe'):
      self.prog1 = self.GetProg('shelxe')
      if self.GetVirtPar('handdet'):
        self.prog2 = self.AddProgCopy(self.prog1, deeper_copy=True)
        self.prog2.SetArg('i')
        lock=threading.Lock()
        run1 = threading.Thread(target=self.RunHand, args=('1',lock))
        run1.start()
        self.RunHand('2',lock)
        run1.join()
      else:
        self.RunHand('1')

  def UpdateInfo(self,line,prog):
    stop_message=None
    if ((hasattr(self,'ccp4i2job') and self.ccp4i2job.testForInterrupt()) or os.path.isfile(self.stop_file)) and \
       (not hasattr(self.prog1,'stopped') or not hasattr(self.prog2,'stopped')):
      stop_message='Stopping early on user request!'
      self.Info(stop_message)
      self.prog1.Stop(),  self.prog2.Stop()
      self.PrintActualLogGraph(stop_message=stop_message)
    if not hasattr(prog,'prev_line'):  prog.prev_line=''
    for stat_nick in self.stat_nicks:
      sv=prog.GetStatGrep(stat_nick, from_str=line if stat_nick!='res_built_per_frag' else prog.prev_line+'\n'+line)
      if sv and (stat_nick!='res_built_per_frag' or filter(None,sv)):
        if stat_nick=='res_built_per_frag':
          sv = [int(l) for l in sv if l]
          stat_nick, sv = 'frag_length', float(sum(sv))/max(len(sv),1)
        if not stat_nick in self.stats_table[prog.hand]:
          self.stats_table[prog.hand][stat_nick]=[]
        self.stats_table[prog.hand][stat_nick].append(sv)
        if stat_nick=='cc':
          self.Info('CC after {0}. model building with hand {1} is {2}'.format(
                      len(self.stats_table[prog.hand][stat_nick]),prog.hand,sv) )
          if not hasattr(getattr(self,'prog'+prog.hand),'stopped') and len(self.stats_table[prog.hand][stat_nick])>2:
            if self.IsParam('threshold_stop') and sv>self.GetParam('threshold_stop'):
              self.prog1.Stop()
              if self.GetVirtPar('handdet'):  self.prog2.Stop()
              stop_message='CC parameter threshold ({0}) reached. Stopping after the current cycle.'.format(self.GetParam('threshold_stop'))
              self.Info(stop_message)
            oh,ioh = (getattr(self,'prog1'),'1')  if prog.hand=='2'  else (getattr(self,'prog2'),'2')
            if self.GetVirtPar('handdet') and self.IsParam('threshold_hand_stop') and not oh.GetStat('untraceable') and \
               sv>self.GetParam('threshold_hand_stop') and not hasattr(oh,'stopped') and sv>self.stats_table[ioh][stat_nick][-1]:
              oh.Stop()
              stop_message='CC parameter threshold ({0}) for hand determination reached. Stopping the other hand after the current cycle.'.format(self.GetParam('threshold_hand_stop'))
              self.Info(stop_message)
        self.PrintActualLogGraph(prog.hand, stat_nick, stop_message=stop_message)
    prog.prev_line=line

  def PrintActualLogGraph(self, hand='1', stat_nick=None, finished=False, stop_message=None):
    prog=self.GetProg(supported=True)
    if finished:
      #result1="Maximum contrast for hand 1 is {0} and for the other hand is {1}".format(max(self.stats_table[0][0]),max(self.stats_table[1][0]))
      #result2="Maximum number of residues built for hand 1 is {0} and for the other hand {1}".format(max(self.stats_table[0][2]),max(self.stats_table[1][2]))
      if not self.prog1.GetStat('untraceable') and not self.prog1.failed: 
        if self.prog1.GetStat('best_cc',accept_none=True):
          result1="Maximum correlation coefficient for hand 1 is {0};".format(self.prog1.GetStat('best_cc'))
        else:
          result1="Map contrast for hand 1 is {0};".format(self.prog1.GetStat('contrast'))
      else:
        result1="Unable to trace map for hand 1;"
      if self.GetVirtPar('handdet'):
        if not self.prog2.GetStat('untraceable') and not self.prog2.failed:
          if self.prog2.GetStat('best_cc',accept_none=True):
            result1+=" max. correlation coef. for the other hand is {0}. Hand {1} is chosen.".format(self.prog2.GetStat('best_cc'),self.chosen_hand)
          else:
            result1+=" map contrast for the other hand is {0}. Hand {1} is chosen.".format(self.prog2.GetStat('contrast'),self.chosen_hand)
        else:
          result1+=" unable to trace map for hand 2."
          if self.chosen_hand:  result1+=" Hand 1 is chosen."
      result2,result3='',''
      if not self.chosen_hand:
        result2="Tracing unsuccessful."
      elif 'res_built' in self.stats_table[self.chosen_hand]:
        result2="Max. number of residues built (for the chosen hand) is {0}".format(max(self.stats_table[self.chosen_hand]['res_built']))
      if self.chosen_hand and self.best_cc and self.GetParam('threshold_failure') and self.best_cc<self.GetParam('threshold_failure'):
        result3='<BR>&emsp;Correlation coefficient too low, the traced model is most likely wrong.'
    if self.opened_loggraph:
      self.GetLogGraphHandle().seek(0,0)
      self.LGInfo(self.GetCCP4Header())
      self.LGInfo(self.GetLogGraphPrograms())
      if self.stats_table:
        for si,stat in enumerate(self.stat_nicks):
          self.LGInfo('\n $TABLE : {0} of each hand per cycle:'.format(self.stat_names[si]))
          if stat in self.stats_table['1'] and stat in self.stats_table['2']:
            self.LGInfo('$GRAPHS    : {0} of each hand per cycle :A:1,2,3:\n$$'.format(self.stat_names[si]))
            self.LGInfo('Cycle   Hand1   Hand2 $$ $$')
            for i in range(min(len(self.stats_table['1'][stat]),len(self.stats_table['2'][stat]))):
              self.LGInfo('{0:3}  {1:8.3f} {2:8.3f}'.format(i+1, self.stats_table['1'][stat][i], 
                                                                 self.stats_table['2'][stat][i]))
          elif self.stats_table['1'] and stat in self.stats_table['1']:
            self.LGInfo('$GRAPHS    : {0} of each hand per cycle :A:1,2:\n$$'.format(self.stat_names[si]))
            self.LGInfo('Cycle   Hand1 $$ $$')
            for i,sh1 in enumerate(self.stats_table['1'][stat]):
              self.LGInfo('{0:3}  {1:8.3f}'.format(i+1, sh1))
          self.LGInfo('$$\n')
        if finished:
          self.LGInfo('\n$TEXT:Result: $$ Correlation coefficient results $$')
          self.LGInfo(result1+'\n$$')
          self.LGInfo('\n$TEXT:Result: $$ Model building results $$')
          self.LGInfo(result2+'\n$$')
    if self.rv_report is not None:
      cc_thr, frag_length_thr = 20., 10.
      if stat_nick in self.stats_table[hand] and stat_nick in ('cc','res_built','frag_length','contrast'):
        if not hasattr(self,'rv_plot_cc'):
          self.rv_plot_cc = self.rv_report.Plot( 'Correlation coef. vs building cycle', "Cycle",  \
                                    "CC of trace against data [%]", block="Statistics vs cycle", legendloc='nw', xmin=0., ymin=0., ymax=cc_thr+1, intx=True )
          self.rv_plot_res_built = self.rv_plot_cc.parent.Plot( 'Residues vs building cycle', \
                                    "Cycle", "Backbone residues built", legendloc='nw', xmin=0., ymin=0., intx=True )
          self.rv_plot_frag_length = self.rv_plot_cc.parent.Plot( 'Chain length vs building cycle', \
                                    "Cycle", "Average built fragment length", legendloc='nw', xmin=0., ymin=0., ymax=frag_length_thr+1, intx=True )
          self.rv_plot_contrast = self.rv_plot_cc.parent.Plot( 'Contrast in the initial density modification', \
                                    "Density modification cycle", "Map contrast", legendloc='nw', xmin=0., ymin=0., intx=True )
          self.rv_cc, self.rv_res_built, self.rv_frag_length, self.rv_contrast = {}, {}, {}, {}
        if stat_nick and hasattr(self, 'rv_'+stat_nick) and (stat_nick!='contrast' or 'cc' not in self.stats_table[hand]):
          rv_line_stat, rv_plot_stat = getattr(self, 'rv_'+stat_nick), getattr(self, 'rv_plot_'+stat_nick)
          if not hand in rv_line_stat:
            rv_line_stat[hand] = rv_plot_stat.PlotLine( ["x"], ["hand"+hand], color="#4BB2C5" if hand=='1' else "#EAA228" )
          stat = self.stats_table[hand][stat_nick]
          if len(stat)==1:
            rv_line_stat[hand].CustomTick( "x", 0, "0" )
          if len(stat)==10:
            rv_line_stat[hand].Reset(reset_data=False)
          rv_line_stat[hand].Data( len(stat), stat[-1], x_type=int, custom_x_tick=len(stat) if len(stat)<10 else False )
          if len(stat)==1:
            rv_line_stat[hand].CustomTick( "x", 2, "2" )
          if '1' in rv_line_stat and ('2' in rv_line_stat or not self.GetVirtPar('handdet') or self.prog2.GetStat('untraceable')):
            if stat_nick+'_thr' in vars():
              stats2=self.stats_table['2'][stat_nick] if '2' in rv_line_stat else []
              rv_plot_stat.SetProperties(ymax=max([vars()[stat_nick+'_thr'],]+self.stats_table['1'][stat_nick]+stats2)+1)
            if stat_nick=='cc':
              if not hasattr(self,'rv_cc_shadow'):
                self.rv_cc_shadow = self.rv_plot_cc.PlotLine( ["x",(-0.01,)], ["Low CC area",(cc_thr,)], fill=1 )
              self.rv_cc_shadow.Data( len(stat)+1, cc_thr )
            if stat_nick=='frag_length':
              if not hasattr(self,'rv_frag_shadow'):
                self.rv_frag_shadow = self.rv_plot_frag_length.PlotLine( ["x",(-0.01,)], ["Small fragment length",(frag_length_thr,)], fill=1 )
              self.rv_frag_shadow.Data( len(stat)+1, frag_length_thr )
          crvapi.Flush()
      if stop_message:
        self.rv_report.Text("<i><small>Stop:</small></i>")
        self.rv_report.Text('&emsp;'+stop_message)
      if finished:
        self.rv_report.Text("<i><small>Result:</small></i>")
        self.rv_report.Text('&emsp;'+result1), self.rv_report.Text('&emsp;'+result2, flush=True)
        if result3:  self.rv_report.Text('&emsp;'+result3, flush=True)
        if not self.parent_process or not self.parent_process.ccp4i2:
          if self.chosen_hand and self.GetParam('bigcyc')>0:
            self.rv_report.Text("<BR><i><small>Output:</small></i>")
            self.rv_files1 = self.rv_report.DataFile(self.out.Get('model').GetFileName(), "xyz", "<small>Built model (chosen hand)</small>")
            if self.GetVirtPar('handdet') and getattr(self,'prog'+str(1+abs(int(self.chosen_hand)-2))).out.Get('model',typ=('partial','partial+substr')):
              self.rv_files2 = self.rv_report.DataFile(getattr(self,'prog'+str(1+abs(int(self.chosen_hand)-2))).out.Get('model',typ=('partial','partial+substr')).GetFileName(), "xyz", "<small>Built model (other hand)</small>")
          # shelxe PHS on request from Andrea, esp. for ccp4-online
          self.rv_files3 = self.rv_report.DataFile(self.out.Get(filetype='phs').GetFileName('phs'), "text", "<small>Phases in PHS (chosen hand)</small>")
          # shelxe input parameters on request from George, esp. for ccp4-online
          with open(prog.GetScrFileName()) as f:
            par=f.read().partition('shelxe')
            if par[1]:
              self.rv_report.Text("<BR><i><small>Input:</small></i> <pre><xx-small> {0} {1} </xx-small></pre>".format(par[1],par[2]))
          #self.rv_files3 = self.rv_report.DataFile(prog.GetScrFileName(), "text", "<small>SHELXE input parameters</small>")
        #self.rv_files.DataFile(self.out.Get('mapcoef',typ='combined').GetFileName(), "hkl:map", flush=True)


  def RunPostprocess(self,restore=True,*args,**kwargs):
    if self.prog1.failed and (not self.GetVirtPar('handdet') or self.prog2.failed):
      common.Error(self.prog1.failed)
    self.chosen_hand = None  if self.prog1.GetStat('untraceable') or self.prog1.failed  else '1'
    enant_mtz=None
    if self.GetVirtPar('handdet'):
      self.score = [ self.prog1.GetStat('best_cc',accept_none=True), self.prog2.GetStat('best_cc',accept_none=True) ]
      if not any(self.score):
        self.score = [ self.prog1.GetStat('contrast',accept_none=True), self.prog2.GetStat('contrast',accept_none=True) ]
      self.score = [ s  if s is not None else 0.0  for s in self.score ]
      if not self.prog2.GetStat('untraceable') and not self.prog2.failed and self.score[0]<self.score[1]:
        self.chosen_hand = '2'
      # change enantio if needed;  we rely on shelxe PDB here
      if self.prog2.out.Get('model',typ=('partial','partial+substr')) and \
         os.path.isfile(self.prog2.out.Get('model',typ=('partial','partial+substr')).GetFileName('pdb')) and \
         self.inp.Get('fsigf',filetype='mtz'):
        spgr1 = self.inp.Get('fsigf',filetype='mtz').GetSpacegroup(self)
        pdbcur=self.GetOrAddProg('pdbcur')
        with open(self.prog2.out.Get('model',typ=('partial','partial+substr')).GetFileName('pdb')) as f:
          spgr2 = pdbcur.GetStatGrep('spgr',from_str=f.read())
        if "".join(spgr1.split()) != "".join(spgr2.split()):
          spgr_num1 = self.inp.Get('fsigf',filetype='mtz').GetSpacegroupNumber(self)
          from processes import handdet
          if [list(ep) for ep in handdet.handdet.enant_pairs if spgr_num1 in ep]:
            self.prog2.out.Get('mapcoef',filetype='phs',typ='best').spgr = spgr2
            enant_mtz = handdet.handdet.ChangeHand2Spacegroup(self,spgr1,spgr2,self.inp.GetAll('fsigf'))
    if self.chosen_hand:
      prog = getattr(self,'prog'+self.chosen_hand)
      self.best_cc = prog.GetStat('best_cc',accept_none=True)
      if prog.out.Get('model',typ='substr'):
        self.out.Set(prog.out.Get('model',typ='substr'))
      self.out.Set(prog.out.Get('mapcoef',filetype='phs',typ='best'))
      self.spacegroup_change=False
      if enant_mtz and self.chosen_hand=='2':
        self.out.Set(enant_mtz)
        self.spacegroup_change=enant_mtz[0].GetSpacegroup(self)
      if prog.out.Get('model',typ=('partial','partial+substr')):
        self.out.Add(prog.out.Get('model',typ=('partial','partial+substr')))
      self.other_phas = getattr(self,'prog1' if self.chosen_hand=='2' else 'prog2') if self.GetVirtPar('handdet') else None
      if self.other_phas and self.other_phas.GetStat('untraceable'):
        self.other_phas.out.ClearAll(propagate=False)
    self.PrintActualLogGraph( finished=True )
    if not self.chosen_hand:
      common.Warning('Unable to trace the model. Stopping.')
      self.out.ClearAll(propagate=False)
      self.other_phas = self.prog2 if self.GetVirtPar('handdet') else None
      if self.other_phas: self.other_phas.out.ClearAll(propagate=False)
      if self.prog1.failed and (not self.GetVirtPar('handdet') or self.prog2.failed):
        common.Error('Model tracing failed with error message: {0}'.format(self.prog1.failed))
      # not returning an error at the moment due to gui's - might be changed later
      else:
        common.Error('Unable to trace the model, map too noisy.', nosuccess=True)
    if self.best_cc and self.GetParam('threshold_failure') and self.best_cc<self.GetParam('threshold_failure'):
      common.Error('CC too low, the traced model is most likely wrong.', nosuccess=True)
    if self.GetProg('shelxe') and self.chosen_hand:
      prog.FixOutput()
      if prog.out.Get('model',typ=('partial+substr')):
        sepsub=self.AddProcess('sepsubstrprot',propagate_inp=False)#,propagate_out=False)
        sepsub.inp.Set(self.out.Get('model',typ='partial+substr',filetype='pdb'))
        sepsub.Run()
    process.RunPostprocess(self,restore,*args,**kwargs)
