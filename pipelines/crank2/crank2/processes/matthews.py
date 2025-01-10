#!/usr/bin/python

import math
import sys

from .. import common
from ..process import process
from ..program import program


class matthews(process):
  name="estimation of Matthews coefficient, solvent content and number of monomers"
  short_name="Matthews estimation"
  supported_progs=["gcx","matthews_coef"]


  def TreatInOutPar(self, set_all_par=False):
    if not self.GetProg(supported=True):
      if self.inp.Get('fsigf',filetype='mtz'):
        self.AddProg('gcx')
      else:
        self.AddProg('matthews_coef')
    process.TreatInOutPar(self,set_all_par)

  def RunPreprocess(self, rundir=None, **kwargs):
    process.RunPreprocess(self,rundir,kwargs)
    prog = self.GetProg(supported=True)
    if prog.nick=='matthews_coef' and not prog.IsKey('molweight'):
      if self.inp.Get('sequence') and self.inp.Get('sequence').GetFileName():
        seqwt=self.AddProg('seqwt')
        seqwt.Run()
        prog.SetKey('molweight',seqwt.GetStat('mol_weight'))
      else:
        common.Warning('Neither sequence nor molecular weight inputted: Matthews coef. estimation may be imprecise.')

  def RunBody(self,*args,**kwargs):
    prog = self.GetProg(supported=True)
    if prog.nick=='matthews_coef':
      prog.Run()
    else:
      try:
        prog.Run()
      except program.ProgramRunError as e:
        # gcx does not accept some more obscure cases, for example carriage returns only as eol
        common.Warning('{}.   (will try matthews_coef now instead)'.format(e))
        if hasattr(sys,'exc_clear'): sys.exc_clear()
        self.programs.remove(prog)
        prog=self.AddProg('matthews_coef')
        self.RunPreprocess()
        prog.Run()

  def RunPostprocess(self,restore=True,*args,**kwargs):
    # let's copy the determined values to input sequence object or model object
    obj = self.inp.Get('sequence')  if self.inp.Get('sequence')  else self.inp.Get(has_residues_mon=True)
    if obj:
      prog = self.GetProg(supported=True)
      # we keep the inputted solvent content value from another container
      # however, at least for now the estimated number of molecules in ASU is not affected by solvent input
      if self.inp.Get(has_solvent_content=True):
        obj.solvent_content=self.inp.Get(has_solvent_content=True).solvent_content
      if prog.nick=='matthews_coef':
        probabs = prog.GetStat('probab_matth_all')
        ind = probabs.index(max(probabs))
        if not obj.solvent_content:
          obj.solvent_content=format(prog.GetStat('solvent_content_all')[ind]*0.01,'.3f')
        obj.monomers_asym=prog.GetStat('monomers_asym_all')[ind]
        if not obj.residues_mon:
          obj.residues_mon=len(obj.GetSequenceString())
      else:
        if not obj.solvent_content:
          obj.solvent_content=prog.GetStat('solvent_content')
        obj.monomers_asym=prog.GetStat('monomers_asym')
        obj.residues_mon=prog.GetStat('residues_mon')
    self.PrintActualGraph(obj)
    process.RunPostprocess(self,restore,*args,**kwargs)

  def PrintActualGraph(self, obj):
    if not self.rv_report:
      return
    prog=self.GetProg(supported=True)
    y=[0]
    try:
      x = list( map( int, filter(None, prog.GetStat('monomers_asym_all') ) ) )
      y = list( map( float, filter(None, prog.GetStat('probab_matth_all') ) ) )
      y2 = list( map( float, filter(None, prog.GetStat('solvent_content_all') ) ) )
    except Exception as e:
      common.Warning('Could not grep Matthews distribution used for report, due to error:'+str(e))
      if hasattr(sys,'exc_clear'): sys.exc_clear()
    else:
     if y and not math.isnan(y[0]):  # a gcx issue with nan in probs...
      if prog.nick=='matthews_coef':
        y2 = list( map( lambda q: q*0.01, y2 ) )
      if len(x)>10:
        # remove the smallest values if too many
        x_ext,y_ext,y2 = list(zip( * filter( lambda q: q[1]>0.001, zip(x,y,y2) ) ))
        y2 = list(y2)
        added = 0
      else:
        # add first and last with 0 prob. for a nicer display
        x_ext = [x[0]-1]+x+[x[-1]+1]
        y_ext = [0]+y+[0]
        added = 1
      self.rv_plot = self.rv_report.Plot( 'Matthews coef. probability distribution', 'Number of monomers', 'Matthews coef. probability', block="Matthews probabilities", legendloc='ne', xmin=0, ymin=0., intx=True )
      self.rv_plotline_monom = self.rv_plot.PlotLine( ["x",x_ext], ["Probability of monomers",y_ext], style='bars', custom_x_tick=True )
      self.rv_plotline_solvent = self.rv_plot.PlotLine( ["x2",list(map(lambda x:x+1+added, range(len(x))))], ["Corresponding solvent content",y2], flush=True )
    self.rv_report.Text("<i><small>Result:</small></i>")
    self.rv_report.Text("{0} monomers and solvent content of {1} will be assumed.".format(obj.monomers_asym,obj.solvent_content))
    if y and max(y)<0.95 and max(y)>0:
      other_monom = [str(x) for x,y in zip(x,y) if y>0.01 and x!=obj.monomers_asym]
      self.rv_report.Text('However, there is a significant likelihood that this assumption might be incorrect and it is suggested to test rerunning the job with other likely number of monomers (eg {}) and corresponding solvent contents.'.format(','.join(other_monom)))
