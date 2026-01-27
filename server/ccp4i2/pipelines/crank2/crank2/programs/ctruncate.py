#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class ctruncate(program):
  name="CTRUNCATE"
  binary="ctruncate"
  modif='-'

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    obj,N=self.inp.GetTargetObjects('SAD', no_model=True, inten=True, fsftyp='mtz', no_warn=True, accept_anom_nat=True)
    if obj is None:
      # convert only I->F
      obj = self.inp.Get('fsigf', filetype='mtz', typ='average', col='i')
      if not obj:
        common.Error('Nothing (suitable) inputted for I->F conversion by {0}'.format(self.name))
      self.SetArg('colin', obj.GetFullLabel('i','sigi'))
      self.objs=(obj,)
    else:
      # convert I+/I->F+/F-
      self.SetArg('colano', obj['f+'].GetFullLabel( 'i', 'sigi', other_objs=(obj['f-'],) ))
      self.objs=(obj['f+'],obj['f-'])
    self.SetArg('hklin', self.objs[0].GetFileName('mtz'))
    if self.inp.Get('sequence') and self.inp.Get('sequence').GetFileName():
      self.SetArg('seqin', self.inp.Get('sequence').GetFileName())

  def DefineOutput(self):
    for obj in self.objs:
      outo=self.out.Add(obj)
      self.out.SetFileToChild(outo, self.outfilename['mtz'], 'mtz')
      # ctruncate hardcodes the label names...
      if outo.GetType()=='average':
        outo.SetLabel( 'f', 'F' )
        outo.SetLabel( 'sigf', 'SIGF' )
      elif outo.GetType()=='plus':
        outo.SetLabel( 'f', 'F(+)' )
        outo.SetLabel( 'sigf', 'SIGF(+)' )
      elif outo.GetType()=='minus':
        outo.SetLabel( 'f', 'F(-)' )
        outo.SetLabel( 'sigf', 'SIGF(-)' )
        # add average if not inputted
        if not [o for o in self.objs if o.GetType()=='average' and o.xname==obj.xname and o.dname==obj.dname]:
          aver = self.out.AddNew( 'fsigf', self.outfilename['mtz'], typ='average', xname=obj.xname, dname=obj.dname )
          aver.SetLabel( 'f', 'FMEAN' )
          aver.SetLabel( 'sigf', 'SIGFMEAN' )
          if 'native' in obj.custom:
            aver.custom.append('native')
      else:
        common.Error('Unexpected data type for {0}.'.format(self.name))
      self.outobj=outo

  def TreatOutput(self):
    self.SetArg('hklout', self.outobj.GetFileName('mtz'))
