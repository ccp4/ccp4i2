#!/usr/bin/python
import os,sys
from program import program
import common

class truncate(program):
  name="TRUNCATE"
  binary="truncate"
  modif=''

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    self.objs=[]
    # always convert I->F (if inputted)
    obj = self.inp.Get('fsigf', filetype='mtz', typ='average', col='i')
    if obj:
      self.AddToKey('labin', ('IMEAN='+obj.GetLabel('i'),'SIGIMEAN='+obj.GetLabel('sigi')) )
      self.objs.append(obj)
    obj,N=self.inp.GetTargetObjects('SAD', no_model=True, inten=True, fsftyp='mtz', no_warn=True, accept_anom_nat=True)
    if obj:
      # convert I+/I->F+/F-
      self.AddToKey('labin', ('I(+)='+obj['f+'].GetLabel('i'),'SIGI(+)='+obj['f+'].GetLabel('sigi')) )
      self.AddToKey('labin', ('I(-)='+obj['f-'].GetLabel('i'),'SIGI(-)='+obj['f-'].GetLabel('sigi')) )
      self.SetKey('anomalous','yes')
      self.objs.extend([obj['f+'],obj['f-']])
    if not self.objs:
      common.Error('Nothing (suitable) inputted for I->F conversion by {0}'.format(self.name))
    self.SetArg('hklin', self.objs[0].GetFileName('mtz'))

  def DefineOutput(self):
    for obj in self.objs:
      outo=self.out.AddCopy(obj)
      self.out.SetFileToChild(outo, self.outfilename['mtz'], 'mtz')
      outo.SetLabel( ['f','sigf'] )

  def TreatOutput(self):
    self.SetArg('hklout', self.out.Get('fsigf').GetFileName('mtz'))
    for outo in self.out.GetAll('fsigf'):
      if outo.GetType()=='average':
        self.AddToKey( 'labout', ('F='+outo.GetLabel('f'),'SIGF='+outo.GetLabel('sigf')) )
      elif outo.GetType()=='plus':
        self.AddToKey( 'labout', ('F(+)='+outo.GetLabel('f'),'SIGF(+)='+outo.GetLabel('sigf')) )
      elif outo.GetType()=='minus':
        self.AddToKey( 'labout', ('F(-)='+outo.GetLabel('f'),'SIGF(-)='+outo.GetLabel('sigf')) )
      else:
        common.Error('Unexpected data type for {0}.'.format(self.name))
