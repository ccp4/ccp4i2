#!/usr/bin/python

import sys

from .. import common
from ..program import program


class mtzMADmod(program):
  name="mtzMADmod"
  binary="mtzMADmod"
  ccp4_parsing=True

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    # collect all f+/f- pairs
    self.fpm=[]
    for fp in self.inp.GetAll(filetype='mtz',typ='plus',col='sigf'):
      fm=self.inp.Get(filetype='mtz',typ='minus',col='sigf',dname=fp.dname,xname=fp.xname)
      if fm and fm.dname not in [fpm[1].dname for fpm in self.fpm]:
        self.fpm.append((fp,fm))
    for i,(fp,fm) in enumerate(self.fpm):
      self.AddToKey( 'LABI', ( 'F{0}(+)={1}'.format(i+1,fp.GetLabel('f')), 'SIGF{0}(+)={1}'.format(i+1,fp.GetLabel('sigf')), 
                             'F{0}(-)={1}'.format(i+1,fm.GetLabel('f')), 'SIGF{0}(-)={1}'.format(i+1,fm.GetLabel('sigf')), ) )
    try:
      self.SetArg('hklin', fp.GetFileName('mtz'))
    except (UnboundLocalError,):
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      common.Error('No (suitable) mtz input to {0}'.format(self.name))

  def DefineOutput(self):
    outmtz=self.outfilename['mtz']
    self.fad = []
    for fp,fm in self.fpm:
      aver  = self.out.AddNew( 'fsigf', outmtz, typ='average', xname=fp.xname, dname=fp.dname )
      delta = self.out.AddNew( 'fsigf', outmtz, typ='delta-anom', xname=fp.xname, dname=fp.dname )
      fad_flat = set().union(*self.fad)
      aver.SetLabel(['f','sigf'], bad_lbl_obj=fad_flat)
      delta.SetLabel(['f','sigf'], bad_lbl_obj=fad_flat)
      self.fad.append((aver,delta))

  def TreatOutput(self):
    self.SetArg('hklout', self.outfilename['mtz'])
    for i,(aver,delta) in enumerate(self.fad):
      self.AddToKey( 'LABO', ( 'F{0}={1}'.format(i+1,aver.GetLabel('f')), 'SIGF{0}={1}'.format(i+1,aver.GetLabel('sigf')), 
                             'D{0}={1}'.format(i+1,delta.GetLabel('f')), 'SIGD{0}={1}'.format(i+1,delta.GetLabel('sigf')) ) )
