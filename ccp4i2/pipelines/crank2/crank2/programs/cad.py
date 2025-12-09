#!/usr/bin/python
import os,sys
from program import program
import common

class cad(program):
  name="CAD"
  binary="cad"
  never_merge=True
  ccp4_parsing=True

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    # create the list of input filenames
    self.filelist = self.GetFileList(self.inp.GetAll())
    # from an unknown reason, cad sometimes crashes if the renaming is not done with first dataset(s), thus reordering
    ###if self.process and hasattr(self.process,'xdrename'):
    ###  for f in [f for f in self.filelist[::-1] if f in self.process.xdrename]:
    ###    if [1 for xdo,xdn in self.process.xdrename[f].items() if xdo[0]!=xdn[0] or xdo[1]!=xdn[1]]:
    ###      self.filelist.insert(0, self.filelist.pop(self.filelist.index(f)))
    self.lbl,self.lbl_uniq=self.GetUniqueLabels(self.inp.GetAll(),self.filelist)
    for fi,f in enumerate(self.filelist):
      self.SetArg('HKLIN{0}'.format(fi+1), f)
      if self.process and hasattr(self.process,'xdrename') and f in self.process.xdrename:
        for xdo,xdn in self.process.xdrename[f].items():
          if xdo[0]!=xdn[0] or xdo[1]!=xdn[1]:
            self.SetKey('DREN FILE_NUMBER {0}'.format(fi+1), '\'{0} {1}\' \'{2} {3}\''.format(xdo[0],xdo[1],xdn[0],xdn[1]))
      for i,l in enumerate(self.lbl[f]):
        obj = self.inp.Get(filename=f, label=l)
        if obj.xname:
          self.AddToKey('XNAM FILE_NUMBER {0}'.format(fi+1), 'E{0}={1}'.format(i+1,obj.xname))
        if obj.dname:
          self.AddToKey('DNAM FILE_NUMBER {0}'.format(fi+1), 'E{0}={1}'.format(i+1,obj.dname))
        self.AddToKey('LABI FILE_NUMBER {0}'.format(fi+1), 'E{0}={1}'.format(i+1,l))

  def DefineOutput(self):
    out_mtz=self.out.AddNew('datafile', self.outfilename['mtz'], 'mtz')
    for fi,f in enumerate(self.filelist):
      for i,l in enumerate(self.lbl_uniq[f]):
        obj = self.inp.Get(filename=f, label=l)
        self.AddToKey('LABO  FILE_NUMBER {0}'.format(fi+1), 'E{0}={1}'.format(i+1,l))

  def TreatOutput(self):
    self.SetArg('HKLOUT', self.out.datafile[-1].GetFileName('mtz'))
