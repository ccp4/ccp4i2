#!/usr/bin/python
import os,sys
from program import program
import common

class coot(program):
  name="COOT"
  binary="coot"

  def Init(self):
    self.outfilename = { 'pdb': self.name+'.pdb', \
                         'mtz': self.name+'.mtz' \
                       }
    self.inpfilename = { 'py': self.name+'.py' }
    if os.name == 'nt' and os.environ.get('CCP4_MASTER'):  # may change in the future?  coot.bat causes pop-ups
      self.binary=os.path.join(os.environ.get('CCP4_MASTER'),'WinCoot','bin','coot-bin')

  def TreatInput(self):
    self.ClearAllArgs()
    self.SetArg('--no-graphics')
    self.SetArg('--script')
    self.SetArg(self.inpfilename['py'])

  def TreatOutput(self):
    self.outmodel = self.out.AddCopy(self.inp.Get('model'))
    self.out.SetFileToChild(self.outmodel, self.outfilename['pdb'], 'pdb')
