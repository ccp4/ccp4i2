#!/usr/bin/python

from .. import common
from ..program import program


class peakmax(program):
  name="PEAKMAX"
  binary="peakmax"
  never_merge=True
  ccp4_parsing=True
  stat={}
  stat['numpeaks'] = common.stats(name='Number of peaks above the threshold',
            regexp=r"the top\s+(\d+)\s+are selected for output")
  stat['heights'] = common.stats(name='peak heights', multiple=True,
            regexp=r"\s+\d+\s+\d+\s+\d+\s+(\d+\.\d+)\s+\d+\s+\d+\s+\d+\s+\d+\.\d+")


  def TreatInput(self):
    self.inpmap=self.inp.Get('mapcoef',filetype='map')
    if not self.inpmap:
      common.Error('No map inputted to program {0}.'.format(self.name))
    self.SetArg('MAPIN',self.inpmap.GetFileName('map'))

  def TreatParams(self):
    if self.process.IsNonFalseVirtPar('rms_threshold') and (not self.GetKey('THRE') or not 'RMS' in self.GetKey('THRE')):
      self.SetKey( 'THRE', ['RMS', self.process.GetVirtPar('rms_threshold')] )
    if self.process.IsNonFalseVirtPar('max_new_atoms') and not self.GetKey('NUMP'):
      self.SetKey( 'NUMPEAKS', self.process.GetVirtPar('max_new_atoms') )
    if self.process.IsNonFalseVirtPar('bfactor') and not self.GetKey('BFAC'):
      self.SetKey( 'BFAC', self.process.GetVirtPar('bfactor') )
    if self.process.IsNonFalseVirtPar('occupancy'):
      if self.GetKey('BFAC') and len(self.GetKey('BFAC',as_list=True))==1:
        self.AddToKey( 'BFAC', self.process.GetVirtPar('occupancy') )
      else:
        common.Warning("Occupancy setting from parent process ignored by {0}.".format(self.name))
    program.TreatParams(self)

  def DefineOutput(self):
    outpdb=self.name+'.pdb'
    self.out.AddNew('model',typ='unknown',xname=self.inpmap.xname,filename=outpdb,filetype='pdb')

  def TreatOutput(self):
    self.SetArg('XYZOUT', self.out.Get('model').GetFileName('pdb'))
