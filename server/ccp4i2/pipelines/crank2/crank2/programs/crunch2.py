#!/usr/bin/python
import os,sys
from program import program
import common

class crunch2(program):
  name="CRUNCH2"
  binary="crunch2"
  ccp4_parsing=True
  stat={}
  stat['patfom'] = common.stats(name='FOM', regexp=r"The figure of merit is:\s+(\S+)")
  stat['try'] = common.stats(name='try number', regexp=r"Start trial number\s+(\d+)")
  stat['cyc_in_try'] = common.stats(name='number of cycles used in the specified trial', 
    regexp=r"Trial number\s+{0} is stopped at cycle\s+(\d+)")
  stat['fom_unopt'] = common.stats(name='FOM unoptimized', \
    regexp=r"EVAL: Correlation coefficient:\s+\S+\s+Fom:\s+(\S+)")
  references = ( "de Graaff RAG, Hilge M, van der Plas JL and Abrahams JP (2001) " +
                 "Matrix methods for solving protein substructures of chlorine and " +
                 "sulfur from anomalous data. Acta Cryst. D57, 1857-1862.", )


  def Interact_output(self, line, popen, empty):
    #if not empty:
    self.process.UpdateInfo(line,popen)

  def Stop(self,popen):
    popen.terminate()
    popen.wait()
    popen.returncode=0

  def TreatInput(self):
    fa = self.inp.Get('fsigf',typ='fa',col='f')
    self.ea = self.inp.Get('fsigf',typ='fa',filetype='drear',inp_cont=fa)
    if self.ea is None:
      common.Error('FA/EA values not inputted to {0}'.format(self.name))
    self.SetArg('HKLIN', self.ea.GetFileName('drear'))
    # while it is told not to affect results much, it is still required else NaNs are returned...
    self.SetKey('SCAT', 15)
    # set cell and spacegroup
    if not self.IsKey('CELL') or not self.IsKey('SYMM'):
      self.SetKey('CELL', ' '.join((str(c) for c in fa.GetCell(self))))
      self.SetKey('SYMM', self.ea.GetSpacegroupNumber(self))
    # expected number of atoms from substr. object
    if self.inp.Get(typ='substr',has_num_atoms=True) and not self.GetKey('NATO'):
      self.SetKey('NATO', self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms)
    # patterson input
    if self.inp.Get('model',typ='patterson'):
      self.SetArg('modelin', self.inp.Get('model',typ='patterson').GetFileName())
      self.SetKey('ICOO',1)
    if not self.GetKey('TRY'):
      common.Error('{0} requires number of trials at input.'.format(self.name))
    if not self.GetKey('NATO'):
      common.Error('{0} requires expected number of atoms at input.'.format(self.name))

  def TreatParams(self):
    if self.process.GetVirtPar('num_trials') and not self.GetKey('TRY'):
      self.SetKey('TRY', (1,self.process.GetVirtPar('num_trials')))
    if self.process.GetParam('num_atoms') and not self.GetKey('NATO'):
      self.SetKey('NATO', self.process.GetParam('num_atoms'))
    if self.process.GetParam('high_res_cutoff') and not self.GetKey('SHEL'):
      self.SetKey('STLM', 1.0/(2.0*self.process.GetParam('high_res_cutoff')))
    program.TreatParams(self)

  def DefineOutput(self):
    # hardcoded and not uniquely given output... will just assume the first trial here and update later...
    outfile='trial1a.pdb'
    if self.inp.Get('model',typ='substr',xname=self.ea.xname):
      mod=self.out.AddCopy(self.inp.Get('model',typ='substr',xname=self.ea.xname))
      self.out.SetFileToChild(mod, outfile, 'pdb')
    else:
      common.Warning('Substructure atom types and crystal name not known!')
      self.out.AddNew('model', outfile, filetype='pdb', typ='substr')

  def GetTrialPdb(self,trial):
    # returns input trial's pdb name (incl. path) 
    # must be called after DefineOutput(), otherwise there is an error.
    return self.out.Get('model').GetFileName().replace('trial1a.pdb','trial'+str(trial)+'a.pdb')
