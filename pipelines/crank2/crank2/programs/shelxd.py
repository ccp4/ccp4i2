#!/usr/bin/python
import os,sys,shutil
from program import program
import common

class shelxd(program):
  name="SHELXD"
  binary="shelxd"
  modif='-'
  arg_divided=False
  ccp4_keys_parsing=True
  never_merge=True
  stat={}
  stat['ccweak'] = common.stats(name='correlation coef. (weak reflections mean)', 
    regexp=r"CC All\/Weak\s+\S+\s+\/\s+(\S+?)\,?\s+")
  stat['cc'] = common.stats(name='correlation coef. (all reflections mean)', 
    regexp=r"CC All\/Weak\s+(\S+)\s+\/\s+\S+?\,?\s+")
  stat['try'] = common.stats(name='try number', regexp=r"Try\s+(\d+),")
  stat['patfom'] = common.stats(name='Patterson FOM', 
    regexp=r"\s+PATFOM\s+(\S+)")
  stat['cutoff'] = common.stats(name='used high resolution cutoff', regexp=r"\s*SHEL\s+dmax\s+\S+\s+dmin\s+(\S+)")
  stat['expired'] = common.stats(name='expired - update from SHELX web needed', regexp=r'version has expired - please update it', accept_none=True)
  stat['error'] = common.stats(name='error', regexp=r' \*\* CANNOT', accept_none=True)
  stat['command_line_switches'] = common.stats(name='command line switches outputted', regexp=r'Command line switches', accept_none=True)
  # this is regexp for the output RES, not log
  stat['occup'] = common.stats(name='atom occupancies', regexp=r'\S+\s+\d\s+-?\d\.\d+\s+-?\d\.\d+\s+-?\d\.\d+\s+(\d\.\d+)\s+\d\.\d+', multiple=True)
  references = ( "Schneider TR, Sheldrick GM (2002) Substructure solution with SHELXD. " +
                 "Acta Cryst. D58, 1772-1779.", )


  def CheckBinary(self,silent=False):
    bin_issue=self.name+' binary problem: '
    ok=program.CheckBinary(self, try_run=True, silent=silent)
    if not silent:
      if self.GetStat('expired'):
        common.Error(bin_issue+self.stat['expired'].name)
      if not self.GetStat('command_line_switches'):
        common.Error(bin_issue+' command line switches not outputted - please check the binary')
    return ok

  def Interact_output(self, line, popen, empty):
    if self.GetStatGrep('error',from_str=line):
      common.Error('{0} failed with error message: {1}'.format(self.name,line))
    self.process.UpdateInfo(line)

  def Stop(self,popen=None):
    with open(self.project+'.fin', 'w') as f:
      f.write('fin.')

  def TreatInput(self):
    self.fa=self.inp.Get('fsigf', typ='fa', filetype='hkl')
    if not self.fa:
      common.Error('FA values not inputted to {0}'.format(self.name))
    ins=self.inp.Get('datafile', filetype='ins')
    if not ins:
      common.Error('.ins file not inputted to {0}'.format(self.name))
    # make sure the inp and hkl file names are as shelxd expects them
    if self.fa.GetFileName('hkl').endswith('.hkl'):
      self.project=self.fa.GetFileName('hkl')[:-4]
    else:
      self.project=self.fa.GetFileName('hkl')
      shutil.copy(self.project, self.project+'.hkl')
      self.fa.GetFile('hkl').name+='.hkl'
    if ins.GetFileName('ins')!=self.project+'.ins':
      shutil.copy(ins.GetFileName('ins'), self.project+'.ins')
      ins.GetFile('ins').name+='.ins'
    # we need to do the relative path manually here as such a file does not exist...
    self.SetArg(os.path.relpath(self.fa.GetFileName('hkl')[:-4], self.rundir), ignore_modif=True)


  def TreatParams(self):
    if self.process.GetVirtPar('num_trials') and not self.GetKey('NTRY'):
      self.SetKey('NTRY', self.process.GetVirtPar('num_trials'))
    if self.process.GetVirtPar('num_dsul') and not self.IsKey('DSUL'):
        self.SetKey('DSUL', self.process.GetVirtPar('num_dsul'))
    if self.GetKey('DSUL') is not None and self.GetKey('DSUL')<=0:
      self.SetKey('DSUL', False, keep_previous=False)
    if not self.GetKey('FIND'):
      num_atoms = self.process.GetParam('num_atoms')
      if not num_atoms and self.inp.Get(typ='substr',has_num_atoms=True):
        num_atoms = self.inp.Get(typ='substr',has_num_atoms=True).exp_num_atoms
      if num_atoms:
        dsul = self.GetKey('DSUL') if self.GetKey('DSUL') else 0
        self.SetKey('FIND', num_atoms-min(dsul,num_atoms/2))
    if self.GetKey('DSUL') and self.IsKey('FIND') and self.GetKey('FIND')<self.GetKey('DSUL'):
      if self.process.IsInputtedParam('num_dsul') or not self.process.IsParam('num_dsul'):
        common.Error('Wrong input to {0}, FIND < DSUL: {1} peaks cannot be resolved into {2} disulphides. Please adjust DSUL/FIND or num_dsul/num_atoms parameters of substrdet.'.format(self.name,self.GetKey('FIND'),self.GetKey('DSUL')))
      else:
        self.SetKey('DSUL', self.GetKey('FIND'), keep_previous=False)
    if self.process.GetParam('high_res_cutoff') and not self.GetKey('SHEL'):
      self.SetKey('SHEL', (999., self.process.GetParam('high_res_cutoff')))
    if self.process.IsParam('min_dist_atoms') and not self.IsKey('MIND'):
      if self.process.IsParam('min_dist_symm_atoms'):
        self.SetKey('MIND', (self.process.GetParam('min_dist_atoms'),self.process.GetParam('min_dist_symm_atoms')))
      else:
        self.SetKey('MIND', self.process.GetParam('min_dist_atoms'))
    if self.process.GetParam('num_threads') and not self.GetArg('t'):
      self.SetArg('t', self.process.GetParam('num_threads'))
    program.TreatParams(self)


  def DefineOutput(self):
    #if self.inp.Get('model'):
    #  atomtypes,atomtype1=self.inp.Get('model').GetAtomTypes(),self.inp.Get('model').GetAtomType()
    #  xname=self.inp.Get('model').GetCrystalName()
    #  self.out.AddNew('model', self.project+'.pdb', filetype='pdb', typ='substr', \
    #                  atomtypes=atomtypes, atomtype1=atomtype1, xname=xname)
    if self.inp.Get('model',typ='substr',xname=self.fa.xname):
      mod=self.out.AddCopy(self.inp.Get('model',typ='substr',xname=self.fa.xname))
      self.out.SetFileToChild(mod, self.project+'.pdb', 'pdb')
    else:
      common.Warning('Substructure atom types and crystal name not known!')
      mod=self.out.AddNew('model', self.project+'.pdb', filetype='pdb', typ='substr', xname=self.fa.xname)
    self.out.AddFileToChild(mod, self.project+'.res', 'res')
    # output hkl from shelxc to be used by shelxe
    if self.inp.Get('fsigf',filetype='hkl',typ='average',try_convert=False):
      self.out.AddCopy(self.inp.Get('fsigf',filetype='hkl',typ='average',try_convert=False))
