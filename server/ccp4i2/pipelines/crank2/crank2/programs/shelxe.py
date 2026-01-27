#!/usr/bin/python
import os,sys,shutil
from ..program import program
from .. import common

class shelxe(program):
  name="SHELXE"
  binary="shelxe"
  modif='-'
  arg_divided=False
  never_merge=True
  stat={}
  stat['cycle'] = common.stats(name='cycle', regexp=r"dens.mod.\s+cycle\s+(\d+)")
  stat['contrast'] = common.stats(name='contrast',
                                  regexp=r"<wt>\s+=\s+\S+,\s+Contrast\s+=\s+(\S+?)\,")
  stat['connect'] = common.stats(name='connectivity',
                                 regexp=r"<wt>\s+=\s+\S+,\s+Contrast\s+=\s+\S+?\,\s+Connect.\s+=\s+(\S+)")
  stat['res_built'] = common.stats(name='number of residues built', regexp=r"\s+(\d+)\s+residues left after pruning")
  stat['res_built_per_frag'] = common.stats(name='number of residues built per each fragment', regexp=r"residues left after pruning, divided into chains as follows:\s*\n"+r'(?:\s+\S\:\s+(\d+))?'*40)
  stat['cc'] = common.stats(name='correlation coefficient after building',
                            regexp=r'CC for partial structure against native data =\s+(\S+)\s*\%')
  stat['best_cc'] = common.stats(name='best correlation coefficient after building', 
                                 regexp=r'Best trace \(cycle\s+\S+\s+with CC\s+(\S+)\%\)')
  stat['untraceable'] = common.stats(name='unable to trace the map', regexp=[r'(Unable to trace map)',r'(CC is less than zero - giving up)'], accept_none=True)
  stat['expired'] = common.stats(name='expired - update from SHELX web needed', regexp=r'version has expired - please update it', accept_none=True)
  stat['options_list'] = common.stats(name='options outputted', regexp=r'list of SHELXE options', accept_none=True)
  references = ( "Sheldrick GM (2010) Experimental phasing with SHELXC/D/E:" +
                 " combining chain tracing with density modification.  Acta Cryst. D66, 479-485.", )

  def CheckBinary(self,silent=False):
    bin_issue=self.name+' binary problem: '
    ok=program.CheckBinary(self, try_run=True, silent=silent)
    if not silent:
      if self.GetStat('expired'):
        common.Error(bin_issue+self.stat['expired'].name)
      if not self.GetStat('options_list'):
        common.Error(bin_issue+' options list not outputted - please check the binary')
    return ok

  def Stop(self,popen=None):
    stop_file = self.project+'_i.fin'  if self.GetArg('i')  else self.project+'.fin'
    with open(stop_file, 'w') as f:
      f.write('fin.')
    self.stopped=True

  def Interact_output(self, line, popen, empty):
    self.process.UpdateInfo(line,self)

  def GetNewInpFilePath(self, inpfile, suffix):
    max_length=75
    inpfile=os.path.basename(inpfile)
    if len(inpfile)>max_length:
      inpfile = inpfile[:max_length-4] + inpfile[-4:]
    new_inpfile = os.path.join( self.rundir, inpfile )
    if not new_inpfile.endswith(suffix):
      new_inpfile += suffix
    return new_inpfile

  def CopyInpFile(self, inpfile, new_inpfile, filetype, inpobj=None):
    if os.path.realpath(inpfile)!=os.path.realpath(new_inpfile):
      shutil.copy(inpfile, new_inpfile)
      if inpobj:
        self.inp.AddFileToChild(inpobj, new_inpfile, filetype)

  def TreatInput(self):
    # either substructure or MR model or phases must be inputted
    self.phases, self.res, self.pdb, project_fa_prev = None, None, None, ''
    phases_needed = self.process and self.process.nick in ('mb','dm')
    if not phases_needed:
      self.res=self.inp.Get('model', typ='substr', filetype='res', inp_cont=self.inp.Get('datafile', filetype='ins'))
      if self.res:
        # we copy over the files to the current rundir.
        new_resfile = self.GetNewInpFilePath(self.res.GetFileName('res'), '.res')
        self.project_fa, project_fa_prev = new_resfile[:-4], self.res.GetFileName('res')[:-4]
        self.CopyInpFile(self.res.GetFileName('res'), new_resfile, 'res', self.res)
      # SAD after MR case - this will need to be adjusted to copy to the rundir etc!
      else:
        self.pdb=self.inp.Get('model', typ=('partial','partial+substr'), filetype='pdb')
        if self.pdb:
          new_pdbfile = self.GetNewInpFilePath(self.pdb.GetFileName('pdb'), '.pda')
          self.project, project_prev = new_pdbfile[:-4], self.pdb.GetFileName('pdb')[:-4]
          self.CopyInpFile(self.pdb.GetFileName('pdb'), new_pdbfile, 'pdb', self.pdb)
      if not self.res and not self.pdb:
        common.Error('Neither substructure model (in .res) nor partial model (in .pda) inputted to {0}'.format(self.name))
    # fnat can be either inputted as an object or at least the corresponding filename must exist
    self.fnat=self.inp.Get('fsigf', typ='average', filetype='hkl', is_native=True, try_convert=False)
    if not self.fnat and (not project_fa_prev or not os.path.isfile(project_fa_prev[:-3]+'.hkl')):
      self.fnat=self.inp.Get('fsigf', typ='average', filetype='hkl')
      if not self.fnat:
        common.Error('F/SIGF (native) not inputted to {0}'.format(self.name))
    fnatfile = self.fnat.GetFileName('hkl')  if self.fnat else  project_fa_prev[:-3]+'.hkl'
    new_fnatfile = self.GetNewInpFilePath(fnatfile, '.hkl')
    self.project, project_prev = new_fnatfile[:-4], fnatfile[:-4]
    self.CopyInpFile(fnatfile, new_fnatfile, 'hkl', self.fnat)
    if phases_needed:
      self.phases=self.inp.Get('mapcoef', typ=("combined","best"), filetype='phs', inp_cont=self.inp.Get('fsigf',typ='average',xname=self.fnat.GetCrystalName(),col='f'))
      if not self.phases:
        self.phases=self.inp.Get('mapcoef', typ="weighted", filetype='phs', inp_cont=self.inp.Get('fsigf',typ='average',xname=self.fnat.GetCrystalName(),col='f'))
      if not self.phases:
        common.Error('Phases not inputted or could not be used by {0}'.format(self.name))
      phasfile = self.phases.GetFileName('phs')
      new_phasfile = self.GetNewInpFilePath(self.project, '.phi')
      self.CopyInpFile(phasfile, new_phasfile, 'phs', self.phases)
    # fa can be either inputted as an object or at least the corresponding filename must exist
    else:
      fa=self.inp.Get('fsigf', typ='fa', filetype='hkl', try_convert=False)
      if not fa and not os.path.isfile(project_fa_prev+'.hkl'):
        fa=self.inp.Get('fsigf', typ='fa', filetype='hkl')
        if not fa:
          common.Error('FA values not inputted to {0}'.format(self.name))
      fafile = fa.GetFileName('hkl')  if fa else  project_fa_prev+'.hkl'
      new_fafile = self.GetNewInpFilePath(fafile, '.hkl')
      self.CopyInpFile(fafile, new_fafile, 'hkl', fa)
    # ins can be either inputted as an object or at least the corresponding filename must exist
    ins=self.inp.Get('datafile', filetype='ins')
    if not ins and not os.path.isfile(project_fa_prev+'.ins'):
      common.Error('.ins file not inputted to {0}'.format(self.name))
    insfile = ins.GetFileName('ins')  if ins else  project_fa_prev+'.ins'
    new_insfile = self.GetNewInpFilePath(insfile if not self.phases else self.project, '.ins')
    self.CopyInpFile(insfile, new_insfile, 'ins', ins)
    # we need to do the relative path manually here as such a file does not exist...
    self.project_path = os.path.relpath(self.project, self.rundir)
    if not phases_needed:
      self.project_path_fa = os.path.relpath(self.project_fa, self.rundir)
      self.SetArg( self.project_path, ignore_modif=True )
      self.SetArg( self.project_path_fa, ignore_modif=True )
    else:
      self.SetArg( self.project_path+'.phi', ignore_modif=True )
    # sequence
    if self.inp.Get('sequence') and self.inp.Get('sequence').GetSequenceString():
      with open(os.path.join(self.rundir,self.project+".seq"),'w') as f:
        f.write('>'+self.project+'\n'+self.inp.Get('sequence').GetSequenceString()+'\n')
    # solvent content - if supplied with an input container
    # it seems that sometimes works better with no solvent specified even if correct (fibro)?
    solv_obj = self.inp.Get(has_solvent_content=True)
    if solv_obj and not self.IsArg('s'):
      self.SetArg('s', solv_obj.solvent_content)


  def TreatParams(self):
    if self.process.IsParam('solvent_content') and not self.IsArg('s'):
      self.SetArg('s', self.process.GetParam('solvent_content'))
    if self.process.IsParam('dmcyc') and not self.IsArg('m'):
      self.SetArg('m', self.process.GetParam('dmcyc'))
    if self.process.IsParam('thorough_build') and self.process.GetParam('thorough_build'):
      if not self.process.IsInputtedParam('bigcyc') and not self.IsArg('a'):
        self.SetArg('a',15)
      if not self.IsArg('t'):
        self.SetArg('t',3)
      if not self.IsArg('q'):
        self.SetArg('q')
    if self.process.IsParam('bigcyc') and not self.IsArg('a'):
      self.SetArg('a', self.process.GetParam('bigcyc'))
    if self.process.IsParam('substr_in_native') and self.process.GetParam('substr_in_native') and not self.IsArg('h'):
      self.SetArg('h')
    if self.process.nick=='phdmmb' and not self.IsArg('a'):
      self.SetArg('a')
    if self.process.nick=='mb' and not self.IsArg('a'):
      self.SetArg('a',1)
    if self.process.nick=='mb' and not self.IsArg('m'):
      self.SetArg('m',0)
    if self.process.nick=='mb' and not self.IsArg('q'):
      self.SetArg('q',True)
    if self.process.nick=='dm' and not self.IsArg('m'):
      self.SetArg('m')
    #!!!!!!!
    #self.SetArg('t',3)
    #self.SetArg('a',3)
    if self.process.nick in ('phdmmb','mb') and not self.IsArg('O'):
      self.SetArg('O',True)
    program.TreatParams(self)


  def DefineOutput(self):
    out_path = self.project_path
    if self.GetArg('i'):
      out_path += '_i'
    #x_name = next(o for o in (self.res,self.pdb,self.phases) if o).xname
    if self.process.nick=='dm':
      self.out.AddNew('mapcoef', out_path+'.phs', filetype='phs', typ='densmod')
    else:
      self.out.AddNew('mapcoef', out_path+'.phs', filetype='phs', typ='best')
    if self.res:
      self.out.AddNew('model', out_path+'.hat', filetype='hat', typ='substr', xname=self.res.xname)
    # we take the input substr. pdb as output instead of the refined one!
    # George suggests the shelxd substr. is generally preferable over the refined one
    # we'd need a hat->pdb convertor to use the refined one anyway
    substrin = self.inp.Get('model',typ='substr',filetype='pdb')
    if substrin and self.res:
      self.out.AddCopy(substrin)
    if self.GetArg('a'):
      x_name = self.fnat.xname if self.fnat else self.inp.Get('fsigf',typ='average').xname
      if self.IsArg('h') and self.res and not self.IsArg('O'):  #currently, -O option means that -h ignored!  discussed with Isabel and looks like this will remain - adjust once changed
        partial=self.out.AddNew('model', out_path+'.pdb', typ='partial+substr', filetype='pdb', xname=x_name, atomtypes=substrin.GetAtomTypes())
      else:
        partial=self.out.AddNew('model', out_path+'.pdb', typ='partial', filetype='pdb', xname=x_name)
      if self.fnat and 'native' in self.fnat.custom:
        partial.custom.append('native')


  def FixOutput(self):
    # variuos workarounds for issues with shelxe files...
    # this routine is supposed to be called from the parent process (phdmmb, mb, dm)
    proc=self.process
    if not proc:
      return
    # shelxe sometimes leaves lines with only zeros and asterixes in the phs file - removing these here
    phs=self.out.Get('mapcoef', filetype='phs')
    if phs:
      lines, removed = [], False
      with open(phs.GetFileName('phs')) as f:
        for line in f:
          if '***' in line:
            removed=True
          else:
            lines.append(line)
      if removed:
        newfile = phs.GetFileName('phs')+'_fixed.phs'
        with open(newfile,'w') as g:
          for line in lines:
            g.write(line)
        phs_new=proc.out.AddCopy(phs)
        proc.out.SetFileToChild(phs_new,newfile,filetype='phs')
    if self.out.Get('model',typ=('partial+substr','partial'),filetype='pdb'):
      # "fix" the shelxe bug with substr. outputted in Z chain sometimes conflicting with protein Z chain
      # and similarly, z chain from shelxe may contain multiple protein chains with the same residue id
      # this can be removed once shelxe fixes the problem
      pdbset=proc.AddProg('pdbset')
      pdbset.inp.Set(self.out.Get('model',typ=('partial+substr','partial'),filetype='pdb'))
      pdbset.SetKey('renumber 1 1 TO 1000 chain Z')
      pdbset.SetKey('renumber 1 1 TO 9999 chain z')
      pdbset.Run()
      if self.out.Get('model',typ=('partial+substr'),filetype='pdb'):
        # "fix" the substructure in the output pdb file
        fixsubstr=proc.AddProcess('fixsubstrpdb')
        fixsubstr.inp.Set(pdbset.out.Get('model'))
        fixsubstr.Run()
      proc.programs.remove(pdbset)
