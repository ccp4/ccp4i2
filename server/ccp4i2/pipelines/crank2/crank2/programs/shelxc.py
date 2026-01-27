#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class shelxc(program):
  name="SHELXC"
  binary="shelxc"
  ccp4_parsing=True
  stat={}
  stat['data_set']   = common.stats(name='SHELX shortcut for the data set',multiple=True,
                              regexp=[r' Reflections read from '+r'(\S+)'+r" file"]) 
  stat['resol_bins'] = common.stats(name='resolution in bin',multiple=True,
                               regexp=[r" Resl.\s+Inf."+r'(?:\s+(\d+\.\d+))?'*15] ) 
  stat['dano_bins'] = common.stats(name='anomalous difference in bin',multiple=True,
                              regexp=[r" <d\"/sig>"+r'(?:\s+(\d+\.\d+))?'*15] ) 
  stat['Ioversigma'] = common.stats(name='Intensity over sigma in bin',multiple=True,
                              regexp=[r" <I/sig>"+r'(?:\s*(-?\d+\.\d?))?'*15] ) 
  stat['complete']  = common.stats(name='percentage completeness in bin',multiple=True,
                              regexp=[r" %Complete"+r'(?:\s+(\d+\.\d+))?'*15] ) 
  stat['resol_bins_cor'] = common.stats(name='resolution in bin for correlations',
                               regexp=[r" Resl.\s+Inf"+r'(?:\s+-\s+(\d+\.\d+))?'*15] ) 
  stat['correl']    = common.stats(name='signed anomalous difference correlation in bin',
                              regexp=[r"{0}/{1}"+r'(?:\s+(-?\d+\.\d+))?'*15] )
  stat['cc1/2'] = common.stats(name='half set correlation coefficient', multiple=True, 
                              regexp=[r" CC\(1/2\)"+r'(?:\s*(-?\d+\.\d?))?'*15])
  stat['resol'] = common.stats(name='data resolution', multiple=True, 
                              regexp=[r"highest resolution\s+(\S+)\s+Angstroms"])
  stat['version'] = common.stats(name='version', regexp=r"Version (\d+/\d+)")
  stat['increase_maxm_error'] = common.stats(name='increase MAXM', regexp=r' \*\*.+increase MAXM', accept_none=True)
  stat['corruption_error'] = common.stats(name='input file corrupted', regexp=r' \*\*.+corrupted at line', accept_none=True)

  references = ( "Sheldrick GM (2008) A short history of SHELX. " +
                 "Acta Cryst. A64, 112-122.", )

  def Init(self):
    self.SetArg('shelx')
    self.default_exper='RICE'

  def TreatInput(self):
    try:
      exper=self.process.GetVirtPar('target')
    except common.CrankError:
      # dummy option, just to get the shelxc logs/conversion
      exper=self.default_exper
    # get required data for the target; try I's first and then F's
    fi_opt=''
    obj,N=self.inp.GetTargetObjects(exper, inten=True, no_model=True, no_warn=True)
    if obj is None:
      obj,N=self.inp.GetTargetObjects(exper, inten=False, no_model=True)
      fi_opt='-f'
    f_conv = ['amplitudes']  if fi_opt  else []
    if obj is None:
      common.Error('The input objects for target {0} of {1} could not be retrieved'.format(exper,self.name))
    self.obj=obj
    # set cell and spacegroup
    if not self.IsKey('cell') or not self.IsKey('spag'):
      inp_obj = 'f+' if 'f+' in obj else 'f'
      if inp_obj=='f+':  obj['f-'].GetCell(self)  # just to make sure that cad merges F- and F+ (it suffers by a bug/feature of quitting if only one of F+/F- is input)
      if not self.IsKey('cell'):
        self.cell=obj[inp_obj].GetCell(self)
        self.SetKey('cell', ' '.join((str(c) for c in self.cell)))
      if not self.IsKey('spag'):
        self.spgr=obj[inp_obj].GetSpacegroup(self)
        self.SetKey('spag', self.spgr)
    # set atomtype
    if not self.IsKey('sfac') and 'mod' in obj and obj['mod'] and obj['mod'].atomtypes:
      self.SetKey('sfac', obj['mod'].GetAtomType())
    # SAD, MAD and SIRAS input
    # shelxc does not deal with long path strings so we use relative path at input
    if exper=='SIRAS' or exper=='SAD':
      self.SetKey(exper, (fi_opt, obj['f+'].GetFileName(('HKL','sca','hkl'),conv_parent=self,inp_cont=obj['f-'],conv_opts=f_conv)))
      self.CheckFilePathLength(exper,fi_opt)
    elif exper=='MAD':
      keys=['PEAK','INFL','HREM','LREM']
      mad_assign={}
      # try to see whether the assignment is defined in dataname or custom attributes
      for f in ['f'+str(n)+'+' for n in range(1,N+1)]:
        for key in keys:
          if key in [ c.upper() for c in obj[f].custom ] or key in obj[f].dname.upper():
            mad_assign[key] = f
      # if not then just define by order
      for f in ['f'+str(n)+'+' for n in range(1,N+1)]:
        if len(mad_assign)==len(keys):
          break
        if f not in mad_assign.values():
          key=[k for k in keys if k not in mad_assign][0]
          mad_assign[key] = f
      for key,val in mad_assign.items():
        self.SetKey(key, (fi_opt, obj[val].GetFileName(('HKL','sca','hkl'),conv_parent=self,inp_cont=obj[val[:-1]+'-'],conv_opts=f_conv)))
        self.CheckFilePathLength(key,fi_opt)
    else:
      # a "dummy" option to obtain shelx logfile etc
      self.SetKey('SAD', (fi_opt, obj['f'].GetFileName(('HKL','sca','hkl'),conv_parent=self,conv_opts=f_conv)))
      self.CheckFilePathLength('SAD',fi_opt)
    if 'fn' in obj:
      self.SetKey('NAT', (fi_opt, obj['fn'].GetFileName(('HKL','sca','hkl'),conv_parent=self,conv_opts=f_conv)))
      self.CheckFilePathLength('NAT',fi_opt)

  def CheckFilePathLength(self,key,fi_opt):
    max_length=75
    if len(self.GetKey(key)[1])>max_length:
      import shutil
      loc_file = os.path.basename(self.GetKey(key)[1])
      if len(loc_file)>max_length:
        loc_file = loc_file[:max_length-4] + loc_file[-4:]
      shutil.copy(self.GetKey(key)[1],loc_file)
      self.SetKey(key, (fi_opt, loc_file), keep_previous=False)


  def TreatParams(self):
    if self.process.nick in ('faest','substrdet') and self.process.GetParam('high_res_cutoff') and not self.GetKey('SHEL'):
      self.SetKey('SHEL', (999., self.process.GetParam('high_res_cutoff')))


  def DefineOutput(self):
    # the shelxc output filenames are hardcoded...
    project=self.arg_list[0]
    # SAD/SIRAS only output yet. MAD will need to be added later.
    o = None
    if   'fn' in self.obj:  o = self.obj['fn']
    elif  'f' in self.obj:  o = self.obj['f']
    elif 'f+' in self.obj:  o = self.obj['f+']
    if o: 
      self.out.AddNew('fsigf', project+'.hkl', filetype='hkl', typ='average', xname=o.xname, dname=o.dname, cell=o.cell, spgr=o.spgr)
    if 'f+' in self.obj:
      fa=self.out.AddNew('fsigf', project+'_fa.hkl', filetype='hkl', typ='fa', xname=self.obj['f+'].xname, dname=self.obj['f+'].dname, cell=self.obj['f+'].cell, spgr=self.obj['f+'].spgr)
      fp=self.out.AddNew('fsigf', project+'_sad.cif', filetype='cif', typ='plus', xname=self.obj['f+'].xname, dname=self.obj['f+'].dname, cell=self.obj['f+'].cell, spgr=self.obj['f+'].spgr)
      fm=self.out.AddNew('fsigf', project+'_sad.cif', filetype='cif', typ='minus', xname=self.obj['f+'].xname, dname=self.obj['f+'].dname, cell=self.obj['f+'].cell, spgr=self.obj['f+'].spgr)
      # this might need to be modified in case shelxc may output amplitudes
      fp.SetLabel( ['i','sigi'] ),  fm.SetLabel( ['i','sigi'] )
    ins=self.out.AddNew('datafile', project+'_fa.ins', typ='ins')


  def Run(self,*args,**kwargs):
    program.Run(self,*args,**kwargs)
    # shelx programs do not return errors...
    maxm=10
    while self.GetStat('increase_maxm_error') and maxm<500:
      if maxm==10:
        common.Warning('MAXM for SHELXC will be increased to fit data into memory.')
      else:
        maxm+=10
      self.SetKey('MAXM', maxm, keep_previous=False)
      program.Run(self,*args,**kwargs)
    if self.GetStat('increase_maxm_error') and maxm>500:
      common.Warning('Could not increase MAXM for SHELXC to fit data into memory.')
