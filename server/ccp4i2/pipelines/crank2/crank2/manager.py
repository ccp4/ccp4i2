import os,sys
import pkgutil,shutil,copy
from . import common,processes,ccp4i2crank
from .process import process, crvapi
from .program import program

class crank(process):
  name='CRANK2'
  ccp4i2=None
  jscofe=None
  # parameters
  supported_params = copy.copy(process.supported_params)
  supported_params['target'] = common.parameter( desc='Experiment/refinement target', typ=str, cap=True )
  supported_params['replace_met_mse'] = common.parameter( desc='Replace methionines by selenomethionines in input model(s) with SE substructure', typ=bool )
  supported_params['force_xdnames'] = common.parameter( desc='Replace the inputted xname/dname in the generated mtz files', typ=bool )
  # all processes
  supported_procs = [ n for _,n,_ in pkgutil.iter_modules( [os.path.dirname(processes.__file__)] ) ]
  # filehandle of the loggraphfile - LGInfo() will send output to this filehandle if not None
  loggraphfilehandle=None
  loggraph_actfilehandle=None
  x_rename, d_rename = {}, {}
  # ccp4 version is determined at the beginning at saved here (for loggraph printing)
  ccp4_version=None
  # peak and inflection point wavelengths for some atoms 
  # used to guess the peak/infl/remote assign. - this can be used as the last resort for shelxc
  # also used the other way around, ie use the theoretical values if dtype (peak/infl etc) was specified
  w_pi = { "SE": [(0.9793,0.97944)], "HG": [(1.0092,1.0098)], "AU": [(1.0400,1.0405),(0.9026,0.9028)], \
           "PT": [(1.0721,1.0725),(0.9339,0.9341)], "ZN": [(1.2836,1.2841)], "BR": [(0.9199,0.9201)] }
  # references
  references = ( 'Skubak P, Pannu NS (2013) Automatic protein structure solution from weak X-ray data.' + 
                 'Nature Communications 4, 2777.',
                 'Pannu NS, Waterreus WJ, Skubak P, Sikharulidze I, Abrahams JP and ' +
                 'de Graaff RAG (2011) Recent advances in the CRANK software suite ' +
                 'for experimental phasing. Acta Cryst. D67, 331-337.', )

  def TreatInOutPar(self, set_all_par=False, emulate=False, initialized=True):
    if self.GetParam('replace_met_mse') is None and not self.inp.Get('model',('substr','partial+substr')):
      for mod in self.inp.GetAll('model',has_atomtypes=True):
        if 'SE' in mod.GetAtomTypes():
          self.SetParam('replace_met_mse', True)
    # the target setting requires initialization.  all dependend on target should come after initialization!
    if initialized:
      # if S substr. than setting S subst. for native too by default (can be disabled eg by inputting empty native model without anom. atoms)
      if self.GetParam('target'):
        self.CheckSSADnative()
      else:
      # automatic determination of the experiment/target from the data.  
        for targ in ('MAD','SIRAS','SAD'):
          if self.inp.GetTargetObjects(targ,inten=None,no_warn=True)[0]:
            self.SetParam('target',targ)
            break
        else:
          common.Warning("Experiment type/target could not be determined from the input data.")
        # taking care of SAD/SIRAS subtleties
        if self.GetParam('target'):
          self.CheckSSADnative()
        if targ=='SIRAS':
          (obj,N),sfi=self.inp.GetTargetObjects(targ,inten=False),'f'
          if obj is None:
            (obj,N),sfi=self.inp.GetTargetObjects(targ,inten=True),'i'
          #!!!!! temporary till refmac SIRAS crash is fixed...
          #if 'phas' in [p.nick for p in self.processes]:
          #  self.SetParam('target','SAD')
          if obj and self.inp.Get('model',typ=('substr','partial+substr'),xname=obj['fn'].GetCrystalName()):
            substr_not_in_native_input=self.inp.Get(is_native=True,has_num_atoms=True) and self.inp.Get(is_native=True,has_num_atoms=True).exp_num_atoms==0
            if not substr_not_in_native_input:
              self.SetParam('target','SAD')
          elif obj and obj['f+'].GetFileName('mtz') and obj['fn'].GetFileName('mtz'):
            sft=self.AddProg('sftools')
            sft.SetKey('read', ['"'+obj['f+'].GetFileName('mtz')+'"','col','"'+obj['f+'].GetLabel(sfi)+'"'])
            sft.SetKey('read', ['"'+obj['fn'].GetFileName('mtz')+'"','col','"'+obj['fn'].GetLabel(sfi)+'"'])
            sft.SetKey('correl', ('col','1','2') ),   sft.SetKey('quit\nY')
            sft.Run()
            corr=sft.GetStat('correl')
            self.programs.remove(sft)
            if corr>0.999:
              self.SetParam('target','SAD')
        if not emulate and self.GetParam('target'):
          self.Info('{0} phasing experiment.'.format(self.GetParam('target')))
    process.TreatInOutPar(self,set_all_par)

  def CheckSSADnative(self):
    obj,N=self.inp.GetTargetObjects(self.GetParam('target'),inten=None)
    if obj and 'fn' in obj and 'mod' in obj and obj['mod'].GetAtomType()=='S' and (not 'modn' in obj or not obj['modn']):
      self.inp.AddNew('model',xname=obj['fn'].GetCrystalName(),atomtype1='S',typ='substr',custom='native')

  def RunPreprocess(self,emulate=False,*args,**kwargs):
    if hasattr(self,'logout') and not hasattr(self,'stdout_save'):
      self.Info('Output will be redirected to {0}'.format(self.logout))
      self.stdout_save = sys.stdout
      sys.stdout = open(self.logout,'w',1)
    global crvapi
    if hasattr(self,'disable_rvapi') and crvapi:
      self.Info('Disabling rvapi output.')
      crvapi = False
    if crvapi:
      if hasattr(self,'rvapi_uri_prefix'):
        crvapi.uri_prefix = self.rvapi_uri_prefix
      if not self.ccp4i2 and hasattr(self,'rvapi_document') and hasattr(self,'rvapi_viewer') and self.rvapi_viewer=='0':
        self.run.jscofe, self.jscofe = True, True
      crvapi.tree = False  if hasattr(self,'rvapi_no_tree')  else True
      crvapi.separate_steps = True  if hasattr(self,'rvapi_separate_steps')  else False
    self.TreatInOutPar(emulate=emulate,initialized=False)
    process.RunPreprocess(self,emulate=emulate,crvapi=crvapi,*args,**kwargs)
    self.FixInpData()
    self.DetermineUnknownModel(emulate=emulate)
    self.GetAnomScatCoefs()
    if not emulate:
      self.FixInpModels()

  def RunBody(self, emulate=False, set_all_par=False, *args, **kwargs):
    """Runs CRANK2 job."""
    self.run.inp = self.inp
    error=self.Initialize(self,emulate=emulate)
    self.run.TreatInOutPar(set_all_par=set_all_par, emulate=emulate)
    if error: 
      if not emulate:
        common.Error(error)
      else:
        common.Warning(error+' Stopping preliminary!')
    # the emulation goes through the TreatInOutPar of all steps without actually running them
    # if set_all_par is True then for all the steps all the parameters are set and printed
    for ip,p in enumerate(self.processes):
      #  at this point we are dealing with "dummy" preparation process instances
      #  (this is in order to write desired XML files by conversion from dummy instances)
      # first pass output from the previous step to input of the current step
      #if not emulate:
      self.PassPrevStepsOutput(p,ip,emulate)
      # set target for the process if needed
      self.SetTargetDefault(p,emulate=emulate)
      # pass previous phasing/dm jobs to the handdet process
      if p.nick=='handdet':
        self.PassToHandDet(p,ip)
      # write XML of process p and create its run instance by reading it
      p_run = self.SubProcess_WriteXML_CreateRun(p,ip,emulate=emulate)
      # this will assure that if eg cell params are found but no spacegroup then -d will return them
      if error and emulate:
        continue
      # any other general process initialization if needed
      self.Initialize(p_run)
      if not emulate:
        # run the p_run process
        self.RunSubProcess(p_run,ip)
        # pass the chosen hand mtz's/pdb's to all the subsequent steps
        if p.nick=='handdet':
          self.PassChosenHand(p_run,ip)
        # pass corrected data from afro
        if p.nick=='faest' and p.GetProg(supported=True) and p.GetProg(supported=True).nick=='afro':
          self.PassLocscl(p_run,ip)
      else:
        p_run.TreatInOutPar(set_all_par)
    # emulation cleanup
    if emulate:
      os.chdir(self.rundir)
      for p_run in self.run.processes:
        if os.path.isdir(p_run.rundir) and os.path.realpath(os.path.realpath(p_run.rundir)).startswith(self.rundir):
          shutil.rmtree(p_run.rundir)


  def CheckBinaries(self):
    # checks whether the binaries needed exist (and their versions if such a check is defined for program)
    # first check sftools (it is needed to figure out what other programs will be used)
    program.from_name('sftools',None).CheckBinary()
    # let us construct a copy of the crank used for the run and run TreatInOutPar() 
    # for each step so that the default programs are assigned
    self.emulate = self.GetProcessCopy()
    #crank = self.AddProcessCopy(self)
    #self.processes.remove(crank)
    self.emulate.run.check_binaries=True
    self.emulate.RunBody(emulate=True)
    checked=[]
    for p in self.emulate.run.GetFlatSubTree():
      if hasattr(p,'binary'):
        program.TreatParams(p)
        # make sure each program is included only once
        if p.binary not in checked:
          p.CheckBinary()
          checked.append(p.binary)

  def GetCCP4Version(self):
    # check/get ccp4 library version
    if not self.ccp4_version:
      bp3=program.from_name('bp3',None)
      bp3.SetRunDir('checkbinary')
      bp3.runname='{0}_ccp4_version'.format(self.name)
      bp3.ExternalRun( [bp3.binary, '-i'], [] )
      self.ccp4_version = bp3.GetStat('ccp4_version')
    return self.ccp4_version

  def GetVersion(self):
    version_file = os.path.join(os.path.dirname(common.__file__), 'VERSION')
    with open(version_file) as f:
      version = f.readline()
    return version

  def DetermineUnknownModel(self,emulate=False):
    for p in self.GetFlatSubTree():
      # try to determine type of 'unknown' pdb from user input
      model_to_fix=p.inp.Get('model',typ='unknown',filetype='pdb')
      if model_to_fix:
        at_obj=p.inp.Get('model',has_atomtypes=True)
        if not at_obj: at_obj=self.inp.Get('model',has_atomtypes=True)
        if at_obj:
          pdbcur=self.AddProg('pdbcur',propagate_out=False,propagate_inp=False)
          pdbcur.inp.Add(model_to_fix)
          for at in at_obj.GetAtomTypes():
            pdbcur.SetKey('delatom', "[{0}]".format(at))
            # some programs output eg Cu,Se and PDBCUR is case sensitive
            pdbcur.SetKey('delatom', "[{0}]".format(at[:-1]+at[-1].lower()))
          if self.GetParam('replace_met_mse') and 'S' not in at_obj.GetAtomTypes():
            pdbcur.SetKey('delatom', "(MET)/SD[S]".format(at))
          pdbcur.Run()
          anom_atoms_present = pdbcur.GetStat('atoms_deleted')
          pdbcur.inp.Set(pdbcur.out.model)
          pdbcur.SetKey('delatom', "[C]", keep_previous=False)
          pdbcur.Run()
          prot_atoms_present =  pdbcur.GetStat('atoms_deleted')
          if anom_atoms_present:
            model_to_fix.SetType('substr')
            if prot_atoms_present:
              model_to_fix.SetType('partial+substr')
          elif prot_atoms_present:
            model_to_fix.SetType('partial')
          else:
            common.Error('Neither substructure nor partial model present in file {0}.'.format(model_to_fix.GetFileName('pdb')))
          self.Info('Type of {0} model set to {1}.'.format(model_to_fix.GetFileName('pdb'), model_to_fix.GetType()))
          if emulate:
            shutil.rmtree(pdbcur.rundir)
          self.programs.remove(pdbcur)
        else:
          action = common.Warning if emulate else common.Error
          action('Could not classify {0} as substructure atom types not specified.'.format(model_to_fix.GetFileName('pdb')))

  def GuessDtype(self,dname,xname,atomtype=None,wl=None,safe_guess=False):
    # some heuristics to "guess" the dtype
    dtype_strings = (("peak",),("infl","inflection"),("hrem","hrm","high","high-remote"),("lrem","lrm","low","low-remote"))
    dtype, dnum, num = None, None, 0
    # safe guess - exact match of the dname or custom
    for ik,key in enumerate(dtype_strings):
      if dname in key or self.inp.Get('fsigf',xname=xname,dname=dname,custom=key):
        return key[0],ik
    if safe_guess:
      return None,None
    # try to guess from the dname
    for idt,dts in enumerate(dtype_strings):
      for dt in dts:
        if dt in dname.lower():
          num+=1
          if num==1:
            dtype, dnum = dts[0], idt
          else:
            dtype, dnum = None, None
            break
    # try to guess from the wavelength and atomtype
    if not dtype and atomtype in self.w_pi:
      edge_dist = [ abs(wl-edge[0]) for edge in self.w_pi[atomtype] ]
      ind = edge_dist.index(min(edge_dist))
      pi_dist = [ abs(wl-pi) for pi in self.w_pi[atomtype][ind] ]
      dnum = ind = pi_dist.index(min(pi_dist))
      dtype = "peak" if ind==0 else "infl"
      if pi_dist[ind]>0.0075:
        dtype = "hrem" if ind==0 else "lrem"
        dnum+=2
    return dtype, dnum

  def GetAnomScatCoefs(self):
    for p in self.GetFlatSubTree():
      # get anomalous scattering coefficients from wavelength if not supplied by user
      models_to_fix = p.inp.GetAll( 'model', has_atomtypes=True )
      for model_to_fix in models_to_fix:
        atomtype = model_to_fix.GetAtomType()
        xname=model_to_fix.GetCrystalName()
        # this only does one atom type (possible to generalize later of course)
        if None in [ sv  for v in model_to_fix.atomtypes[atomtype]  for sv in v[:3]] and (self.inp.Get('fsigf',col='f',xname=xname) or self.inp.Get('fsigf',col='i',xname=xname)):
          crossec=self.AddProg('crossec')
          crossec.inp.Set(model_to_fix)
          # we assign theoretical wavelength if wavelength is not known but dname/custom specifies the MAD "type"
          # in fact, we want to prefer wavelength from dname/custom for f',f'' generation of peak, infl
          if atomtype in self.w_pi:
            crossec.BackupAnyPars()
            try:
              crossec.TreatInput()
            except common.CrankError as e:
              if hasattr(sys,'exc_clear'): sys.exc_clear()
            crossec.RestoreAnyPars()
            for dn in crossec.dnames:
              if not crossec.wavels or dn not in list(zip(*crossec.wavels))[0]:
                dtype,idt = self.GuessDtype(dn,xname,safe_guess=True)
                if dtype:
                  ass_obj = self.inp.Get('fsigf',xname=xname,dname=dn)
                  ass_obj.wavel = self.w_pi[atomtype][0][idt]  if idt<2  else  self.w_pi[atomtype][0][idt-2]-0.004*(idt-2.5)
                  common.Warning("Wavelength for dataset {0} not known, using theoretical value of {1}".format(dn,ass_obj.wavel))
              elif dn in ('infl','inflection','peak'):
                ass_obj = crossec.inp.AddCopy( self.inp.Get('fsigf',xname=xname,dname=dn) )
                ass_obj.wavel = self.w_pi[atomtype][0][0]  if dn=='peak'  else  self.w_pi[atomtype][0][1]
                self.Info("Using theoretical wavelength of {1} for f',f'' {0} estimation".format(dn,ass_obj.wavel))
          try:
            crossec.Run()
          except common.CrankError as e:
            common.Warning("Could not determine f',f'' due to error: {0}".format(e))
            if hasattr(sys,'exc_clear'): sys.exc_clear()
          else:
            # saving the determined f',f'' in the input model for all datasets where it was not set before
            for i,(dn,wl) in enumerate(crossec.wavels):
              fp,fpp,d,a=model_to_fix.Getfpfpp(atomtype,dn)
              if not fp and not fpp:
                # some heuristics to "guess" the dtype
                dtype,idt = self.GuessDtype(dn,xname,atomtype,wl)
                model_to_fix.SetAtomType( atomtype, *(crossec.GetStat('anom_scat_coef',atomtype)[i]), dname=dn, guessed_dtype=dtype )
                fp,fpp,d,a = model_to_fix.Getfpfpp(atomtype,dn)
                if atomtype.upper()=='SE' and dtype == 'peak' and fpp and fpp<1.5: # workaround for crossec providing wrong values for SE peak if wavelength is at the edge
                  model_to_fix.SetAtomType( atomtype, fp=-8.0, fpp=4.0, dname=dn, guessed_dtype=dtype )
                  self.Info("Will use f''=4 for SE peak.")
                elif atomtype.upper()=='SE' and dtype in ('infl','inflection') and fpp and fpp<1.0: # workaround for crossec providing wrong values for SE inflection point
                  model_to_fix.SetAtomType( atomtype, fp=-9.0, fpp=2.0, dname=dn, guessed_dtype=dtype )
                  self.Info("Will use f''=2 for SE inflection point.")
                elif [p for p in self.GetFlatSubTree() if p.nick in ('refatompick','comb_phdmmb')]:  # adjust for SAD - numberical reasons (too small fpp may lead to refmac crashes but also in general, it's hardly possible to solve by SAD otherwise)
                  if atomtype.upper()=='SE' and fpp and fpp<1.0:
                    model_to_fix.SetAtomType( atomtype, fp=fp, fpp=1.0, dname=dn, guessed_dtype=dtype )
                    self.Info("Will use f''=1")
                  elif fpp and fpp<0.5:
                    model_to_fix.SetAtomType( atomtype, fp=fp, fpp=0.5, dname=dn, guessed_dtype=dtype )
                    self.Info("Will use f''=0.5")
              self.Info("Determined f'={} and f''={} for data {}".format(fp,fpp,d))
          self.programs.remove(crossec)
      ok=False if models_to_fix and [p for p in self.GetFlatSubTree() if p.nick in ('refatompick','comb_phdmmb')] else True
      for model_to_fix in models_to_fix:  # report an issue if we have no fp,fpp for the SAD function
        if not ok:
          for at in model_to_fix.GetAtomTypes():
            for data in model_to_fix.atomtypes[at]:
              if data[0] is not None and data[1] is not None: #fp and fpp
                ok=True
                break
      if not ok:
        common.Error('Anomalous scattering coeficients could not be determined.  Please input them or the wavelength.')

  def FixInpModels(self):
    for p in self.GetFlatSubTree():
      # try fixing substructure pdb from user input
      model_to_fix=p.inp.Get('model',typ='substr',filetype='pdb')
      if model_to_fix:
        fixsub=self.AddProcess('fixsubstrpdb',propagate_out=False,propagate_inp=False)
        fixsub.inp.Add(model_to_fix)
        fixsub.AddProg('pdbcur')
        fixsub.GetProg().outfilename['pdb'] = "fixed_"+p.nick+".pdb"
        fixsub.Run()
        p.inp.model.remove(model_to_fix)
        p.inp.Add(fixsub.out.Get('model'), propagate=False)
        self.processes.remove(fixsub)
      # estimate number of substr. atoms
      model_to_fix=p.inp.Get('model',typ=('substr','partial+substr'),filetype='pdb')
      if model_to_fix and not model_to_fix.exp_num_atoms:
        model_to_fix.GuessNumSubstrAtomsFromSeq(self)
      if self.GetParam('replace_met_mse'):
        for mod in p.inp.GetAll('model',typ=('partial','partial+substr'),has_atomtypes=True):
          if 'SE' in mod.GetAtomTypes():
            pdbcur=self.AddProg('pdbcur',propagate_out=False,propagate_inp=False)
            pdbcur.inp.Add(mod)
            pdbcur.SetKey('renresidue', "(MET) 'MSE'")
            pdbcur.SetKey('renatom', "SD[S] 'SE'")
            pdbcur.SetKey('renelement', "SE[S] 'SE'")
            pdbcur.outfilename['pdb'] = "mse_repl_"+p.nick+".pdb"
            pdbcur.Run()
            p.inp.model.remove(mod)
            p.inp.Add(pdbcur.out.Get('model'), propagate=False)
            self.programs.remove(pdbcur)

  def FixInpData(self):
    # makes sure that xname/dname assignments match those in the mtz
    objs = {}
    for o in self.inp.GetAll('fsigf',filetype='mtz',try_convert=False):
      objs.setdefault(o.GetFileName('mtz'),[]).append(o)
    if objs:
      mergemtz=self.GetProcess('mergemtz')
      if not mergemtz:
        mergemtz=self.AddProcess('mergemtz',propagate_out=False,propagate_inp=False)
      mergemtz.SetParam('cons_check_suff')
      simul = not self.GetParam('force_xdnames')
      issue, self.x_rename, self.d_rename, save_xd = "", {}, {}, []
      for i,f in enumerate(objs):
        cons = None  if f not in mergemtz.xdrename  else not mergemtz.xdrename[f]
        mergemtz.inp.Set(objs[f])
        mergemtz.runname = 'merge'+str(i)
        mergemtz.Run(outmtz=os.path.basename(f).rsplit('.',1)[0]+'_cad{0}.mtz'.format(i),simulate=simul,consistent=cons)
        if not simul:
          for o in objs[f]:
            self.inp.SetFileToChild(o,mergemtz.out.Get(filetype='mtz').GetFileName('mtz'))
        else:
          substr_not_in_native_input=self.inp.Get(xname='native',has_num_atoms=True) and self.inp.Get(xname='native',has_num_atoms=True).exp_num_atoms==0
          for (xmtz,dmtz),(xinp,dinp) in mergemtz.xdrename[f].items():
            save_xd.append( (xmtz,dmtz,xinp,dinp) )
            if xinp not in self.x_rename:
              self.x_rename[xinp]=xmtz
            elif self.x_rename[xinp]!=xmtz:
              issue='Inputted xname {0} corresponds to multiple mtz xnames ({1},{2}). '.format(xinp,xmtz,self.x_rename[xinp])
            if dinp not in self.d_rename:
              self.d_rename[dinp]=dmtz
            elif self.d_rename[dinp]!=dmtz:
              issue='Inputted dname {0} corresponds to multiple mtz dnames ({1},{2}). '.format(dinp,dmtz,self.d_rename[dinp])
            for di in [ di  for (xm,dm,xi,di) in save_xd  if xm==xmtz and dm==dmtz and dinp!=di ]:
              issue='MTZ xname,dname ({2},{3}) corresponds to multiple inputted dnames {0},{1}. '.format(di,dinp,xmtz,dmtz)
            if substr_not_in_native_input:
              xinp_other=[i for i,m in self.x_rename.items() if m==xmtz and i!=xinp]
              xinp_multi=[i[0] for i in mergemtz.xdrename_multiple[f][(xmtz,dmtz)] ]  if (xmtz,dmtz) in mergemtz.xdrename_multiple[f] else []
              if (xinp_other and (xinp=='native' or 'native' in xinp_other)) or \
                 (xinp_multi and 'native' in xinp_multi):
                issue='Crystal mtz name {0} is the same for native and derivative but it is specified that native contains no substructure. '.format(xmtz)
            if issue:
              if self.IsParam('force_xdnames'):
                common.Error(issue)
              else:
                common.Warning(issue+' Forcing all the input xnames+dnames!')
                self.SetParam('force_xdnames',True)
                self.FixInpData()
                return
      self.processes.remove(mergemtz)
      if simul:
        # if we do not adjust xname/dname in mtz then we have to adjust xname/dname in all the input objects
        for p in [self,]+self.GetAllSubPrograms()+self.GetAllSubProcesses():
          for inp_obj in p.inp.GetAll():
            if 'xname' in inp_obj.supported_attributes:
              if inp_obj.GetCrystalName() in self.x_rename:
                if inp_obj.GetCrystalName()=='native' and 'custom' in inp_obj.supported_attributes:
                  inp_obj.custom.append(inp_obj.GetCrystalName())
                inp_obj.SetCrystalName(self.x_rename[inp_obj.GetCrystalName()])
            if 'dname' in inp_obj.supported_attributes:
              if inp_obj.GetDataName() in self.d_rename:
                if inp_obj.GetDataName() in ('peak','infl','inflection','hrem','high-remote','lrem','low-remote') and 'custom' in inp_obj.supported_attributes:
                  inp_obj.custom.append(inp_obj.GetDataName())
                inp_obj.SetDataName(self.d_rename[inp_obj.GetDataName()])
            if 'd_name' in inp_obj.supported_attributes:
              for at,vs in inp_obj.GetAtomTypes().items():
                for v in vs:
                  if v[2] in self.d_rename:
                    inp_obj.SetAtomType( at,v[0],v[1],self.d_rename[v[2]] )
    # replace sigma_i < 0.0 by MNF  (do we want to do this for sigma_f as well?)
    objs = {}
    for o in self.inp.GetAll('fsigf',col='sigi',typ=('plus','minus','average'),filetype='mtz'):
      objs.setdefault(o.GetFileName('mtz'),[]).append(o)
    if objs:
      for i,f in enumerate(objs):
        sft=self.AddProg('sftools')
        sft.SetKey('read ("{0}", "\nY")'.format(f))
        sft.SetKey('quit\nY')
        sft.Run()
        if sft.GetStat('zero_sigma',accept_none=True):
          outmtz = objs[f][0].GetFileName('mtz',trim_path=True).rsplit('.',1)[0]+'_sft{0}.mtz'.format(i)
          for o in objs[f]:
            sft.SetKey('absent', ('col',o.GetLabel('i'),o.GetLabel('sigi'),'if','col',o.GetLabel('sigi'),'<=','0.0'), insert_index=-1)
            self.inp.SetFileToChild(o,os.path.join(sft.rundir,outmtz))
          sft.SetKey('write', outmtz, insert_index=-1)
          sft.Run()
        self.programs.remove(sft)
    # check for an empty sequence
    if self.inp.Get('sequence') and not self.inp.Get('sequence').GetSequenceString():
      common.Error('Check the input sequence file - it appears to be empty.')


  def SetTargetDefault(self, proc, emulate=False):
    # exclude processes that never need/support targets first
    if proc.nick not in ['phas','faest','substrdet','handdet','dmfull','phcomb','ref','comb_phdmmb','mbref','refatompick']:
      return
    if proc.nick=='dmfull' and self.GetProg(supported=True):
      return
    crank_targ = self.GetParam('target')  if self.GetParam('target')  else None
    if not crank_targ and hasattr(self,'run'):  crank_targ = self.run.GetParam('target')
    # if overall target/experiment is set
    if crank_targ and not proc.IsParam('target'):
      # now set the default target for the current process
      if proc.nick in ['phas','faest','substrdet'] or crank_targ=='RICE' or \
         (crank_targ=='SAD' and ((proc.nick!='ref' and proc.nick!='mbref') or \
          (not proc.inp.Get('fsigf',is_native=True) and self.inp.Get('model',has_fpp=True)) )) or \
         (crank_targ=='SIRAS' and (proc.nick=='dmfull' or proc.nick=='phcomb')):
         ### should disable SIRAS for DM (line above)? slower and suboptimal results eg for toxd,x15
        proc.SetParam('target', crank_targ)
      else:
        proc.SetParam('target', 'MLHL')
    # if target/experiment is not set
    if not crank_targ and not proc.IsParam('target'):
      if proc.nick in ['ref','phcomb','mbref']:
        proc.SetParam('target', 'RICE')
      elif not emulate:
        common.Error('No target specified for process of {0}'.format(proc.name))


  def PassToHandDet(self, handdet, ip):
    tobepassed=("phas","dmfull","refatompick")
    for tbp in tobepassed:
      if not handdet.GetParam('ignore_passed_phas') or tbp=='dmfull':
        for pp in [ i for i in range(ip) if self.processes[i].nick==tbp ]:
          handdet.processes.append(self.processes[pp])
          for o in self.run.processes[pp].out.GetAll(stored_order=True):
            handdet.processes[-1].out.AddCopy(o)
        #handdet.processes[-1].programs = self.run.processes[i]
      #for pp in [ self.run.processes[i] for i in range(ip) if self.run.processes[i].nick==tbp ]:
      #  handdet.processes.append(pp)

  def PassChosenHand(self, handdet_run, ip):
    for pp in [ self.processes[i] for i in range(ip+1,len(self.processes)) ]:
      # we could invert any previously inputted models to subsequent steps
      # but I would need to create createotherhand process to it this way.
      #for mod in pp.inp.GetAll('model'):
      #  handdet.CreateHand2Model(mod, pp.inp.Get('fsigf',xname=mod.xname).GetSpacegroup(handdet))
      #pp.inp.SetCopy()
      # so we simply replace all mtz's instead - if something else is needed it should be added/adjusted here
      pp.inp.SetCopy(handdet_run.out.fsigf,propagate=True)
      if hasattr(pp,'inp2'):
        pp.inp2.SetCopy(handdet_run.out2.fsigf,propagate=True)

  def PassLocscl(self, afro_run, ip):
    for pp in [ self.processes[i] for i in range(ip+1,len(self.processes)) ]:
      for fsf in afro_run.out.fsigf:
        pp.inp.AddCopy(fsf,propagate=True)

  def SubProcess_WriteXML_CreateRun(self,p,ip,emulate=False):
    # write the process input xml
    p_xml=os.path.join(self.rundir, p.nick+'.xml')
    crank_p=copy.copy(self)
    crank_p.processes=[p,]
#    if emulate:
#      p_xml = crank_p.Data2XML()
#    else:
#      common.WriteXML(crank_p, p_xml)
    common.WriteXML(crank_p, p_xml)
    # "undummify" by re-creating the process from the written xml
    # pass extra attrs set by user from parse as eg disable_mtz_label_prefix needs to be set at the process creation
    from . import parse
    extra_attrs = { attr:getattr(self,attr) for attr in parse.parse(dummy=True).share_with_crank if hasattr(self,attr) }
    crank_p=self.from_xml(p_xml, rundir=self.rundir, extra_attrs=extra_attrs)
    p_run = crank_p.processes[0]
    self.run.processes[ip] = p_run
    p_run.parent_process=self.run
    # initialize the 'actual' crank (meta)attributes not saved in xml
    if ip==0:
      self.run.prep=self
      for at in ('logfilehandle','loggraphfilehandle','opened_log','opened_loggraph','ccp4_version',
                 'rv_report','emulate','d_rename','x_rename') + \
                parse.parse(dummy=True).share_with_crank:
        if hasattr(self,at):
          setattr(self.run, at, getattr(self,at))
    p_run.SetRunDir(os.path.join(self.rundir,str(ip)+'-'+p.nick))
    return p_run

  def Initialize(self,p_run, emulate=False):
    # any other general process initialization
    #convert I+/I- to I average in case no I average was supplied
    for op in p_run.inp.GetAll('fsigf',col='i',typ='plus',filetype=('sca','hkl','HKL'),try_convert=False):
      if not p_run.inp.Get('fsigf',typ='average',col='i',xname=op.GetCrystalName(),dname=op.GetDataName()):
        o_new=p_run.inp.AddNew('fsigf',op.GetFileName(),typ='average',filetype=set(op.GetFileTypes()).intersection(('sca','hkl','HKL')).pop(), \
                 xname=op.GetCrystalName(),dname=op.GetDataName(),cell=op.cell,spgr=op.spgr,resol=op.resol,propagate=False)
        o_new.SetLabel( ['i','sigi'] )
    num=0
    for op in p_run.inp.GetAll('fsigf',col='i',typ='plus',filetype='mtz',try_convert=False):
      if not p_run.inp.Get('fsigf',typ='average',col='i',xname=op.GetCrystalName(),dname=op.GetDataName()):
        om = p_run.inp.Get('fsigf',typ='minus',col='i',xname=op.GetCrystalName(),dname=op.GetDataName())
        if om:
          num+=1
          sft=p_run.AddProg('sftools',propagate_out=False)
          sft.runname+='_I'+str(num)
          sft.SetRunDir()
          outname = op.GetFileName(trim_path=True).rsplit('.',1)[0] + '_Isft{}.mtz'.format(num)
          o_new = p_run.out.AddNew('fsigf',os.path.join(sft.rundir,outname),typ='average',xname=op.GetCrystalName(),dname=op.GetDataName(),propagate=False)
          o_new.SetLabel( ['i','sigi'] )
          o_new.custom.extend([c for c in om.custom if c in ('native','hrem','lrem','peak','infl')])
          i, si = o_new.GetLabel('i'), o_new.GetLabel('sigi')
          sft.SetKey('read', ['"'+op.GetFileName('mtz')+'"','col','"'+op.GetLabel('i')+'"','"'+op.GetLabel('sigi')+'"'] )
          sft.SetKey('read', ['"'+om.GetFileName('mtz')+'"','col','"'+om.GetLabel('i')+'"','"'+om.GetLabel('sigi')+'"'] )
          if op.GetLabel('f') and op.GetLabel('sigf') and om.GetLabel('f') and om.GetLabel('sigf'):
            sft.SetKey('read', ['"'+op.GetFileName('mtz')+'"','col','"'+op.GetLabel('f')+'"','"'+op.GetLabel('sigf')+'"'] )
            sft.SetKey('read', ['"'+om.GetFileName('mtz')+'"','col','"'+om.GetLabel('f')+'"','"'+om.GetLabel('sigf')+'"'] )
          sft.SetKey('calc', ['J','col',i,'=','col','1','col','3','+','2','/'] )
          sft.SetKey('calc', ['Q','col',si,'=','col','2','2','**','col','4','2','**','+','0.5','**','2','/'])
          sft.SetKey('select','centro'),  sft.SetKey('calc', ['Q','col',si,'=','col','2'])
          sft.SetKey('select','all'), sft.SetKey('select', ['col','1']),  sft.SetKey('select', ['not','col','3'])
          sft.SetKey('calc', ['J','col',i,'=','col','1']), sft.SetKey('calc', ['Q','col',si,'=','col','2'])
          sft.SetKey('select','all'), sft.SetKey('select', ['col','3']),  sft.SetKey('select', ['not','col','1'])
          sft.SetKey('calc', ['J','col',i,'=','col','3']), sft.SetKey('calc', ['Q','col',si,'=','col','4'])
          sft.SetKey('select','all'), sft.SetKey('write',outname), sft.SetKey('Y'), sft.SetKey('quit'), sft.SetKey('Y')
          sft.Run(check_bin=True)
          p_run.inp.Add(o_new)
          op_new, om_new = p_run.inp.AddCopy(op), p_run.inp.AddCopy(om)
          #op_new.InitLabels(),  om_new.InitLabels()
          #op_new.SetLabel( ['i','sigi'], [op.GetLabel('i'),op.GetLabel('sigi')] )
          #om_new.SetLabel( ['i','sigi'], [om.GetLabel('i'),om.GetLabel('sigi')] )
          self.inp.SetFileToChild(op_new, os.path.join(sft.rundir,outname), 'mtz')
          self.inp.SetFileToChild(om_new, os.path.join(sft.rundir,outname), 'mtz')
          p_run.programs.remove(sft)
    # try to obtain basic info about datasets (doing it only once per each xname+dname)
    cell_use, spgr_use, todo_cellspg = None, None, []
    for xn in p_run.inp.GetCrystalsList():
      for dn in p_run.inp.GetDataList():
        cell, spgr, datasets = None, None, p_run.inp.GetAll('fsigf',xname=xn,dname=dn)
        for o in datasets:
          if not cell:
            cell = o.GetCell(self,accept_none=True)
          if cell:
            o.cell, cell_use = cell, cell
          else:
            todo_cellspg.append(o)
          if not spgr or o.spgr or o.spgr_num:  # the last two make sure to handle eg spgr change at handdet
            spgr = o.GetSpacegroup(self,accept_none=True)
          if spgr:
            o.spgr, spgr_use = spgr, spgr
          else:
            todo_cellspg.append(o)
        if not cell and datasets:
          common.Warning('Cell parameters for crystal {0}, data {1} not known, will try to assign it cell from other dataset if possible.'.format(xn,dn))
        if not spgr and datasets:
          common.Warning('Spacegroup for crystal {0}, data {1} not known, will try to assign spacegroup from other dataset if possible.'.format(xn,dn))
    # take native cell+spacegroup from derivative if it was not found or vice versa
    for o in todo_cellspg:
      if not o.GetCell(self,accept_none=True):  o.cell = cell_use
      if not o.GetSpacegroup(self,accept_none=True):  o.spgr = spgr_use
    #convert I's to F's in case no F was supplied
    num=0
    wilson_scale={}
    if (self.processes and not [1 for p in self.processes if p.nick=='phdmmb']) or p_run.nick in ('mbref','ref'): # for the shelx pipeline, only convert to F's if needed in the mbref step (George's request)
     for o in p_run.inp.GetAll('fsigf',col='i',typ=('plus','average'),filetype='mtz'):
      if not p_run.inp.Get('fsigf',col='f',typ=o.GetType(),xname=o.GetCrystalName(),dname=o.GetDataName(),filetype='mtz') and \
         (not o.GetType()=='average' or not p_run.inp.Get('fsigf',col='i',typ='plus',xname=o.GetCrystalName(),dname=o.GetDataName(),filetype='mtz')):
        num+=1
# sftools conversion disabled as of now - it seems truncate works better at least for thio (?)
#        outname = o.GetFileName(trim_path=True).rsplit('.',1)[0] + '_Fsft'+str(num)+'.mtz'
#        o.SetLabel( ['f','sigf'] )
#        lab_types = ('G','L')  if o.GetType()=='plus'  else  ('F','Q')
#        sft=p_run.AddProg('sftools',propagate_out=False)
#        sft.runname+='_F'+str(num)
#        sft.SetRunDir()
#        sft.SetKey('read', ['"'+o.GetFileName('mtz')+'"','col','"'+o.GetLabel('i')+'"','"'+o.GetLabel('sigi')+'"'] )
#        sft.SetKey('i2f',['col','"'+o.GetLabel('i')+'"','"'+o.GetLabel('sigi')+'"'])
#        sft.SetKey('set',['label','col','3','4']), sft.SetKey(o.GetLabel('f')), sft.SetKey(o.GetLabel('sigf'))
#        sft.SetKey('set',['type','col','3','4']), sft.SetKey(lab_types[0]), sft.SetKey(lab_types[1])
#        if o.GetType()=='plus':
#          om = p_run.inp.Get('fsigf',typ='minus',col='i',xname=o.GetCrystalName(),dname=o.GetDataName())
#          om.SetLabel( ['f','sigf'] )
#          sft.SetKey('read', ['"'+om.GetFileName('mtz')+'"','col','"'+om.GetLabel('i')+'"','"'+om.GetLabel('sigi')+'"'] )
#          sft.SetKey('i2f',['col','"'+om.GetLabel('i')+'"','"'+om.GetLabel('sigi')+'"'])
#          sft.SetKey('set',['label','col','7','8']), sft.SetKey(om.GetLabel('f')), sft.SetKey(om.GetLabel('sigf'))
#          sft.SetKey('set',['type','col','7','8']), sft.SetKey(lab_types[0]), sft.SetKey(lab_types[1])
#        sft.SetKey('write',outname), sft.SetKey('Y'), sft.SetKey('quit'), sft.SetKey('Y')
#        sft.Run(check_bin=True)
#        self.inp.SetFileToChild(o, os.path.join(sft.rundir,outname), 'mtz')
#        if o.GetType()=='plus':
#          self.inp.SetFileToChild(om, os.path.join(sft.rundir,outname), 'mtz')
#        p_run.programs.remove(sft)
        # now put on rough absolute scale
#        if not emulate:
#          if o.GetCrystalName() not in wilson_scale: 
#            wilson = p_run.AddProg('wilson')
#            wilson.inp.Add(o)
#            wilson.Run()
#            wilson_scale[o.GetCrystalName()] = math.sqrt(wilson.GetStat('wilson_scale'))
#            p_run.programs.remove(wilson)
#          wsc=wilson_scale[o.GetCrystalName()]
#          outname = o.GetFileName(trim_path=True).rsplit('.',1)[0] + 'sc.mtz'
#          sft=p_run.AddProg('sftools',propagate_out=False)
#          sft.runname+='_Fsc'+str(num)
#          sft.SetRunDir()
#          sft.SetKey('read', '"'+o.GetFileName('mtz')+'"' )
#          sft.SetKey('calc',['col',o.GetLabel('f'),'=','col',o.GetLabel('f'),wsc,'*'])
#          sft.SetKey('calc',['col',o.GetLabel('sigf'),'=','col',o.GetLabel('sigf'),wsc,'*'])
#          if o.GetType()=='plus':
#            sft.SetKey('calc',['col',om.GetLabel('f'),'=','col',om.GetLabel('f'),wsc,'*'])
#            sft.SetKey('calc',['col',om.GetLabel('sigf'),'=','col',om.GetLabel('sigf'),wsc,'*'])
#          sft.SetKey('write',outname), sft.SetKey('Y'), sft.SetKey('quit'), sft.SetKey('Y')
#          sft.Run(check_bin=True)
#          self.inp.SetFileToChild(o, os.path.join(sft.rundir,outname), 'mtz')
#          if o.GetType()=='plus':
#            self.inp.SetFileToChild(om, os.path.join(sft.rundir,outname), 'mtz')
#          p_run.programs.remove(sft)

        orig_i = [copy.deepcopy(o),]
        om,oa=None,None
        if o.GetType()=='plus':
          om = p_run.inp.Get('fsigf',typ='minus',col='i',xname=o.GetCrystalName(),dname=o.GetDataName(),filetype='mtz')
          oa = p_run.inp.Get('fsigf',typ='average',col='i',xname=o.GetCrystalName(),dname=o.GetDataName(),filetype='mtz')
          if om:  orig_i.append(copy.deepcopy(om))
          if oa:  orig_i.append(copy.deepcopy(oa))
        ctrun=p_run.AddProg('ctruncate',propagate_out=False)
        ctrun.inp.Set(orig_i)
        o_base = o.GetFileName(trim_path=True).rsplit('.',1)[0]
        ctrun.outfilename['mtz'] = o_base+'_'+ctrun.name[:3]+str(num)+'.mtz'
        try:
          ctrun.Run(check_bin=True,return_to_orig_filenames=False)
        except ctrun.ProgramRunError as e:
          if hasattr(sys,'exc_clear'): sys.exc_clear()
          ctrun.out.ClearAll()
          trun=p_run.AddProg('truncate',propagate_out=False)
          trun.inp.Set(o)
          if om: trun.inp.Add(om)
          if oa: trun.inp.Add(oa)
          trun.runname+='_'+str(num)
          trun.outfilename['mtz'] = o_base+'_'+trun.name[:3]+str(num)+'.mtz'
          self.Info('{0} failed with error message: {1}'.format(ctrun.name,e))
          self.Info('Will try {0} instead.'.format(trun.name))
          trun.Run(check_bin=True)
          for o in trun.out.GetAll('fsigf'):
            p_run.inp.Add(o)
          p_run.programs.remove(trun)
        else:
          for o in ctrun.out.GetAll('fsigf'):
            p_run.inp.Add(o)
        p_run.programs.remove(ctrun)
    #convert F+/F- to F average (+FD) in case no F average was not supplied
    for o in p_run.inp.GetAll('fsigf',col='f',typ='plus',filetype='mtz'):
      if not p_run.inp.Get('fsigf',typ='average',col='f',xname=o.GetCrystalName(),dname=o.GetDataName()):
        mtzmadmod=p_run.AddProg('mtzMADmod',propagate_out=False)
        o_base = o.GetFileName(trim_path=True).rsplit('.',1)[0]
        mtzmadmod.outfilename['mtz'] = o_base+'_'+mtzmadmod.name[:7]+'.mtz'
        mtzmadmod.Run(check_bin=True)
        for out in mtzmadmod.out.GetAll(filetype='mtz',typ='average'):
          p_run.inp.Add( out, ind=p_run.inp.fsigf.index(p_run.inp.Get('fsigf',typ='plus',dname=out.dname,xname=out.xname)), propagate=False )
        p_run.programs.remove(mtzmadmod)
    # scale native to derivative
    if p_run.nick=='crank' and not emulate and p_run.inp.Get('fsigf',col='sigf',typ='average',filetype='mtz',is_native=True):
      scaleit=p_run.AddProg('scaleit', propagate_out=False)
      scaleit.revert_if_one=True
      #scaleit.SetKey('refi','scale')
      scaleit.Run(return_to_orig_filenames=False)
      p_run.programs.remove(scaleit)
    # setting solvent content and/or num. of monomers
    if (not self.inp.Get(has_solvent_content=True) or not self.inp.Get(has_monomers_asym=True)) and \
        self.inp.Get('fsigf') and (self.inp.Get('sequence') or self.inp.Get(has_residues_mon=True)) and \
        cell_use and spgr_use:
      matthews=p_run.AddProcess('matthews', propagate_out=False)
      #matthews.no_reporting = True
      matthews.Run()
      p_run.processes.remove(matthews)
    else:
      if self.inp.Get('sequence') and not self.inp.Get(has_residues_mon=True):
        self.inp.Get('sequence').residues_mon=len(self.inp.Get('sequence').GetSequenceString())
      if self.inp.Get(has_solvent_content=True) and self.inp.Get('sequence') and not self.inp.Get('sequence',has_solvent_content=True):
        self.inp.Get('sequence').solvent_content=self.inp.Get(has_solvent_content=True).solvent_content
    if self.inp.Get('sequence') and not self.inp.Get('sequence').seq_monomers:
      self.inp.Get('sequence').GetSequenceString()
    if not cell_use or not spgr_use:
      if not cell_use: error='No cell could be determined!'
      if not spgr_use: error='No spacegroup could be determined!'
      return error
    # get resolution for native only
    resol = None
    for o in p_run.inp.GetAll('fsigf',xname='native',filetype='mtz',try_convert=False)+p_run.inp.GetAll('fsigf',custom='native',filetype='mtz',try_convert=False):
      if not resol:
        resol = o.GetResolution(self,accept_none=True)
      else:
        o.resol = resol

  def RunSubProcess(self,p_run,ip):
    # skip ref if nothing was obtained
    prev_proc=self.run.processes[ip-1]
    if prev_proc.nick in ('mbref','comb_phdmmb') and prev_proc.mb_res_all and prev_proc.mb_res_all[-1][1]==0:
      common.Warning('Skipping process of {0} as no residues were obtained in building.'.format(p_run.name))
      return
    # register process with ccp4i2 if run from ccp4i2 and run the process
    if self.run.ccp4i2:
      ccp4i2crank.RegisterProcessToCCP4i2(self.run.ccp4i2, p_run)
    # run the current process
    error,nosuccess=False,False
    try:
      p_run.Run(from_ccp4i2=self.run.ccp4i2)
    except common.Unsuccessful as e:
      nosuccess=e
      if not self.no_reporting and crvapi:
        summary = self.rv_body.Section("summary","Job Summary",row=0)
        summary.Text('<i>The pipeline stopped preliminary at {}:</i>'.format(p_run.name))
        summary.Text(str(e), flush=True)
    except Exception as e:
      error=e
    if self.run.ccp4i2:
      # register process output with ccp4i2
      ccp4i2crank.RegisterOutputToCCP4i2(p_run,error,nosuccess)
    if nosuccess or error:
      try: #python2
        raise
      except RuntimeError: #python3
        raise error if error else nosuccess #from None

  def PassPrevStepsOutput(self,p,ip,emulate=False):
      # (indirectly) pass shelxc output to the subsequent steps if needed
      prev_proc=self.run.processes[ip-1]
      if hasattr(prev_proc,'out2'):
        # create hand2 inp
        if not hasattr(p,'inp2'):
          from . import inout
          p.inp2 = inout.input_output(is_output=False,parent=p)
          for o in p.inp.GetAll(stored_order=True):
            p.inp2.AddCopy(o,propagate=True)
      if 'shelxc' in [r.nick for r in prev_proc.programs] and not prev_proc.IsNonFalseVirtPar('no_output_to_next_step') and \
         [pr for i,pr in enumerate(self.processes) if i>ip and not [prg for prg in pr.programs if 'shelx' in pr.name]]:
        for hkl in prev_proc.out.GetAll('fsigf',filetype=('hkl','cif'),try_convert=False):
          if not self.inp.Get('fsigf',filetype='mtz',typ=hkl.GetType(),try_convert=False,xname=hkl.GetCrystalName(),dname=hkl.GetDataName()):
            mtz_obj = prev_proc.out.Get('fsigf',filetype='mtz',typ=hkl.GetType(),xname=hkl.GetCrystalName(),dname=hkl.GetDataName())
            if mtz_obj:
              self.inp.Add( mtz_obj )
            else:
              common.Warning('Could not convert shelxc output {0} to mtz!'.format(hkl.GetFileName()))
        self.Initialize(self)
      # pass "unbiased" phases from phasing to building
      if p.nick in ('mbref','comb_phdmmb') and p.GetParam('target')!='SAD':
        handdet=next((r for ir,r in enumerate(self.run.processes) if ir<ip and r.nick in ('handdet') and r.out.Get('mapcoef',typ='best')), None)
        if handdet and handdet.out.Get('mapcoef',typ='best'):
          p.inp.Add( handdet.out.Get('mapcoef',typ='best') )
          if hasattr(p,'inp2') and hasattr(handdet,'out2') and handdet.out2.Get('mapcoef',typ='best'):
            p.inp2.Add( handdet.out2.Get('mapcoef',typ='best') )
        else:
          ph=next((r for ir,r in enumerate(self.run.processes) if ir<ip and r.nick in ('phas','refatompick') and r.out.Get('mapcoef',typ='best')), None)
          if ph and ph.out.Get('mapcoef',typ='best'):
            p.inp.Add( ph.out.Get('mapcoef',typ='best') )
            if hasattr(p,'inp2') and hasattr(ph,'out2') and ph.out2.Get('mapcoef',typ='best'):
              p.inp2.Add( ph.out2.Get('mapcoef',typ='best') )
      # pass (almost) all output containers from previous steps to input of this process p
      if ip>0 and not prev_proc.IsNonFalseVirtPar('no_output_to_next_step'):
        if hasattr(prev_proc,'out2'):
          # pass all for hand 2 (if present)
          for o in prev_proc.out2.GetAll(stored_order=True):
            if o.typ!='dluz' and (o.typ!='weighted' or not [m for m in prev_proc.out2.GetAll('mapcoef') if m.typ!='weighted']):
              p.inp2.AddCopy(o,propagate=True)
        # pass all except the specified exceptions from the previous step
        for o in prev_proc.out.GetAll(stored_order=True):
          # we exclude any luzd output as that is generally not meant to be used by next steps
          # similarly, we use best rather than weighted mapcoef - any exceptions should be specified here
          if o.typ!='dluz' and (o.typ!='weighted' or not [m for m in prev_proc.out.GetAll('mapcoef') if m.typ!='weighted']):
            p.inp.AddCopy(o,propagate=True)
      # pass what was specified by the input
      if not emulate:
        p.inp.DelayedInit(ip,self)
        if hasattr(p,'inp2'): p.inp2.DelayedInit(ip,self)


  def RunPostprocess(self,restore=True,emulate=False,*args,**kwargs):
    if emulate:
      os.chdir(self.previous_cwd)
      try:
        os.remove(self.GetLogFileName())
        os.rmdir(self.rundir)
      except (OSError,):
        if hasattr(sys,'exc_clear'): sys.exc_clear()
    else:
      mapout=self.run.out.Get('mapcoef')
      if mapout:
        outmtz=self.run.out.Get('datafile',filetype='mtz',custom='refmac')
        if not outmtz:
          outmtz=self.run.out.Get('mapcoef',filetype='mtz',inp_cont=self.inp.Get('fsigf',xname=mapout.xname,dname=mapout.dname))
        if outmtz and hasattr(self,'hklout') and outmtz.GetFileName()!=self.hklout:
          if crvapi:
            outmtzstr = crvapi.SetMtzFWT(self,outmtz)
          else:
            outmtzstr = outmtz.GetFileName()
          shutil.copy(outmtzstr, self.hklout)
          # preparing the maps for jsCoFe - should be not needed anymore - remove code after making sure
          #if not self.ccp4i2:
            #outmap=self.run.out.Get('mapcoef',filetype='map',typ='weighted',inp_cont=self.inp.Get('fsigf',xname=mapout.xname,dname=mapout.dname))
            #if outmap: shutil.copy(outmap.GetFileName('map'), self.hklout+'.map')
            #outdmap=self.run.out.Get('mapcoef',filetype='map',typ='diff',conv_opts=['no_rewrite'],inp_cont=self.inp.Get('fsigf',xname=mapout.xname,dname=mapout.dname))
            #if outdmap: shutil.copy(outdmap.GetFileName('map'), self.hklout+'_diff.map')
      outpdb=self.run.out.Get(filetype='pdb',typ=('partial+substr','partial'))
      if not outpdb:  outpdb=self.run.out.Get(filetype='pdb')
      if outpdb and hasattr(self,'xyzout') and outpdb.GetFileName()!=self.xyzout:
        shutil.copy(outpdb.GetFileName(), self.xyzout)
      outpdb=self.run.out.Get(filetype='pdb',typ='substr')
      if outpdb and hasattr(self,'xyzsubout') and outpdb.GetFileName()!=self.xyzsubout:
        shutil.copy(outpdb.GetFileName('pdb'), self.xyzsubout)
      process.RunPostprocess(self,restore)
    self.Verdict()
    if hasattr(self,'logout'):
      sys.stdout = self.stdout_save

  def Verdict(self):
    """returns verdict in form of three statements as currently experimentally used for Cloud"""
    # nothing fancy, just a very simple ad-hoc approach is used as of now.
    mb=next((p for p in reversed(self.run.processes) if p.nick in ('comb_phdmmb','mbref')),None)
    if not mb or not mb.result_str or not hasattr(mb,'ph') or not hasattr(mb.ph,'R') or not mb.ph.R:
      return
    files=["verdict1","verdict2","verdict3"]
    files=[os.path.join(self.rundir,f) for f in files]
    score=min(1.,max(0.,(1.0-7.0*(mb.ph.R-0.2)**2)))
    with open(files[0],'w') as g:
      g.write( str(score) )
    with open(files[1],'w') as g:
      g.write(mb.result_str)
    with open(files[2],'w') as g:
      if score>0.9:
        g.write("The built model should be of a good quality.  It is suggested to inspect and further improve it in Coot.")
      elif score>0.8:
        g.write("It may be still possible to further improve the quality of the model by a follow-up Crank2 job or to proceed with manual corrections in Coot.")
      elif score>0.5:
        g.write("The model is incomplete but the solution looks promising. A follow-up Crank2 job, with a larger number of building cycles and possibly changing parameters such as solvent content, is suggested to further improve the model.")
      else:
        substrdet=next((p for p in self.run.processes if p.nick=='substrdet'),None)
        handdet=next((p for p in self.run.processes if p.nick=='handdet'),None)
        dmfull=next((p for p in self.run.processes if p.nick=='dmfull'),None)
        done=False
        if handdet and dmfull and substrdet and handdet.guess is not None and dmfull.guess is not None and substrdet.guess is not None:
          total=sum((handdet.guess,substrdet.guess,dmfull.guess))
          if total>=5:
            g.write("The anomalous substructure looks promising.  You can try to rerun Crank2 changing parameters important after substructure determination such as solvent content, model building parameters etc")
            done=True
          elif total<=1:
            g.write("Judging from scores from the structure solution process, the obtained substructure could be wrong. You can try to rerun Crank2 again, changing substructure determination parameters such as resolution cutoff, higher number of trials, the program used etc")
            done=True
        if not done:
          if score>0.33:
            g.write("It is unclear whether the built model is useful or not.  You can still try a follow-up Crank2 job, with a larger number of building cycles and possibly changing parameters such as solvent content.")
          else:
            g.write("The built model is not likely to be useful and it is not clear whether the substructure has been found.  You can try a new Crank2 job modifying the input parameters for substructure determination, solvent content, model building parameters etc.")
