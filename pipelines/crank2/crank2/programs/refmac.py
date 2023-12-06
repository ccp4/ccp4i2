#!/usr/bin/python
import os,sys,shutil,time
from program import program
import common

class refmac(program):
  name="REFMAC5"
  binary="refmac5"
  labelout_prefix="REFM_"
  ccp4_parsing=True
  stat={}
  stat['fom'] = common.stats(name='FOM', regexp=r"Overall figure of merit\s+=\s+(\S+)")
  stat['correl'] = common.stats(name='correlation coefficient', regexp=r"Overall correlation coefficient\s+=\s+(\S+)")
  stat['rfact'] = common.stats(name='R factor', regexp=r"Overall R factor\s+=\s+(\S+)", multiple=True)
  stat['rfree'] = common.stats(name='R-free factor', regexp=r"Free R factor\s+=\s+(\S+)", multiple=True)
  stat['version'] = common.stats(name='version', regexp=r"refmac5;\s+version\s+(\S+)")
  stat['lg_fom'] = common.stats(name='loggraph FOM info', 
                   regexp=r"\$TABLE: Cycle\s+\S+\s+Fom\(<cos\(DelPhi\)>")
  references = ( "Murshudov GN, Skubak P, Lebedev AA, Pannu NS, Steiner RA, Nicholls RA, " +
                 "Winn MD, Long F and Vagin AA (2011) REFMAC5 for the refinement of " +
                 "macromolecular crystal structures. Acta Cryst. D67, 355-367.", )

  def Init(self):
    self.outfilename = { 'pdb': self.name+'.pdb', \
                         'pdbnat': self.name+'_nat.pdb', \
                         'mtz': self.name+'.mtz' \
                       }
    self.basiclib_list = [ "SE","S","FE","ZN","HG","BR","CL","CU","I","ZR","PT","U","XE","AG","ARS","AL","AU","CA","CD","CO","CR","EU","F","K","LI","MG","MN","NA","NI","OS","PB","TA","W","UNX" ]

  def CheckBinary(self):
    program.CheckBinary(self)
    self.SetRunDir('checkbinary')
    self.runname='{0}_version'.format(self.name)
    self.ExternalRun( [self.binary, '-i'], [] )
    ver=self.GetStat('version')
    if ver in ('5.7.0028','5.7.0029','5.7.0030','5.7.0031'):
      common.Warning('Version {0} of {1} is known to contain bugs leading to significantly worse results.'.format(ver,self.name))
      common.Warning('It is STRONGLY RECOMMENDED TO STOP NOW and update {0} to 5.7.0032 or higher.'.format(self.name))
      time.sleep(3)
    if not ver.startswith('5.'):
      common.Warning('Version {0} of {1} might not be supported. It is recommended to use the current CCP4 version.'.format(ver,self.name))
      time.sleep(3)

  def Interact_output(self, line, popen, empty):
    #if not empty:
    self.process.UpdateInfo(line,popen)

  def TreatInput(self):
    exper=self.process.GetParam('target',capitalized=True)
    obj,N=self.inp.GetTargetObjects(exper,modtyp='pdb',fsftyp='mtz',accept_anom_nosubstr=self.process.IsNonFalseParam('no_substr'))
    if obj is None:
      common.Error('The input objects for target {0} of {1} could not be retrieved'.format(exper,self.name))
    # model(s) to be refined
    if not obj['mod']:
      common.Error('No model to be refined by {0}'.format(self.name))
    # SAD and SIRAS input
    if exper=='SIRAS':
      self.fn=obj['fn']
      self.SetKey( 'labin', ('FN='+obj['fn'].GetLabel('f'),'SIGFN='+obj['fn'].GetLabel('sigf')) )
    if exper=='SAD' or exper=='SIRAS': 
      self.SetArg( 'hklin', obj['f+'].GetFileName('mtz') )
      self.AddToKey('labin', ('F+='+obj['f+'].GetLabel('f'),'SIGF+='+obj['f+'].GetLabel('sigf')) )
      self.AddToKey('labin', ('F-='+obj['f-'].GetLabel('f'),'SIGF-='+obj['f-'].GetLabel('sigf')) )
      ind=0
      for at in obj['mod'].GetAtomTypes(getlist=True):
        fp,fpp,dn,att=obj['mod'].Getfpfpp(at,obj['f+'].GetDataName())
        if fp is not None and fpp is not None:
          self.SetKey( 'ANOM FORM', (at,fp,fpp), insert_index=ind )  # inserting to the beginning of the file - workaround for refmac failing if there was a 0.0,0.0 before
          ind+=1
    # MLHL and RICE input
    elif exper=='MLHL':
      self.SetArg( 'hklin', obj['f'].GetFileName('mtz') )
      self.SetKey(  'labin', ('FP='+obj['f'].GetLabel('f'),'SIGFP='+obj['f'].GetLabel('sigf')) )
      if obj['hl']:
        self.AddToKey('labin', ('HLA='+obj['hl'].GetLabel('hla'),'HLB='+obj['hl'].GetLabel('hlb')))
        self.AddToKey('labin', ('HLC='+obj['hl'].GetLabel('hlc'),'HLD='+obj['hl'].GetLabel('hld')))
      else:
        self.AddToKey('labin', ('PHIB='+obj['phfom'].GetLabel('ph'),'FOM='+obj['phfom'].GetLabel('fom')))
    elif exper=='RICE':
      self.SetArg( 'hklin', obj['f'].GetFileName('mtz') )
      self.SetKey( 'labin', ('FP='+obj['f'].GetLabel('f'),'SIGFP='+obj['f'].GetLabel('sigf')) )
    # the map we want to combine with (at the moment, phase comb. key is DM in refmac but any map can be used)
    if self.GetKey('DM') or self.GetKey('DM2'):
      comb_typ=self.process.GetVirtPar('comb_with_typ')
      if comb_typ!='no_phcomb':
        self.SetKey('scpa','1')
        ph = self.inp.Get('mapcoef',typ=comb_typ,col='ph',filetype='mtz')
        if ph:
          self.AddToKey('labin', ('FPART1='+ph.GetLabel('f'),'PHIP1='+ph.GetLabel('ph')))
        else:
          common.Error("No phase information to combine with for {0}".format(self.name))
    # just saving so that output can reach the objects easily
    self.mdl,self.mdlnat,self.f = obj['mod'],None,None
    for f in [f  for f in ('f+','fn','f') if f in obj]:
      self.f=obj[f]
    # xyzin
    if 'modn' in obj and exper=='SIRAS':
      self.mdlnat=obj['modn']
      self.SetArg('xyzin',obj['modn'].GetFileName('pdb'))
      self.SetArg('xyzin2',obj['mod'].GetFileName('pdb'))
    else:
      self.SetArg('xyzin',obj['mod'].GetFileName('pdb'))
    # free assignment
    free=self.inp.Get('exclude',typ='freeR',col='free')
    if free:
      self.AddToKey( 'labin', "FREE={0}".format(free.GetLabel('free')) )
    # read in luzzatti D parameters
    if self.inp.Get('datafile',typ='dluz') and self.rundir:
      if not self.inp.Get('datafile',typ='dluz').GetFileName().startswith(self.rundir):
        shutil.copy(self.inp.Get('datafile',typ='dluz').GetFileName(), self.rundir)
      self.SetKey('REFI','RLUZZ')
    self.obj=obj

  def TreatParams(self):
    if self.process.IsNonFalseVirtPar('cycles') and not self.IsKey('ncyc'):
      self.SetKey('ncyc', self.process.GetVirtPar('cycles'))
    if self.process.IsNonFalseVirtPar('beta') and not self.IsKey('dmbl'):
      self.SetKey('dmbl', self.process.GetVirtPar('beta'))
    if self.process.nick=='phcomb' or (self.process.nick=='ref' and self.process.GetVirtPar('phcomb')):
      if not self.GetKey('DM2') and not self.GetKey('DM'):
        if self.process.GetVirtPar('target')=='SIRAS' or self.process.nick=='ref':
          self.SetKey('DM2')
        else:
          numat=0
          try:
            import gemmi
            struct=gemmi.read_structure(self.inp.Get('model',typ='substr').GetFileName())
            for s in struct[0].all():
              numat+=1
          except ImportError:
            pass
          if numat>3:
            self.SetKey('DM2')
          else:
            self.SetKey('DM')
      if not self.IsKey('ncyc'):
        self.SetKey('ncyc',0)
    if self.process.nick!='ref':
      self.SetKey('solv','no')
    if (self.process.nick!='ref' or self.process.GetVirtPar('phcomb')): 
      if not any('SUBS' in v for v in self.GetKey('refi',allval=True,capitalized=True)):
        self.SetKey('refi',('subs','yes'))
      if self.process.nick=='phcomb' and not any('TYPE' in v for v in self.GetKey('refi',allval=True,capitalized=True)):
        self.SetKey('refi',('type','unre'))
    if self.process.nick in ('ref','phas'):
      at = self.inp.Get(has_atomtypes=True)
      if self.process.nick=='phas' and at and at.GetAtomType()=='S': #impose auto disuph. links in phasing
        self.SetKey('make',('link','y'))
        ###self.SetKey('dist',3)
      if self.GetParam('allow_basiclib') is not False and at and at.GetAtomType() in self.basiclib_list:
        dummy_env = os.path.join( os.getenv('CRANK2'), 'bin', 'dummy' ) if os.getenv('CRANK2') else None
        basiclib = None
        if not dummy_env or not os.path.isdir(dummy_env):
          dummy_env = os.path.join( os.getenv('CRANK'), 'bin', 'dummy' ) if os.getenv('CRANK') else None
        if os.path.isdir(dummy_env) and os.path.isdir(os.path.join(dummy_env,'basiclib')):
          basiclib = os.path.join(dummy_env,'basiclib')
        if not basiclib and os.getenv('CCP4') and os.path.isdir( os.path.join( os.getenv('CCP4'), 'share', 'ccp4i', 'crank', 'bin', 'dummy', 'basiclib' ) ):
          basiclib = os.path.join( os.getenv('CCP4'), 'share', 'ccp4i', 'crank', 'bin', 'dummy', 'basiclib' )
        if not basiclib and os.getenv('CCP4') and os.path.isdir( os.path.join( os.getenv('CCP4'), 'share', 'ccp4i2', 'pipelines', 'crank2', 'crank2', 'basiclib' ) ):
          basiclib = os.path.join( os.getenv('CCP4'), 'share', 'ccp4i2', 'pipelines', 'crank2', 'crank2', 'basiclib' )
        if basiclib:
          self.SetArg('libin', os.path.join(basiclib,'mon_lib_list_basic.cif'))
          self.SetKey('make', ('newl','noex'))
          self.env['CLIBD_MON'] = os.path.join(basiclib,'')
    self.SetKey('REFI','WLUZZ')
    program.TreatParams(self)

  def DefineOutput(self):
    besttyp='combined'
    if self.process.nick=='phas':
      besttyp='best'
    self.outmapcc = self.out.AddNew( 'mapcoef', self.outfilename['mtz'], typ=besttyp, xname=self.f.xname, dname=self.f.dname)
    if self.process.nick=='ref':
      self.outmapcc.SetLabel( ['ph','fom','hla','hlb','hlc','hld'] )
    else:
      self.outmapcc.SetLabel( ['f','ph','fom','hla','hlb','hlc','hld'] )
    self.outmapcw = self.out.AddNew( 'mapcoef', self.outfilename['mtz'], typ='weighted', xname=self.f.xname, dname=self.f.dname)
    self.outmapcw.SetLabel( ['f','ph'] )
    self.outmapdiff = self.out.AddNew( 'mapcoef', self.outfilename['mtz'], typ='diff', xname=self.f.xname, dname=self.f.dname)
    self.outmapdiff.SetLabel( ['f','ph'] )
    self.outmapcad, self.outmapcan = None, None
    if self.process.GetVirtPar('target') in ('SAD','SIRAS') and self.obj['mod'].GetType() in ('partial+substr','substr'):
      self.outmapcad = self.out.AddNew( 'mapcoef', self.outfilename['mtz'], typ='anom-diff', xname=self.f.xname, dname=self.f.dname )
      self.outmapcad.SetLabel( ['f','ph'] )
    elif [v for v in self.GetKey('anom',as_list=True) if v and str(v).lower().startswith('mapo')] and self.obj['mod'].GetType()=='partial':
      # diff anom diff equal to anom diff
      self.outmapcan = self.out.AddNew( 'mapcoef', self.outfilename['mtz'], typ='anom-diff', xname=self.f.xname, dname=self.f.dname )
      self.outmapcan.SetLabel( ['f','ph'] )
    self.outfsigf = self.out.AddNew( 'fsigf', self.outfilename['mtz'], typ='average-derived', xname=self.f.xname, dname=self.f.dname, custom=self.f.custom)
    self.outfsigf.SetLabel( ['f','sigf'] )
    self.outmodel = self.out.AddCopy(self.mdl)
    self.out.SetFileToChild(self.outmodel, self.outfilename['pdb'], 'pdb')
    #self.outmodel = self.out.AddNew( 'model', self.outfilename['pdb'], typ=self.mdl.GetType(), atomtypes=self.mdl.GetAtomTypes() )
    if self.mdlnat:
      self.outmodelnat = self.out.AddCopy(self.mdlnat)
      self.out.SetFileToChild(self.outmodelnat, self.outfilename['pdbnat'], 'pdb')
      #self.outmodelnat = self.out.AddNew( 'model', self.outfilename['pdbnat'], typ=self.mdlnat.GetType(), atomtypes=self.mdlnat.GetAtomTypes() )
    if any(val for val in self.GetKey('REFI',allval=True,capitalized=True) if 'WLUZ' in val):
      self.out.AddNew( 'datafile', typ='dluz', filename='luzzd' )
    self.out.AddNew( 'datafile', filename=self.outfilename['mtz'], custom='refmac', filetype='mtz' )

  def TreatOutput(self):
    self.SetKey('phout')
    if self.process.nick=='ref':
      self.SetKey( 'labout', 'PHCOMB='+self.outmapcc.GetLabel('ph') )
    elif self.process.nick=='phcomb':
      self.SetKey( 'labout', ('FB='+self.outmapcc.GetLabel('f'),'PHCOMB='+self.outmapcc.GetLabel('ph')) )
    else:
      self.SetKey( 'labout', ('FB='+self.outmapcc.GetLabel('f'),'PHIB='+self.outmapcc.GetLabel('ph')) )
    self.AddToKey( 'labout', 'FOM='+self.outmapcc.GetLabel('fom') )
    if self.process.nick=='phas':
      self.AddToKey( 'labout', ('HLA='+self.outmapcc.GetLabel('hla'),'HLB='+self.outmapcc.GetLabel('hlb')) )
      self.AddToKey( 'labout', ('HLC='+self.outmapcc.GetLabel('hlc'),'HLD='+self.outmapcc.GetLabel('hld')) )
    else:
      self.AddToKey( 'labout', ('HLACOMB='+self.outmapcc.GetLabel('hla'),'HLBCOMB='+self.outmapcc.GetLabel('hlb')) )
      self.AddToKey( 'labout', ('HLCCOMB='+self.outmapcc.GetLabel('hlc'),'HLDCOMB='+self.outmapcc.GetLabel('hld')) )
    self.AddToKey( 'labout', ('FWT='+self.outmapcw.GetLabel('f'),'PHWT='+self.outmapcw.GetLabel('ph')) )
    self.AddToKey( 'labout', ('DELFWT='+self.outmapdiff.GetLabel('f'),'PHDELWT='+self.outmapdiff.GetLabel('ph')) )
    if self.outmapcad:
      self.AddToKey( 'labout', ('DELFAN='+self.outmapcad.GetLabel('f'),'PHDELAN='+self.outmapcad.GetLabel('ph')) )
    if self.outmapcan:
      self.AddToKey( 'labout', ('FAN='+self.outmapcan.GetLabel('f'),'PHAN='+self.outmapcan.GetLabel('ph')) )
    outf = ('FP','SIGFP')
    if self.f and ('FN='+self.f.GetLabel('f') in self.GetKey('labin') or 'FP='+self.f.GetLabel('f') in self.GetKey('labin')):
      outf = (self.f.GetLabel('f'),self.f.GetLabel('sigf'))
    self.AddToKey( 'labout', ( outf[0]+'='+self.outfsigf.GetLabel('f'), outf[1]+'='+self.outfsigf.GetLabel('sigf')) )
    self.SetArg( 'hklout', self.outmapcc.GetFileName('mtz') )
    if self.mdlnat:
      self.SetArg( 'xyzout', self.outmodelnat.GetFileName('pdb') )
      self.SetArg( 'xyzout2', self.outmodel.GetFileName('pdb') )
    else:
      self.SetArg( 'xyzout', self.outmodel.GetFileName('pdb') )
