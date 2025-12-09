#!/usr/bin/python
import os,sys
from program import program
import common

class multicomb(program):
  name="Multicomb"
  binary="multicomb"
  labelout_prefix="MULT_"
  ccp4_parsing=True
  stat={}
  stat['fom'] = common.stats(name='FOM', regexp=r"Overall MEAN FOM is\s+(\S+)",multiple=True)
  stat['correl'] = common.stats(name='correlation coefficient', regexp=r"Overall Correlation is\s+(\S+)")
  references = ( "Skubak P, Waterreus WJ and Pannu NS (2010) " +
                 "Multivariate phase combination improves automated crystallographic" +
                 "model building.  Acta Cryst. D66, 783-788.", )


  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    exper=self.GetKey('targ')
    if exper=='PSAD':
      exper='SAD'
    obj,N=self.inp.GetTargetObjects(exper)
    if obj is None:
      common.Error('The input objects for target {0} of {1} could not be retrieved'.format(exper,self.name))
    # SAD input
    if exper=='SAD': 
      # crystal and dataset name are compulsory for Multicomb...
      self.SetKey('XTAL', str(obj['f+'].xname))
      self.SetKey('DNAME', str(obj['f+'].dname))
      self.SetArg( 'hklin', obj['f+'].GetFileName('mtz') )
      self.SetKey(  'colu', ('F+='+obj['f+'].GetLabel('f'), 'SF+='+obj['f+'].GetLabel('sigf')) )
      self.AddToKey('colu', ('F-='+obj['f-'].GetLabel('f'), 'SF-='+obj['f-'].GetLabel('sigf')) )
      self.SetKey('modl',obj['mod'].GetFileName('pdb'))
      for at in obj['mod'].GetAtomTypes(getlist=True):
        fp,fpp,dn,att=obj['mod'].Getfpfpp(at,obj['f+'].dname)
        if fp is not None and fpp is not None:
          self.SetKey('FORM', '{0} FP={1} FPP={2}'.format(at,fp,fpp))
      self.f=obj['f+']
    # MLHL input
    elif exper=='MLHL':
      # crystal and dataset name are compulsory for Multicomb...
      self.SetKey('XTAL', str(obj['f'].xname))
      self.SetKey('DNAME', str(obj['f'].dname))
      self.SetArg( 'hklin', obj['f'].GetFileName('mtz') )
      self.SetKey(  'colu', ('F='+obj['f'].GetLabel('f'),'SF='+obj['f'].GetLabel('sigf')) )
      self.AddToKey('colu', ('HLA='+obj['hl'].GetLabel('hla'),'HLB='+obj['hl'].GetLabel('hlb')))
      self.AddToKey('colu', ('HLC='+obj['hl'].GetLabel('hlc'),'HLD='+obj['hl'].GetLabel('hld')))
      self.f=obj['f']
    else:
      common.Error("Wrong target {0} specified for {1}".format(self.GetKey('targ'),self.name))
    # the calcuted map to be combined with
    ph = self.inp.Get('mapcoef',typ='densmod',filetype='mtz',col='ph')
    if ph:
      self.AddToKey('colu', ('FC1='+ph.GetLabel('f'),'PC1='+ph.GetLabel('ph')))
    else:
      common.Error("No phase information to combine with for {0}".format(self.name))

  def TreatParams(self):
    if self.process.GetVirtPar('cycles'):
      self.SetKey('cycl', self.process.GetVirtPar('cycles'))
    if self.process.GetVirtPar('beta'):
      self.SetKey('beta', self.process.GetVirtPar('beta'))
    if not self.process.GetVirtPar('target') and not self.GetKey('targ'):
      common.Error('Target/experiment has to be specified for program {0}'.format(self.name))
    # just direct translation without any checks here
    if not self.GetKey('targ'):
      if self.process.GetVirtPar('target')=='SAD':
        self.SetKey('targ','PSAD')
      else:
        self.SetKey('targ','MLHL')
    else:
      self.CapitalizeKey('targ')
    program.TreatParams(self)

  def DefineOutput(self):
    outmtz=self.outfilename['mtz']
    #outpdb=self.nick+'.pdb'
    outpdb='heavy-1.pdb'
    substr=self.inp.Get('model',typ='substr')
    # perhaps the crystal name should be obtained from input f rather than model?
    if substr:
      xname=substr.xname
    if substr and self.GetKey('targ')!='MLHL':
      self.outmodel = self.out.AddNew( 'model', outpdb, typ='substr', atomtypes=substr.GetAtomTypes(), \
                      xname=xname )
    self.outmapcc = self.out.AddNew( 'mapcoef', outmtz, typ='combined', xname=xname )
    self.outmapcc.SetLabel( ['f','ph','fom','hla','hlb','hlc','hld'] )
    self.outmapcw = self.out.AddNew( 'mapcoef', outmtz, typ='weighted', xname=xname )
    self.outmapcw.SetLabel( ['f','ph'] )
    self.outfsigf = self.out.AddNew( 'fsigf', self.outfilename['mtz'], typ='average-derived', xname=self.f.xname, dname=self.f.dname, custom=self.f.custom)
    self.outfsigf.SetLabel( ['f','sigf'] )

  def TreatOutput(self):
    self.SetKey( 'labout', ('FB='+self.outmapcc.GetLabel('f'),'PHIB='+self.outmapcc.GetLabel('ph')) )
    self.AddToKey( 'labout', 'FOM='+self.outmapcc.GetLabel('fom') )
    self.AddToKey( 'labout', ('HLA='+self.outmapcc.GetLabel('hla'),'HLB='+self.outmapcc.GetLabel('hlb')) )
    self.AddToKey( 'labout', ('HLC='+self.outmapcc.GetLabel('hlc'),'HLD='+self.outmapcc.GetLabel('hld')) )
    self.AddToKey( 'labout', ('FWT='+self.outmapcw.GetLabel('f'),'PHWT='+self.outmapcw.GetLabel('ph')) )
    self.AddToKey( 'labout', ( 'F='+self.outfsigf.GetLabel('f'), 'SIGF='+self.outfsigf.GetLabel('sigf')) )
    self.SetArg( 'hklout', self.outmapcc.GetFileName('mtz') )
    #self.SetArg( 'xyzout', self.outmodel.GetFileName('pdb') )
