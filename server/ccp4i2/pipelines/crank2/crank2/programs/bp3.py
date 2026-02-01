#!/usr/bin/python
import os,sys
from ..program import program
from .. import common

class bp3(program):
  name="Bp3"
  binary="bp3"
  labelout_prefix="BP3_"
  ccp4_parsing=True
  stat={}
  stat['fom'] = common.stats(name='FOM', regexp=r"The overall FOM is\s+(\S+)",multiple=True)
  stat['luzzati'] = common.stats(name='average Luzzati parameter', regexp=r"the average anomalous Luzzati error is\s+(\S+)")
  stat['version'] = common.stats(name='version', regexp=r"bp3;\s+version\s+(\S+)")
  stat['ccp4_version'] = common.stats(name='version', regexp=r"CP4 software suite: library version\s+(\S+)")
  references = ( "Pannu NS, Waterreus WJ,  Skubak P, Sikharulidze I, Abrahams JP and "+
                 "de Graaff RAG (2011) Recent advances in the CRANK software suite "+
                 "for experimental phasing. Acta Cryst. D67, 331-337.", )

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    exper=self.GetKey('targ')  if self.GetKey('targ') not in ('MSRS','UNCO')  else 'SIRAS'
    obj,N=self.inp.GetTargetObjects(exper)
    if obj is None:
      common.Error('The input objects for target {0} of {1} could not be retrieved'.format(exper,self.name))
    # SAD input
    if exper in ('SAD','SIRAS','MAD'):
      # crystal and dataset name are compulsory for BP3...
      if exper=='SIRAS' and 'fn' in obj and obj['fn']:
        self.SetKey('XTAL', str(obj['fn'].xname))
        self.SetKey('DNAME', str(obj['fn'].dname))
        self.SetKey(  'colu', ('F='+obj['fn'].GetLabel('f'), 'SF='+obj['fn'].GetLabel('sigf')) )
      self.SetKey('XTAL', str(obj['f+'].xname))
      for i in range(N):
        istr=''  if i==0  else str(i+1)
        fpl, fmi = 'f{}+'.format(istr), 'f{}-'.format(istr)
        self.SetKey('DNAME', str(obj[fpl].dname))
        self.SetArg( 'hklin', obj[fpl].GetFileName('mtz') )
        self.SetKey(  'colu', ('F+='+obj[fpl].GetLabel('f'), 'SF+='+obj[fpl].GetLabel('sigf')) )
        self.AddToKey('colu', ('F-='+obj[fmi].GetLabel('f'), 'SF-='+obj[fmi].GetLabel('sigf')) )
        for at in obj['mod'].GetAtomTypes(getlist=True):
          fp,fpp,dn,att=obj['mod'].Getfpfpp(at,obj[fpl].dname)
          if fp is not None and fpp is not None:
            self.SetKey('FORM', '{0} FP={1} FPP={2}'.format(at,fp,fpp))
      self.SetKey('modl',obj['mod'].GetFileName('pdb'))
    if exper=='MAD':
      self.SetKey('PHASE')

  def TreatParams(self):
    if self.process.IsNonFalseVirtPar('cycles'):
      self.SetKey('cycl', self.process.GetVirtPar('cycles'))
    if self.process.IsNonFalseVirtPar('beta'):
      self.SetKey('beta', self.process.GetVirtPar('beta'))
    if not self.process.GetVirtPar('target') and not self.GetKey('targ'):
      common.Error('Target/experiment has to be specified for program {0}'.format(self.name))
    if not self.GetKey('targ'):
      if self.process.GetVirtPar('target') in ('SAD','MAD'):
        self.SetKey('targ', self.process.GetVirtPar('target'))
      elif self.process.GetVirtPar('target')=='SIRAS':
        #self.SetKey('targ', 'MSRS')
        self.SetKey('targ', 'UNCO')
    else:
      self.CapitalizeKey('targ')
    program.TreatParams(self)

  def DefineOutput(self):
    # most of these seem to be unfortunately hardcoded (?) - no control over it from here.
    outmtz, outmtz_oh = self.outfilename['mtz'], self.outfilename['mtz']+'-oh'
    outi = '2' if self.GetKey('targ') in ('MSRS','UNCO') else '1'
    outpdb, outpdb_oh  = 'heavy-{0}.pdb'.format(outi), 'heavy-oh-{0}.pdb'.format(outi)
    atomtypes,xname={},None
    if self.inp.Get('model'):
      atomtypes=self.inp.Get('model').GetAtomTypes()
      xname=self.inp.Get('model').GetCrystalName()
    if not self.GetKey('noha'):
      self.out.AddNew('model', outpdb_oh, filetype='pdb', typ='substr', custom='otherhand', \
                      atomtypes=atomtypes, xname=xname)
      self.outmapcc_oh = self.out.AddNew( 'mapcoef', outmtz_oh, typ='best', xname=xname, custom='otherhand' )
      self.outmapcc_oh.SetLabel( ['f','ph','fom','hla','hlb','hlc','hld'] )
    self.out.AddNew('model', outpdb, filetype='pdb', typ='substr', \
                    atomtypes=atomtypes, xname=xname)
    self.outmapcc = self.out.AddNew( 'mapcoef', outmtz, typ='best', xname=xname )
    self.outmapcc.SetLabel( ['f','ph','fom','hla','hlb','hlc','hld'] )
    self.outmapcad = None
    if self.process.GetVirtPar('target') in ('SAD','SIRAS') and self.GetKey('targ')!='UNCO':
      self.outmapcad = self.out.AddNew( 'mapcoef', self.outfilename['mtz'], typ='anom-diff', xname=xname )
      self.outmapcad.SetLabel( ['f','ph'] )


  def TreatOutput(self):
    self.SetKey( 'labout', ('FB='+self.outmapcc.GetLabel('f'),'PHIB='+self.outmapcc.GetLabel('ph')) )
    self.AddToKey( 'labout', 'FOM='+self.outmapcc.GetLabel('fom') )
    self.AddToKey( 'labout', ('HLA='+self.outmapcc.GetLabel('hla'),'HLB='+self.outmapcc.GetLabel('hlb')) )
    self.AddToKey( 'labout', ('HLC='+self.outmapcc.GetLabel('hlc'),'HLD='+self.outmapcc.GetLabel('hld')) )
    if self.outmapcad:
      self.AddToKey( 'labout', ('FDIFF='+self.outmapcad.GetLabel('f'),'PDIFF='+self.outmapcad.GetLabel('ph')) )
    self.SetArg( 'hklout', self.outmapcc.GetFileName('mtz') )
    #self.SetArg( 'xyzout', self.outmodel.GetFileName('pdb') )
