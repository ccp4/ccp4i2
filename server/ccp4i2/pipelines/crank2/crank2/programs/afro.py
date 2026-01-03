#!/usr/bin/python
import os,sys
from program import program
import common

class afro(program):
  name="Afro"
  binary="afro"
  labelout_prefix="AFRO_"
  ccp4_parsing=True
  stat={}
  stat['version'] = common.stats(name='version', regexp=r"afro;\s+version\s+(\S+)")
  stat['ccp4_version'] = common.stats(name='version', regexp=r"CCP4 software suite: library version\s+(\S+)")
  stat['resol_bins'] = common.stats(name='resolution in bin', 
     regexp=r'Bin  LoRes  HiRes    Res   STOL2    Refls     F/SigmaF   AnoRefls    Dano/SigDano \$\$\s\$\$'+r'(?:\s+\d+\s+\S+\s+\S+\s+(\d+\.\d+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s)?'*99+r'\$\$')
  stat['dano_bins'] = common.stats(name='anomalous difference in bin', 
#     regexp=r'Bin  LoRes  HiRes    Res   STOL2    Refls     F/SigmaF   AnoRefls    Dano/SigDano \$\$\s\$\$'+r'(?:\s+\d+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+\.\d+)\s)?'*99+r'\$\$')
     regexp=r'IsoRefls   Diso/SigDiso   AnoRefls    Dano/SigDano \$\$\s\$\$'+r'(?:\s+\d+\s+\S+\s+\S+\s+\s+\S+\s+\S+\s+\S+\s+\S+\s+(\d+\.\d+)\s)?'*99+r'\$\$')
  references = ( "Pannu NS, Waterreus WJ,  Skubak P, Sikharulidze I, Abrahams JP and "+
                 "de Graaff RAG (2011) Recent advances in the CRANK software suite "+
                 "for experimental phasing. Acta Cryst. D67, 331-337.", )

  def Init(self):
    self.outfilename = { 'mtz': self.name+'.mtz' }

  def TreatInput(self):
    exper=self.GetKey('targ')
    self.obj,N=self.inp.GetTargetObjects(exper)
    if self.obj is None:
      common.Error('The input objects for target {0} of {1} could not be retrieved'.format(exper,self.name))
    # SAD input
    if exper=='SAD': 
      # crystal and dataset name are compulsory for Afro...
      self.SetKey('XTAL', self.obj['f+'].xname)
      if self.obj['mod'].GetAtomType():
        self.SetKey('ATOM', self.obj['mod'].GetAtomType())
      else:
        common.Error('Substructure atom type required by {}'.format(self.name))
      if self.obj['mod'].exp_num_atoms:
        self.SetKey('NUMB', self.obj['mod'].exp_num_atoms)
      else:  # this should normally not happen but if it does from some reason then starting from 10 and refining should be fine
        self.SetKey('NUMB', 10)
      if self.process.IsNonFalseVirtPar('bfactor') and not self.GetKey('BISO'):
        self.SetKey( 'BISO', self.process.GetVirtPar('bfactor') )
      else:
        self.SetKey('BISO', 25)
      # disabling no hires cutoff as it is done elsewhere
      self.SetKey('HIRE',0)
      self.SetKey('DNAME', self.obj['f+'].dname)
      if self.process.GetVirtPar('high_res_cutoff'):
        self.SetKey('RESO',(999.,self.process.GetVirtPar('high_res_cutoff')))
      self.SetArg( 'hklin', self.obj['f+'].GetFileName('mtz') )
      self.SetKey(  'colu', ('F+='+self.obj['f+'].GetLabel('f'), 'SF+='+self.obj['f+'].GetLabel('sigf')) )
      self.AddToKey('colu', ('F-='+self.obj['f-'].GetLabel('f'), 'SF-='+self.obj['f-'].GetLabel('sigf')) )
      for at in self.obj['mod'].GetAtomTypes(getlist=True):
        fp,fpp,dn,att=self.obj['mod'].Getfpfpp(at,self.obj['f+'].dname)
        if fp is not None and fpp is not None:
          self.SetKey('FORM', '{0} FP={1} FPP={2}'.format(at,fp,fpp))
      # if EXCL was inputted before then move here (so that EXCL can be passed by command line - works for SAD only)
      excl=['EOUT',]
      if self.GetKey('EXCL'):
        excl.extend(self.GetKey('EXCL',allval=True))
        self.UnsetParam('EXCL',is_arg=False,is_key=True)
      self.SetKey('EXCL',excl)
      #self.SetKey('EXCL','EOUT')
      self.SetKey('FOUT')

  def TreatParams(self):
    if not self.GetKey('targ'):
      if self.process.GetVirtPar('target')=='SAD':
        self.SetKey('targ','SAD')
    else:
      self.CapitalizeKey('targ')
    program.TreatParams(self)

  def DefineOutput(self):
    xname=self.obj['f+'].GetCrystalName()
    dname=self.obj['f+'].GetDataName()
    outmtz=self.outfilename['mtz']
    self.outfa = self.out.AddNew( 'fsigf', outmtz, typ='fa', xname=xname, dname=dname )
    self.outfa.SetLabel( ['f','sigf','e','sige'] )
    self.outfp = self.out.AddNew( 'fsigf', outmtz, typ='plus', xname=xname, dname=dname )
    self.outfp.SetLabel( ['f','sigf'], ignore_prefix=True )  # just using the default afro labels matching those in crank2
    self.outfm = self.out.AddNew( 'fsigf', outmtz, typ='minus', xname=xname, dname=dname )
    self.outfm.SetLabel( ['f','sigf'], ignore_prefix=True )
    self.outf = self.out.AddNew( 'fsigf', outmtz, typ='average', xname=xname, dname=dname)
    self.outf.SetLabel( ['f','sigf'], ignore_prefix=True )


  def TreatOutput(self):
    self.SetKey( 'labout',   ('FA='+self.outfa.GetLabel('f'), 'SGFA='+self.outfa.GetLabel('sigf')) )
    self.AddToKey( 'labout', ('EA='+self.outfa.GetLabel('e'), 'SGEA='+self.outfa.GetLabel('sige')) )
    self.SetArg( 'hklout', self.outfa.GetFileName('mtz') )
    #self.SetArg( 'hklout', self.outea.GetFileName('mtz') )
    #self.SetArg( 'xyzout', self.outmodel.GetFileName('pdb') )
