#!/usr/bin/python
import os,sys
from program import program
import common

class arpwarp_auto(program):
  name="ARP/wARP auto tracing/preparation"
  binary="auto_tracing.sh"
  #labelout_prefix="AARP_"
  references = ( "Perrakis A, Morris RM and Lamzin VS (1999) Automated protein "+
                 "model building combined with iterative structure refinement. "+
                 "Nature Struct. Biol. 6, 458-463)", )

  def TreatInput(self):
    # F,sigF
    self.fsf = self.inp.Get('fsigf',typ='average',col='f')
    if self.fsf:
      self.SetArg('datafile', self.fsf.GetFileName('mtz'))
      self.SetArg('fp', self.fsf.GetLabel('f'))
      self.SetArg('sigfp', self.fsf.GetLabel('sigf'))
    else:
      common.Error('No data inputted to {0}'.format(self.name))
    # PHIB+FOM
    mapc=self.inp.Get('mapcoef',typ=('best','combined'),col=('ph','fom'),filetype='mtz',inp_cont=self.fsf)
    if mapc:
      self.SetArg('phibest', mapc.GetLabel('ph'))
      self.SetArg('fom', mapc.GetLabel('fom'))
    #else:
    #  common.Error('Best phase and/or FOM not inputted to {0}'.format(self.name))
    # FREE
    excl=self.inp.Get('exclude',typ='freeR')
    if excl:
      self.SetArg('freelabin', excl.GetLabel('free'))
    # sequence and information typically derived from it (or inputted)
    seq=self.inp.Get('sequence')
    if seq:
      self.SetArg('seqin', seq.GetFileName())
    monom_obj=self.inp.Get(has_monomers_asym=True)
    if monom_obj and not self.GetArg('cgr'):
      self.SetArg('cgr', monom_obj.monomers_asym)
    if monom_obj and hasattr(monom_obj,'residues_mon') and not self.GetArg('residues'):
      self.SetArg('residues', monom_obj.residues_mon*monom_obj.monomers_asym)
    # F+/F- in case of SAD - just checks they exist
    if self.process.GetParam('target') == 'SAD':
      objs,N = self.inp.GetTargetObjects(self.process.GetParam('target'), modtyp='pdb', fsftyp='mtz')
      if objs is None:
        common.Error('Input data/model for ARP/wARP SAD could not be retrieved.')
      objs['f+'].GetFileName('mtz'), objs['f-'].GetFileName('mtz')
      # SAD as implemented in arpwarp.
      if self.process.GetParam('arp_sad'):
        self.SetArg('phaselabin', ' F+={0} SIGF+={1} F-={2} SIGF-={3} '.format( objs['f+'].GetLabel('f'), 
          objs['f+'].GetLabel('sigf'), objs['f-'].GetLabel('f'), objs['f-'].GetLabel('sigf')))
        self.SetArg('sad', '1')
        self.SetArg('heavyin', objs['mod'].GetFileName())
        f=objs['mod'].GetAtomTypes()[objs['mod'].GetAtomType()]
        if f:
          self.SetArg('wavelength', ' ANOM FORM {0}  {1} {2} '.format(objs['mod'].GetAtomType(),f[0][0],f[0][1]))

  def TreatParams(self):
    if self.process.IsNonFalseVirtPar('bigcyc') and not self.IsArg('buildingcycles'):
      self.SetArg('buildingcycles', self.process.GetVirtPar('bigcyc'))
    program.TreatParams(self)

  def DefineOutput(self):
    # arpwarp uses fixed output files...
    parfname='arp_warp_tracing.par'
    pdbfname='arp_warp_tracing.pdb'
    mtzfname='arp_warp_tracing.mtz'
    if not self.IsArg('jobId'):
      self.SetArg('jobId', '.')
    arp_dir=self.rundir
    if not self.GetArg('jobId').startswith('.'):
      arp_dir=os.path.join(self.rundir, self.GetArg('jobId'))
    if self.GetArg('parfile'):
      parfname=self.GetArg('parfile')
    else:
      out_model = self.out.AddNew( 'model', os.path.join(arp_dir,pdbfname), typ='partial', xname=self.fsf.xname )
      out_map = self.out.AddNew( 'mapcoef', os.path.join(arp_dir,mtzfname), typ='weighted', xname=self.fsf.xname )
      out_map.SetLabel( ['f','ph','fom'], ['FWT','PHWT','FOM'] )
    par_file = self.out.AddNew( 'datafile', os.path.join(arp_dir,parfname), typ='par' )

