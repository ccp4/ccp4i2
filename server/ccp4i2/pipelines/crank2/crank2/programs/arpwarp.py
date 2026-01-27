#!/usr/bin/python
import os,sys,re
from ..program import program
from .. import common

class arpwarp(program):
  name="ARP/wARP"
  binary="warp_tracing.sh"
  #labelout_prefix="ARP_"
  stat={}
  stat['res_built'] = common.stats(name='number of residues built', regexp=r"(\d+)\s+peptides, \d+\s+chains.")
  stat['frag_built'] = common.stats(name='number of fragments built', regexp=r"\d+\s+peptides, (\d+)\s+chains.")
  stat['build_cycle'] = common.stats(name='actual building cycle', regexp=r"Building Cycle (\d+)")
  stat['round_cycle'] = common.stats(name='the "best" round building cycle', regexp=r"Taking the results from Round (\d+)")
  stat['rfact'] = common.stats(name='R factor', regexp=r"After refmac, R =\s+(\S+)\s+\(Rfree =\s+\S+")
  stat['rfree'] = common.stats(name='R-free factor', regexp=r"After refmac, R =\s+\S+\s+\(Rfree =\s+(\d+\.\d+)")
  references = ( "Perrakis A, Morris RM and Lamzin VS (1999) Automated protein "+
                 "model building combined with iterative structure refinement. "+
                 "Nature Struct. Biol. 6, 458-463.", )

  def Interact_output(self, line, popen, empty):
    self.process.UpdateInfo(line)

  def TreatInput(self):
    self.par=self.inp.Get('datafile',typ='par')
    if not self.par:
      common.Error('No par file inputted to {0}'.format(self.name))
    self.SetArg(self.par.GetFileName())
    self.SetArg('1')

  def DefineOutput(self):
    with open(self.par.GetFileName()) as f:
      content=f.read()
      arp_dir=re.findall(r'set\s+WORKDIR\s+=\s+\'?(.+)\'?',content)[-1].strip()
      job=re.findall(r'set\s+JOB_ID\s+=\s+\'?(.+)\'?',content)[-1].strip()
    pdbfname=job+'_warpNtrace.pdb'
    mtzfname=job+'_warpNtrace.mtz'
    out_model = self.out.AddNew( 'model', os.path.join(arp_dir,pdbfname), typ='partial', xname=self.fsf.xname )
    out_map = self.out.AddNew( 'mapcoef', os.path.join(arp_dir,mtzfname), typ='weighted', xname=self.fsf.xname )
    out_map.SetLabel( ['f','ph','fom'], ['FWT','PHWT','FOM'] )
