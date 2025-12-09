#!/usr/bin/python
import os,sys,string
from process import process, crvapi
import common
par=common.parameter

class ref(process):
  name="macromolecular reciprocal space refinement"
  short_name="model refinement"
  supported_progs=["refmac"]
  supported_params = {}
  supported_params['target'] = par( desc='Experiment/refinement target', typ=str, cap=True )
  supported_params['cycles'] = par( desc='Number of model refinement cycles', typ=(int,bool) )
  supported_params['beta'] = par( desc='Bias parameter ("phcomb" only)', typ=float )
  supported_params['phcomb'] = par( desc='Refine with restrains to an inputted F and phase', typ=bool )
  supported_params['low_res_par'] = par( desc='Refine using low resolution specific parameters', typ=bool )
  supported_params['comb_with_typ'] = par( desc='Combine phases with the specified mapcoef type (only if "phcomb")', typ=str )
  supported_params['catch_output'] = par( desc='Catch program output continously', typ=bool )
  supported_params['no_substr'] = par( desc='Refine input partial model excluding substructure atoms', typ=bool )

  def TreatInOutPar(self, set_all_par=False):
    if not self.GetProg(supported=True):
      self.AddProg('refmac')
    if self.GetVirtPar('phcomb') and not self.GetVirtPar('comb_with_typ'):
      self.SetVirtPar('comb_with_typ', 'densmod')
    # determine the resolution of the relevant dataset against which we'll refine...
    obj,N=self.inp.GetTargetObjects(self.GetParam('target'),no_model=True,no_warn=True)
    if obj is None:
      # if no success then we try to assume Rice (this is useful eg in emulate mode when not all exists)
      obj,N=self.inp.GetTargetObjects('RICE',no_model=True,inten=None)
      if not obj:
        common.Error('The input objects for target Rice of {0} could not be retrieved'.format(self.name))
    ref_data = obj['f']  if 'f' in obj and obj['f']  else obj['f+']
    resol=ref_data.GetResolution(pro=self,accept_none=True)
    if (self.GetVirtPar('low_res_par') is not False and resol and resol>3.0):
      self.SetVirtPar('low_res_par')
    # changing the default number of cycles here in case of a separate refinement step 
    # (setting cycles to False restores the program default)
    if self.IsTrueOrNoneParam('cycles') and (not self.parent_process or self.parent_process.nick=='crank'):
      self.SetVirtPar('cycles', 15)
    # disabling hydrogens by default for low res
    if self.GetParam('low_res_par') and self.GetProg('refmac') and not self.GetProg('refmac').IsKey('make'):
      self.GetProg('refmac').SetKey('make',('hydrogens','no'))
    # setting default solvent params - testing only!
    #if self.GetProg(supported=True).nick=='refmac' and not self.GetProg(supported=True).IsKey('solvent'):
    #  self.GetProg(supported=True).SetKey('solvent',('yes','VDWProb','1.4','IONProb','0.8','RSHRink','0.8'))
    if self.GetParam('catch_output') is False:
      self.Info('Continous log output disabled.')
    if self.GetParam('no_substr') is None and self.inp.Get(has_atomtypes=True) and self.GetCrankParent():
      # if this refinement is the last crank subprocess, we do not want substr. atoms by default (finalizing)
      if self.inp.Get(has_atomtypes=True).GetAtomType() in ('S','SE') and self.GetParam('target') in ('MLHL','RICE') and \
         self is self.GetCrankParent().processes[-1]:
        self.SetParam('no_substr')
    elif not self.parent_process or self.parent_process.nick=='crank' or self.GetParam('catch_output'):
      self.GetProg(supported=True).interact_output=True
    process.TreatInOutPar(self,set_all_par)

  def RunPreprocess(self,*args,**kwargs):
    process.RunPreprocess(self,*args,**kwargs)
    modpart=self.inp.Get('model',typ='partial',filetype='pdb')
    if not self.GetParam('no_substr'):
      # merge partial model and substructure model into one if both are supplied for the same crystal
      if modpart and modpart is self.inp.Get('model',typ=('partial+substr','partial'),filetype='pdb'):
        modsubst=self.inp.Get('model',typ='substr',xname=modpart.xname,filetype='pdb')
        if not modsubst and (modpart.xname=='native' or 'native' in modpart.custom) and self.GetParam('target')=='SAD':
          modsubst=self.inp.Get('model',typ='substr',filetype='pdb')
        if modsubst:
          pdbmerge=self.AddProg('pdbmerge',propagate_inp=False,propagate_out=False)
          pdbmerge.inp.Set([modpart,modsubst])
          pdbmerge.SetKey('nomerge')
          pdbmerge.Run()
          self.inp.Add(pdbmerge.out.Get('model'))
          self.programs.remove(pdbmerge)
    elif not modpart or modpart is not self.inp.Get('model',filetype='pdb'):
      # make sure partial model will be used (without substructure)
      if not modpart and self.inp.Get('model',typ='partial+substr',filetype='pdb'):
        sepsub=self.AddProcess('sepsubstrprot',propagate_inp=False,propagate_out=False)
        sepsub.inp.Set(self.inp.Get('model',typ='partial+substr',filetype='pdb'))
        sepsub.Run()
        modpart=sepsub.out.Get('model',typ='partial')
      if modpart:
        self.inp.AddCopy(modpart)
    # a hack to prevent refmac from choking by some chain ID
    # refmac chokes by _ and \ , maybe more will be needed to add here...
    bad_ids = ('\\', '_', '\x7f', '\x80','\x81','\x82','\x83','\x84','\x85','\x86','\x87','\x88','\x2e','\x09','\x0a','\x0b','\x0c','\x0d', '"', '\'', '#', ';', '?')
    mod=self.inp.Get('model',filetype='pdb')
    if mod:
      pdbcur=self.AddProg('pdbcur')
      pdbcur.SetKey('summarise')
      pdbcur.inp.Set(mod)
      pdbcur.Run()
      chain_ids = pdbcur.GetStat('chain_ids')
      for bad_id in set(bad_ids).intersection(chain_ids):
        for new_id in string.printable:
          if new_id not in chain_ids and new_id not in bad_ids and new_id not in pdbcur.bad_chain_ids:
            pdbcur.SetKey('renchain', ("\""+bad_id+"\"", "\""+new_id+"\"") )
            chain_ids.append(new_id)
            break
      if pdbcur.IsKey('renchain'):
        pdbcur.Run()
        mod=self.inp.Add(pdbcur.out.Get('model'))
      if set(pdbcur.bad_chain_ids).intersection(chain_ids).intersection(bad_ids):
        # rename substr. chains if ids that pdbcur does not accept are present
        pdbset=self.AddProg('pdbset')
        pdbset.inp.Set(mod)
        for bad_id in set(pdbcur.bad_chain_ids).intersection(chain_ids).intersection(bad_ids):
          new_id = next(nid for nid in string.printable if nid not in chain_ids and nid not in bad_ids and nid not in pdbcur.bad_chain_ids and nid not in pdbset.bad_chain_ids)
          pdbset.SetKey('chain',(bad_id,new_id))
          chain_ids.append(new_id)
        pdbset.Run()
        self.inp.Add(pdbset.out.Get('model'))
        self.programs.remove(pdbset)
      self.programs.remove(pdbcur)
    # low res par
    if self.GetVirtPar('low_res_par') and self.GetProg().nick=='refmac':
      refmac=self.GetProg()
      self.save_refmac=refmac.BackupAnyPars()
      if not refmac.IsKey('RIDG') or (refmac.GetKey('RIDG') and not 'DIST' in refmac.GetKey('RIDG',as_list=True)):
        refmac.SetKey('RIDG',['DIST','SIGM','0.1'])
      #if not refmac.IsKey('MAPC') or not 'SHAR' in refmac.GetKey('MAPC'):
        #if (self.parent_process.nick!='dmfull' or self.parent_process.GetParam('dmcyc')==self.parent_process.cyc) \
        #   and self.parent_process.nick!='refatompick':
        #  #self.GetProg().SetKey('MAPC',['SHAR',gcx.GetStat('wilson_B')/2])
        #  refmac.SetKey('MAPC', ['SHAR','ALPH','0.2'])
      # Refmac NCSR LOCAL (as of 5.8.0155) crashes if non-numerical insertion codes are used (generated by Buccaneer if there are many chains)...
      if not refmac.IsKey('NCSR'): 
        mon_obj=self.inp.Get('sequence',has_monomers_asym=True)
        if mon_obj:
          mon_obj.GetSequenceString()
        else:
          matthews=self.AddProcess('matthews', propagate_out=False)
          matthews.Run()
          mon_asym=matthews.GetProg().GetStat('monomers_asym')
          self.processes.remove(matthews)
        if (mon_obj and (mon_obj.monomers_asym>1 or mon_obj.seq_monomers>1)) or (not mon_obj and mon_asym>1):
          refmac.SetKey('NCSR', 'LOCAL')

  def RunPostprocess(self,restore=True,*args,**kwargs):
    refmac=self.GetProg()
    if hasattr(self,'save_refmac'):
      refmac.RestoreAnyPars(self.save_refmac)
    self.R, self.Rfree = refmac.GetStat('rfact')[-1], (refmac.GetStat('rfree', accept_none=True) or [None])[-1]
    self.fom=refmac.GetStat('fom')[-1]
    self.report_R, self.report_Rfree = self.R, self.Rfree
    if self.parent_process.nick != 'dmfull':
      self.Info('{0} after refinement is {1}'.format(refmac.stat['rfact'].name, self.R))
      if self.rv_report:
        self.rv_report.Text("<i><small>Result:</small></i>")
        self.rv_report.Text('&emsp;'+'Final R factor is {0}'.format(self.R))
      rfree=refmac.GetStat('rfree',accept_none=True)
      if rfree:
        self.Info('{0} after refinement is {1}'.format(refmac.stat['rfree'].name, rfree[-1]))
        if self.rv_report:
          self.rv_report.Text('&emsp;Final Rfree factor is {0}'.format(rfree[-1]))
    # make sure that there is an output substructure - this may be removed later if not needed
    if self.out.Get('model',typ='partial+substr') and not self.out.Get('model',typ='substr'):
      separate_models=self.AddProcess('sepsubstrprot')
      separate_models.inp.Set(self.out.Get('model',typ='partial+substr'))
      # filter on input substr. chains, to make sure that SE/S built by Buccaneer rebuilding is not included in substr.
      substr=self.inp.Get('model',typ='substr',has_atomtypes=True)
      if substr and substr.GetAtomType in ('S','SE'):
        separate_models.inp.Add(self.inp.Get('model',typ='substr'))
      separate_models.Run()
      self.processes.remove(separate_models)
    if self.rv_report:
      self.rv_plotref = self.rv_report.Plot( 'R factor vs refinement cycle', "Cycle", "R factors", intx=True )
      Rs = refmac.GetStat('rfact')
      cyc = [i for i in range(len(Rs))]
      self.rv_R = self.rv_plotref.PlotLine(["x",cyc],["R",Rs])
      if rfree:
        self.rv_R = self.rv_plotref.PlotLine(["x",cyc],["Rfree",rfree])
      if not self.parent_process or not self.parent_process.ccp4i2:
        out_pdb = self.out.Get('model',typ=('partial','partial+substr')).GetFileName()
        out_mtz = crvapi.SetMtzFWT(self,self.out.Get('mapcoef',typ='combined'))
        #out_map = self.out.Get('mapcoef',typ='weighted',filetype='map',conv_opts=['no_rewrite']).GetFileName('map')
        #out_dmap = self.out.Get('mapcoef',typ='diff',filetype='map',conv_opts=['no_rewrite']).GetFileName('map')
        self.rv_report.Text("<BR><i><small>Output:</small></i>")
        self.rv_files = self.rv_report.DataFile(out_pdb, "xyz", "<small>Refined model and map</small>")
        self.rv_files.DataFile(out_mtz, "hkl:map", flush=True)
        #self.rv_files.DataFile(out_map, "hkl:ccp4_map")
        #self.rv_files.DataFile(out_dmap, "hkl:ccp4_dmap", flush=True)
    process.RunPostprocess(self,restore,*args,**kwargs)

  def UpdateInfo(self,line,popen):
    # copies over the log (to be interpreted as loggraph)
    if self.opened_loggraph:
      self.LGInfo(line)
