#!/usr/bin/python
import os,sys,copy
from process import process
from program import program
import common,data

class faest(process):
  name="FA estimation and/or other substructure detection preparations"
  short_name="FA estimation"
  supported_progs=["afro","shelxc","ecalc"]
  supported_params = {}
  supported_params['target'] = common.parameter( desc='Experiment/refinement target', typ=str, cap=True )
  supported_params['bfactor'] = common.parameter( desc='set B factor of atoms', typ=(float,bool) )

  def TreatInOutPar(self, set_all_par=False):
    #self.GetProg(supported=True).SetKey('MULT')
    # just to prevent a useless warning of parameter not used for shelx
    process.TreatInOutPar(self,set_all_par)
    obj,N=self.inp.GetTargetObjects(self.GetParam('target'), no_warn=True)
    if obj is None:
      obj,N=self.inp.GetTargetObjects(self.GetParam('target'), inten=True)
    # expected number of substr. atoms should have been guessed already in manager?
    if obj and obj['mod'] and obj['mod'].exp_num_atoms is None:
      obj['mod'].GuessNumSubstrAtomsFromSeq(self)
      # passing output to be used by substrdet
      self.out.Add(obj['mod'])
    if self.GetProg('ecalc') or self.GetProg('afro'):
      # this will be used in the emulation mode so that the GUI gets the defaults for substr. cutoff
      # if too slow then this can be moved to RunPreprocess()
      fp_obj=self.inp.Get('fsigf',typ='plus',col='f')
      fm_obj=self.inp.Get('fsigf',typ='minus',col='f')
      fa_obj=self.inp.Get('fsigf',typ='fa',filetype='mtz', inp_cont=self.inp.Get('fsigf',typ='average',filetype='mtz',try_convert=False))
      if not fa_obj:
        fa_obj=self.inp.Get('fsigf',typ='delta-anom',filetype='mtz')
      if not fa_obj and fp_obj and fm_obj:
        mtzmadmod=self.AddProg('mtzMADmod',propagate_out=False)
        mtzmadmod.Run(check_bin=True)
        fa_obj = mtzmadmod.out.Get(filetype='mtz',typ='delta-anom')
        self.inp.Add(fa_obj)
      # exlusion done directly in Afro now
      if self.GetProg('ecalc'):
        for ob in self.ExcludeRefs(fa_obj,fp_obj,fm_obj):
          if ob:  self.inp.Add( ob )
        return
    cp = self.GetCrankParent()
    if not self.GetProg(supported=True):
      if program.from_name('shelxc',None).CheckBinary(silent=True) and (not cp or not cp.jscofe or self.GetParam('target')!='SAD'):
        self.AddProg('shelxc')
      else:  # afro+prasa is default for sad at cloud (only) now
        self.AddProg('afro')
    # pass shelxd keywords to shelxc 
    # first create copies of crank, substrdet, shelxd and run treatpar to make sure all required keys are passed
    # this is not ideal but the very idea of shelxc generating a script for shelxd is far from ideal...
    if self.GetProg('shelxc') and self.GetCrankParent() and self in self.GetCrankParent().processes:
      for sd in [sd for sd in cp.GetProcesses('substrdet') if sd in cp.processes[cp.processes.index(self):]]:
        crank = copy.copy(cp)
        crank.processes = cp.processes[:]
        substrdet = crank.AddProcessCopy( sd, propagate_inp=True )
        # pass the adjustments done to obj['mod'] in this routine (not needed anymore?)
        if obj and obj['mod'] and obj['mod'] not in substrdet.inp.GetAll('model'):
          substrdet.inp.Add(obj['mod'])
        crank.SetTargetDefault(substrdet,emulate=True)
        substrdet.TreatInOutPar(skip_cutoffest=True)
        if substrdet.GetProg('shelxd'):
          shelxd = substrdet.AddProgCopy( sd.GetProg('shelxd'), propagate_inp=True )
          shelxd.TreatParams()
          for key in shelxd.key_list:
            if not self.GetProg('shelxc').GetKey(key):
              self.GetProg('shelxc').AddToKey(key, shelxd.GetKey(key,allval=True))
    # get&set bfactor for afro
    if self.GetProg('afro') and self.IsTrueOrNoneParam('bfactor'):
      gcx=self.AddProg('gcx', propagate_out=False)
      try:
        gcx.Run()
      except (gcx.ProgramRunError,):
        wilsonB = 50.
      else:
        wilsonB = gcx.GetStat('wilson_B')
      self.programs.remove(gcx)
      self.SetParam('bfactor', wilsonB)

#  def RunPreprocess(self,*args,**kwargs):
#    process.RunPreprocess(self,*args,**kwargs)
#    if self.GetProg('ecalc'):
#      fp_obj=self.inp.Get('fsigf',typ='plus')
#      fm_obj=self.inp.Get('fsigf',typ='minus')
#      fa_obj=self.inp.Get('fsigf',typ='fa',filetype='mtz', inp_cont=self.inp.Get('fsigf',typ='average',filetype='mtz',try_convert=False))
#      if not fa_obj:
#        fa_obj=self.inp.Get('fsigf',typ='delta-anom',filetype='mtz')
#      if not fa_obj:
#        mtzmadmod=self.AddProg('mtzMADmod',propagate_out=False)
#        mtzmadmod.Run(check_bin=True)
#        fa_obj = mtzmadmod.out.Get(filetype='mtz',typ='delta-anom')
#        self.inp.Add(fa_obj)
#      self.ExcludeRefs(fa_obj,fp_obj,fm_obj)

  def ExcludeRefs(self,fa_obj,fp_obj,fm_obj):
    # do not try to determine if called only to check binaries...
    if self.GetCrankParent() and hasattr(self.GetCrankParent(),'check_binaries'):  return None,None,None
    if fp_obj and fm_obj and fa_obj:
      sft=self.AddProg('sftools')
      sft_outmtz = 'sft_adjusted.mtz'
      sft.runname=sft.name+'_sig'
      sft.SetKey('read', ('"'+fa_obj.GetFileName('mtz')+'"','col','"'+fa_obj.GetLabel('f')+'"','"'+fa_obj.GetLabel('sigf')+'"'))
      sft.SetKey('read', ('"'+fp_obj.GetFileName('mtz')+'"','col','"'+fp_obj.GetLabel('f')+'"','"'+fp_obj.GetLabel('sigf')+'"'))
      sft.SetKey('read', ('"'+fm_obj.GetFileName('mtz')+'"','col','"'+fm_obj.GetLabel('f')+'"','"'+fm_obj.GetLabel('sigf')+'"'))
      sft.SetKey('calc', ('R','col','_minpm_','=','col','"'+fp_obj.GetLabel('f')+'"','col','"'+fm_obj.GetLabel('f')+'"','min'))
      sft.SetKey('calc', ('R','col','_delpm_','=','col','"'+fp_obj.GetLabel('f')+'"','col','"'+fm_obj.GetLabel('f')+'"','-'))
      sft.SetKey('calc', ('R','col','_delpm_','=','col','_delpm_','abs'))
      sft.SetKey('calc', ('R','col','_delratio_','=','col','_delpm_','col','_minpm_','/'))
      sft.SetKey('select', ('col','_delratio_','<','1.0'))
      sft.SetKey('calc', ('R','col','_sigpratio_','=','col','"'+fp_obj.GetLabel('f')+'"','col','"'+fp_obj.GetLabel('sigf')+'"','/'))
      sft.SetKey('calc', ('R','col','_sigmratio_','=','col','"'+fm_obj.GetLabel('f')+'"','col','"'+fm_obj.GetLabel('sigf')+'"','/'))
      #sft.SetKey('select', ('col','_sigpratio_','>','2.5'))
      #sft.SetKey('select', ('col','_sigmratio_','>','2.5'))
      sft.SetKey('calc', ('R','col','_sigpmratio_','=','col','"'+fp_obj.GetLabel('sigf')+'"','col','"'+fm_obj.GetLabel('sigf')+'"','/'))
      sft.SetKey('select', ('col','_sigpmratio_','>','0.333'))
      sft.SetKey('select', ('col','_sigpmratio_','<','3.'))
      sft.SetKey('write', sft_outmtz)
      sft.SetKey('quit\nY')
      sft.Run()
      fa_obj=sft.out.AddCopy(fa_obj)
      fa_obj.InitLabels({'f':fa_obj.GetLabel('f'),'sigf':fa_obj.GetLabel('sigf')})
      sft.out.SetFileToChild(fa_obj,sft_outmtz)
      fp_obj=sft.out.AddCopy(fp_obj)
      fp_obj.InitLabels({'f':fp_obj.GetLabel('f'),'sigf':fp_obj.GetLabel('sigf')})
      sft.out.SetFileToChild(fp_obj,sft_outmtz)
      fm_obj=sft.out.AddCopy(fm_obj)
      fm_obj.InitLabels({'f':fm_obj.GetLabel('f'),'sigf':fm_obj.GetLabel('sigf')})
      sft.out.SetFileToChild(fm_obj,sft_outmtz)
      self.programs.remove(sft)
    else:
      common.Warning('Exclusion of reflections could not be done due to F+/F-/FA not supplied.')
      return None,None,None
    return fa_obj,fp_obj,fm_obj

  def PlotInterpolate(self,c,r,high_r):
    c_interp, c, r = [], [float(i) for i in c if i is not None], [float(i) for i in r if i is not None]
    for hr in high_r:
      i = next((i for i,cr in enumerate(r[:-1]) if cr>=hr and r[i+1]<=hr), None)
      c_interp.append( c[i]+(c[i+1]-c[i])*(r[i]-hr)/(r[i]-r[i+1]) if i is not None else None )
    return c_interp

  def PrintActualLogGraph(self):
    prog=self.GetProg(supported=True)
    reso,dsigd,cc_half=None,None,None
    if prog.nick == 'afro':
      name=['SAD']
      reso =  [list(filter(None,prog.GetStat("resol_bins")))]
      high_res_anom = reso[0][-1]
      dsigd = [list(filter(None,prog.GetStat("dano_bins")))]
    if prog.nick == 'shelxc':
      name=prog.GetStat("data_set")
      if prog.GetStat("resol_bins"):
        reso = [list(filter(None,r)) for r in prog.GetStat("resol_bins")]
        high_res = min((min(res,key=float) for res in reso if res),key=float)
        high_res_anom = min((min(res,key=float) for res,n in zip(reso,name) if res and n!='NAT'),key=float)
        high_r = [float(i) for i in next(r for r in reso if high_res in r) if i is not None]
        high_r_anom = [float(i) for i in next(r for r in reso if high_res_anom in r) if i is not None]
        dsigd = [list(filter(None,d)) for d in prog.GetStat("dano_bins")]
        dsigd = list(filter(None,dsigd))
        cc_half = [list(filter(None,c)) for c in prog.GetStat("cc1/2")]
        if 'NAT' in name and (not 'SIRA' in name or len(name)>len(dsigd)):
          dsigd.insert( name.index('NAT'),() )
          if cc_half:
            cc_half.insert( name.index('NAT'),() )
        # this line should be moved elsewhere as it is not output specific
        self.cc_half_plot=(reso,cc_half,name)
        if self.opened_loggraph and dsigd:
          self.GetLogGraphHandle().seek(0,0)
          self.LGInfo(self.GetCCP4Header())
          self.LGInfo(self.GetLogGraphPrograms())
          self.LGInfo('\n $TABLE : Anomalous signal to noise vs resolution:')
          self.LGInfo('$GRAPHS    : Dano/Sigdano vs resolution :A:1,2:\n$$')
          self.LGInfo('Resolution Dano/Sigdano $$ $$')
          for i in range(min(len(reso[0]),len(dsigd[0]))):
            self.LGInfo('{0} {1}'.format(reso[0][i],dsigd[0][i]))
          self.LGInfo('$$\n')
    if reso:
        if self.rv_report is not None:
          if dsigd:
            self.rv_plot1 = self.rv_report.Plot( 'Dano/Sigdano', "Resolution [Angstrom]", block="Data stats vs resolution", legendloc='ne', ymin=0. )
            for r,d,n in zip(reso,dsigd,name):
              if n!='NAT':
                if high_res_anom in r:
                  self.rv_plot1.PlotLine( ["x",r], [n,d], custom_x_tick=True )
                else:
                  d_interp = self.PlotInterpolate(d,r,high_r_anom)
                  self.rv_plot1.PlotLine( ["x",[i+1 for i in range(len(d_interp))]], [n,d_interp] )
                if not cc_half:  # very primitive interpretation but should be ok for now; only printed if cchalf not available
                  message="The anomalous signal for {} appears to be very weak or not present at all.".format(n)
                  if any(dv for dv,rv in zip(d,r) if float(rv)<6 and float(dv)>1.0): message="The anomalous signal for {} appears to be weak.".format(n)
                  if any(dv for dv,rv in zip(d,r) if float(rv)<6 and float(dv)>1.2): message="The anomalous signal for {} appears to be rather weak but possibly sufficient.".format(n)
                  if any(dv for dv,rv in zip(d,r) if float(rv)<5 and float(dv)>1.5): message="The anomalous signal for {} appears to be rather strong.".format(n)
                  if any(dv for dv,rv in zip(d,r) if float(rv)<5 and float(dv)>2.0): message="The anomalous signal for {} appears to be strong.".format(n)
                  self.rv_report.Text(message)
            self.rv_shadow1 = self.rv_plot1.PlotLine( ["x",(0,len(r)+1)], ["Very small or no anom. signal",(0.8,0.8)], fill=1 )
            if prog.nick == 'shelxc':
             comp = [list(filter(None,c)) for c in prog.GetStat("complete")]
             if comp:
              self.rv_plot2 = self.rv_plot1.parent.Plot( 'Completeness [%]', "Resolution [Angstrom]", legendloc='nw')
              for r,c,n in zip(reso,comp,name):
                if high_res in r:
                  self.rv_plot2.PlotLine( ["x",r], [n,c], custom_x_tick=True )
                else:
                  c_interp = self.PlotInterpolate(c,r,high_r)
                  self.rv_plot2.PlotLine( ["x",[i+1 for i in range(len(c_interp))]], [n,c_interp] )
             iovers = [list(filter(None,Is)) for Is in prog.GetStat("Ioversigma")]
             if iovers:
              self.rv_plot3 = self.rv_plot1.parent.Plot( 'I/SigI', "Resolution [Angstrom]", legendloc='ne')
              for r,Is,n in zip(reso,iovers,name):
                if high_res in r:
                  self.rv_plot3.PlotLine( ["x",r], [n,Is], custom_x_tick=True )
                else:
                  Is_interp = self.PlotInterpolate(Is,r,high_r)
                  self.rv_plot3.PlotLine( ["x",[i+1 for i in range(len(Is_interp))]], [n,Is_interp] )
             cresol = list(filter(None, prog.GetStat("resol_bins_cor")))
             if cresol:
              self.rv_plot4 = self.rv_plot1.parent.Plot( 'Signed anom diff data correlations', "Resolution [Angstrom]", \
                                "Signed Dano correlation between datasets", legendloc='ne', ymin=0.0)
              for n1 in name:
                for n2 in name:
                  correl = prog.GetStat("correl",(n1,n2),accept_none=True)
                  if correl:
                    correl = list(filter(None,correl))
                    self.rv_plot4.PlotLine( ["x",cresol], [ n1+'/'+n2, correl ], custom_x_tick=True )
              #self.rv_shadow4 = self.rv_plot4.PlotLine( ["x",(-1,len(r)+1)], ["Very small or no anom. signal",(30,30)], fill=1 )
             if cc_half:
              self.rv_plot5 = self.rv_plot1.parent.Plot( 'CCanom1/2', "Resolution [Angstrom]", legendloc='ne', ymin=0.0 )
              for r,cc,n in zip(reso,cc_half,name):
                if n!='NAT':
                  if high_res_anom in r:
                    self.rv_plot5.PlotLine( ["x",r], [n,cc], custom_x_tick=True )
                  else:
                    cc_interp = self.PlotInterpolate(cc,r,high_r_anom)
                    self.rv_plot5.PlotLine( ["x",[i+1 for i in range(len(cc_interp))]], [n,cc_interp] )
                  # very primitive interpretation but should be ok for now
                  message="The anomalous signal for {} appears to be very weak or not present at all.".format(n)
                  if any(cv for cv,rv in zip(cc,r) if float(rv)<6 and float(cv)>25): message="The anomalous signal for {} appears to be weak.".format(n)
                  if any(cv for cv,rv in zip(cc,r) if float(rv)<6 and float(cv)>35): message="The anomalous signal for {} appears to be rather weak but possibly sufficient.".format(n)
                  if any(cv for cv,rv in zip(cc,r) if float(rv)<5 and float(cv)>50): message="The anomalous signal for {} appears to be rather strong.".format(n)
                  if any(cv for cv,rv in zip(cc,r) if float(rv)<4.5 and float(cv)>75): message="The anomalous signal for {} appears to be strong.".format(n)
                  self.rv_report.Text(message)
              self.rv_shadow5 = self.rv_plot5.PlotLine( ["x",(-1,len(r)+1)], ["Very small or no anom. signal",(30,30)], fill=1 )


  def RunPostprocess(self,restore=True,*args,**kwargs):
    self.PrintActualLogGraph()
    # check whether the cif file exists (eg not outputted by older shelx versions)
    prog=self.GetProg(supported=True)
    if prog.nick=='shelxc':
      for o in prog.out.GetAll(filetype='cif'):
        if not os.path.isfile(o.GetFileName('cif')):
          prog.out.Delete(o)
          if self.GetParam('target')=='SAD':
            common.Warning('CIF file with anomalous data not outputted by {0}.'.format(prog.name))
    # check whether any anomalous data is present
    fa_mtz = prog.out.Get('fsigf',typ=('fa'),col='f',filetype='mtz')
    if fa_mtz:
      sft=self.AddProg('sftools')
      sft.SetKey('read','"'+fa_mtz.GetFileName()+'"')
      sft.SetKey('checkhkl')
      sft.Run()
      if sft.GetStat('all_absent',fa_mtz.GetLabel('f'),accept_none=True) or sft.GetStat('all_zero',fa_mtz.GetLabel('f'),accept_none=True):
        common.Error('No anomalous signal found. Check your data.')
    # use CC1/2 cutoff by default if unmerged data inputted and if shelxc was run
    # (perhaps this should be in manager?)
    thres=25.
    crank,substrdet=self.GetCrankParent(),None
    if crank:  substrdet=crank.GetProcess('substrdet')
    if substrdet:  substrdet.TreatInOutPar(skip_cutoffest=True)
    if substrdet and substrdet.GetParam('high_res_cutoff_cchalf'):
      an_i = next((i for i,n in enumerate(self.cc_half_plot[2]) if n!='NAT'),-1)  if \
          hasattr(self,'cc_half_plot') and self.cc_half_plot and filter(None,self.cc_half_plot[1]) else -1
      if an_i>=0:
        reso, cc_half = list(map(float,self.cc_half_plot[0][an_i])), list(map(float,self.cc_half_plot[1][an_i]))
        #reso = map( float,[filter(None,r) for r in prog.GetStat("resol_bins")][0] )
        #cc_half = map( float,[filter(None,c) for c in prog.GetStat("cc1/2")][0] )
        cutoff,cutoff0,cc,cc0 = reso[0],reso[0],cc_half[0],cc_half[0]
        for r,c in zip(reso,cc_half):
          # stop at two subsequent bins below the threshold
          if c>thres:
            cutoff,cutoff0,cc,cc0 = r,r,c,c
          else:
            if abs(cutoff0-cutoff)>1e-5:
              break
            cutoff0,cc0 = r,c
        # linear interpolation
        if abs(cutoff0-cutoff)>1e-5:
          cutoff -= float("{0:.2f}".format( (cc-thres)/(cc-cc0) * (cutoff-cutoff0) ))
        if cutoff>5.5:
          self.Info('CCanom1/2({}) cutoff too large ({}A), will not be set from CCanom1/2.'.format(thres,cutoff))
        else:
          crank.prep.GetProcess('substrdet').SetParam('high_res_cutoff', cutoff)
          self.Info('CCanom1/2({}) cutoff: {}'.format(thres,cutoff))
      else:
        common.Warning('CCanom1/2 could not be determined')
    # this assumes that program is specified.  if not and if prasa is default then this is dealt with directly in substrdet
    if substrdet and prog.nick=='afro' and substrdet.GetProg(supported=True) and substrdet.GetProg(supported=True).nick=='prasa' and not substrdet.IsParam('afroprasa'):
      crank.prep.GetProcess('substrdet').SetParam('afroprasa',True)
    elif substrdet and prog.nick!='afro':
      crank.prep.GetProcess('substrdet').SetParam('afroprasa',False)
    process.RunPostprocess(self,restore,*args,**kwargs)
