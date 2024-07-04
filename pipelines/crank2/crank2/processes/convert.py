#!/usr/bin/python
import os,sys
from process import process
import common


class convert(process):
  name="data format conversion"
  short_name="data conversion"
  supported_progs=["fft", "invfft", "mtz2various", "sftools", "convert2mtz", "f2mtz","shelxc"]
  supported_params = {}
  supported_params['outputformat'] = common.parameter( desc='Output format for mtz2various conversion', typ=str )
  supported_params['inputformat'] = common.parameter( desc='Input format for convert2mtz/f2mtz conversion', typ=str )
  no_reporting=True
  # setup for the conversion - filled automatically when creatin the convert process instance using RegisterConv
  setup = {}
#  opts = {}


  def RegisterConv(converts,func,fromlist,tolist):
    for f in fromlist:
      for t in tolist:
        converts[(f,t)] = func

  def ConvertMTZ2Various(self,oldtyp,newtyp):
    self.SetParam('outputformat',newtyp)
    # lets convert hl to phases if only hl are available and phs output asked
    if newtyp=='phs' and not self.inp.Get('mapcoef',filetype='mtz',col=('fom','ph')):
      hl=self.inp.Get('mapcoef',filetype='mtz',col='hla')
      if hl:
        hl2fom=self.AddProg('chltofom',propagate_out=False)
        hl2fom.Run()
        self.programs.remove(hl2fom)
      else:
        return
    mtz2var=self.AddProg('mtz2various')
    mtz2var.outfilename[newtyp] = self.inp.Get(filetype='mtz').GetFileName(trim_path=True).rsplit('.',1)[0]+'_'+mtz2var.name[:5]+'.'+newtyp
    # automatic scaling does not work correctly in mtz2various for example if max(I-)>max(I+)
    # as a workaround, find and set the scale here to prevent overflows with hkl.
    mtz_ip=self.inp.Get('fsigf',filetype='mtz',col='i',typ='plus')
    mtz_im=self.inp.Get('fsigf',filetype='mtz',col='i',typ='minus')
    if newtyp in ('hkl','sca') and not self.GetProg().GetKey('SCAL') and mtz_ip and mtz_im and not 'amplitudes' in self.opts:
      try:
        import gemmi
      except ImportError: # may occasionally fail if the values are so large that even sftools reports just stars
        sftools=self.AddProg('sftools')
        sftools.SetKey("read", ['"'+mtz_ip.GetFileName('mtz')+'"', 'col', "\"{0}\"".format(mtz_ip.GetLabel('i'))])
        sftools.SetKey("checkhkl")
        sftools.Run()
        maxp = sftools.GetStat('max', mtz_ip.GetLabel('i')[:12])
        minp = sftools.GetStat('min', mtz_ip.GetLabel('i')[:12])
        sftools.SetKey("read", ['"'+mtz_im.GetFileName('mtz')+'"', 'col', "\"{0}\"".format(mtz_im.GetLabel('i'))])
        sftools.SetKey("checkhkl")
        sftools.Run()
        scal = 99999./max(maxp,abs(minp),sftools.GetStat('max', mtz_im.GetLabel('i')[:12]),abs(sftools.GetStat('min', mtz_im.GetLabel('i')[:12])))
        self.programs.remove(sftools)
      else:
        mtz = gemmi.read_mtz_file(mtz_ip.GetFileName('mtz'))
        ip=mtz.column_with_label(mtz_ip.GetLabel('i'))
        im=mtz.column_with_label(mtz_im.GetLabel('i'))
        scal = 99999./max(ip.max_value,abs(ip.min_value),im.max_value,abs(im.min_value))
      self.GetProg().SetKey('SCAL', scal)
  RegisterConv(setup, ConvertMTZ2Various, ['mtz'], ['hkl','sca','phs'])

  def ConvertMTZ2drear(self,oldtyp,newtyp):
    out_filename="convert.drear"
    sftools=self.AddProg('sftools')
    fa=self.inp.Get('fsigf',typ='fa',filetype='mtz')
    sftools.SetKey("read", '"'+fa.GetFileName('mtz')+'"')
    sftools.SetKey("write", (out_filename, "format(3I5,2F9.2,2F9.3)", "col", \
      fa.GetLabel('f'), fa.GetLabel('sigf'), fa.GetLabel('e'), fa.GetLabel('sige')))
    sftools.SetKey("quit")
    self.out.Add(fa)
    self.out.AddFileToChild(fa,out_filename,newtyp)
  RegisterConv(setup, ConvertMTZ2drear, ['mtz'], ['drear'])

  ### from reasons yet unknown does not work for conversion of shelc's fa's
  def ConvertVarious2MTZ(self,oldtyp,newtyp):
    if self.inp.Get('fsigf',typ='fa'):
      self.parent_process.ConvertF2MTZ(oldtyp,newtyp)
    else:
      self.convert_proc.SetParam('inputformat',oldtyp)
      conv2mtz=self.AddProg('convert2mtz')
      conv2mtz.SetKey('no-complete')
      hklin=self.inp.Get('fsigf',filetype='hkl')
      if oldtyp=='hkl' and hklin:
        mtzin=self.inp.Get('fsigf',typ='average',filetype='mtz',try_convert=False)
        if mtzin:
          conv2mtz.SetKey('spacegroup', mtzin.GetSpacegroup(conv2mtz).replace(" ",""))
          conv2mtz.SetKey('cell', ','.join(str(c) for c in mtzin.GetCell(conv2mtz)))
        else:
          self.no_conversion=True
    # automatic scaling does not work in convert2mtz, we will need to find and set the scale here 
    # to prevent overflows with hkl. this can be based on the code from mtz2various.
  #RegisterConv(setup, ConvertVarious2MTZ, ['hkl','sca','phs'], ['mtz'])

  def ConvertF2MTZ(self,oldtyp,newtyp):
    # f2mtz supports average, fa and phs as of now
    self.SetParam('inputformat',oldtyp)
    oo=self.inp.Get(filetype=oldtyp,try_convert=False)
    # check for unmerged hkl
    unmerged=False
    if oldtyp=='hkl' and oo.GetType() in ('average',):
      with open(oo.GetFileName('hkl')) as f:
        prev_line = '               '
        for il,line in enumerate(f):
          if len(line)>20:
            if line[:12]==prev_line[:12] and line[:12]!='   0   0   0':
              unmerged=True
              break
            if il>5000:  break
            prev_line=line
      if unmerged:
        shelxc=self.AddProg('shelxc',propagate_out=False)
        shelxc.rundir=os.path.join(self.rundir,'shel_'+oo.GetCrystalName()+oo.GetDataName())
        shelxc.Run()
        self.inp.Add(shelxc.out.Get('fsigf',filetype='hkl'))
    mtzin=self.inp.Get('fsigf',typ='average',filetype='mtz',col='f',xname=oo.xname,dname=oo.dname,try_convert=False)
    if not mtzin and (oo.GetType()=='fa' or oldtyp=='phs'):
      mtzin=self.parent_process.inp.Get('fsigf',typ='average',filetype='mtz',col='f',xname=oo.xname,dname=oo.dname,try_convert=False)
      if mtzin:
        self.inp.AddCopy(mtzin)
    #if mtzin and oo.GetType()=='fa':
    f2mtz=self.AddProg('f2mtz')
    f2mtz.Run()
    if f2mtz.out.Get(filetype='mtz'):
        sftools=self.AddProg('sftools')
        sftools.inp.Add(f2mtz.out.Get(filetype='mtz'))
        sftools.SetKey("read", '"'+oo.GetFileName(oldtyp)+'"')
        sftools.SetKey("reduce")
        sftools.SetKey("quit")
        self.programs.remove(f2mtz), self.programs.remove(sftools)
    # this is not really correct: the same object shared by unmerged hkl and merged mtz - a distinction should be made!
    if unmerged:
      self.inp.AddFileToChild(oo, f2mtz.out.Get(filetype='mtz').GetFileName('mtz'), 'mtz')
    #elif oldtyp=='hkl':
    #  self.ConvertVarious2MTZ(oldtyp,newtyp)
      #self.no_conversion=True
    #else:
    #  f2mtz=self.AddProg('f2mtz')
  RegisterConv(setup, ConvertF2MTZ, ['hkl','phs'], ['mtz'])

  def ConvertSCAHKL2MTZ(self,oldtyp,newtyp):
    oo=self.inp.Get(filetype=oldtyp,try_convert=False)
    shelxc=self.AddProg('shelxc',propagate_out=False)
    shelxc.rundir=os.path.join(self.rundir,'shel_'+oo.GetCrystalName()+oo.GetDataName())
    if oo.GetType()=='minus': 
      shelxc.ExternalRun( [shelxc.binary,], [] )
      shelxc.default_exper='SAD'
    if oo.GetType()=='average' or (oo.GetType()=='minus' and shelxc.GetStat('version')>='2015/1'):
      shelxc.Run()
      if oo.GetType()=='average':
        hkl_out=shelxc.out.Get('fsigf',filetype='hkl')
        self.inp.AddFileToChild(oo, hkl_out.GetFileName(), 'hkl')
        self.ConvertF2MTZ('hkl','mtz')
      elif oo.GetType()=='minus':
        cif_out=shelxc.out.Get('fsigf',typ='minus',filetype='cif')
        self.inp.AddFileToChild(oo, cif_out.GetFileName(), 'cif')
        oo.SetLabel( ('i','sigi'), (cif_out.GetLabel('i'),cif_out.GetLabel('sigi')) )
        oo_plus=self.inp.Get(typ='plus',filename=oo.GetFileName(oldtyp))
        if oo_plus:
          cif_out=shelxc.out.Get('fsigf',typ='plus',filetype='cif')
          self.inp.AddFileToChild(oo_plus, cif_out.GetFileName(), 'cif')
          oo_plus.SetLabel( ('i','sigi'), (cif_out.GetLabel('i'),cif_out.GetLabel('sigi')) )
        self.ConvertCIF2MTZ('cif','mtz')
    self.programs.remove(shelxc)
  RegisterConv(setup, ConvertSCAHKL2MTZ, ['sca','HKL'], ['mtz'])

  def ConvertMTZ2Map(self,oldtyp,newtyp):
    fft=self.AddProg('fft')
    if 'outfilename' in self.opts:
      fft.outfilename['map'] = self.opts[1+self.opts.index('outfilename')]
  RegisterConv(setup, ConvertMTZ2Map, ['mtz'], ['map'])

  def ConvertMap2MTZ(self,oldtyp,newtyp):
    invfft=self.AddProg('invfft')
    if 'outfilename' in self.opts:
      invfft.outfilename['mtz'] = self.opts[1+self.opts.index('outfilename')]
    if not invfft.inp.Get('fsigf',filetype='mtz',col='f',try_convert=False):
      if self.parent_process:
        f = self.parent_process.inp.Get('fsigf', filetype='mtz', col='f', typ='average', \
          xname=self.inp.Get(filetype=oldtyp).GetCrystalName())
        if f:
          invfft.inp.Add(f)
  RegisterConv(setup, ConvertMap2MTZ, ['map'], ['mtz'])

  def ConvertCIF2MTZ(self,oldtyp,newtyp):
    cif2mtz=self.AddProg('cif2mtz')
  RegisterConv(setup, ConvertCIF2MTZ, ['cif'], ['mtz'])

  def ConvertPDB2RES(self,oldtyp,newtyp):
    ins=self.inp.Get('datafile', filetype='ins')
    if ins:
      self.SetParam('inputformat',oldtyp), self.SetParam('outputformat','frac')
      coordconv=self.AddProg('coordconv')
      coordconv.Run()
      if coordconv.out.Get('model',filetype='frac'):
        at = self.inp.Get('model',filetype=oldtyp,has_atomtypes=True).GetAtomType() if self.inp.Get('model',filetype=oldtyp,has_atomtypes=True) else 'S'
        resfname=ins.GetFileName('ins')[:-4]+'.res'
        with open(resfname,'w') as h:
          with open(ins.GetFileName('ins')) as f:
            for fline in f:
              if fline.startswith('HKLF'):
                with open(coordconv.out.Get('model',filetype='frac').GetFileName('frac')) as g:
                  for ig,gline in enumerate(g):
                    if gline:
                      h.write(at+format(ig,'02d')+' 1 '+gline[6:36]+gline[46:49]+'  0.2\n')
              h.write(fline)
        out=self.out.Add(self.inp.Get('model',filetype=oldtyp))
        self.inp.AddFileToChild(out,resfname,'res')
  RegisterConv(setup, ConvertPDB2RES, ['pdb'], ['res'])

  def ConvertSequence(self,oldtyp,newtyp):
    seq=self.inp.Get(filetype=oldtyp)
    if not seq.seqstr:
      seq.GetSequenceString()
    self.SetRunDir()
    newfile=os.path.join(self.rundir,os.path.basename(seq.GetFileName(oldtyp))+'.'+newtyp)
    with open(newfile,'w') as g:
      for i in range(0,len(seq.seqstr),80):
        g.write(seq.seqstr[i:i+80]+'\n')
    seq.AddFile(newfile,newtyp)
    self.out.Add(seq)
  RegisterConv(setup, ConvertSequence, ['fasta','pir'], ['seq'])


#  def RunPreprocess(self,*args,**kwargs):
#    if self.GetParam('inputformat') and self.GetParam('outputformat'):
#      from_to=(self.GetParam('inputformat'),self.GetParam('outputformat'))
#      self.setup[from_to](self,*from_to)
#      process.RunPreprocess(self,*args,**kwargs)
