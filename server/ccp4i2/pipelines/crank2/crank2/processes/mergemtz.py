#!/usr/bin/python
import os,sys
from process import process
from program import program
import common


class mergemtz(process):
  name="checking mtz files and merging them"
  short_name="merging mtz files"
  supported_progs=["sftools","cad"]
  supported_params={}
  supported_params['cons_check_suff'] = common.parameter( desc='no merge needed for a single input mtz with dname/xname consistent', typ=bool )
  no_reporting=True

  def Init(self):
    self.xdrename, self.xdrename_multiple = {}, {}

  def RunBody(self,*args,**kwargs):
    # switch between CAD and SFtools; SFtools does not support xname/dname setting!
    use_cad=True
    consistent=kwargs.pop('consistent',None)
    consistency_check_sufficient=kwargs.pop('cons_check_suff',self.GetParam('cons_check_suff'))
    keep_nodata=kwargs.pop('keep_nodata',None)
    prog=kwargs.pop('prog',program(None))
    simulate=kwargs.pop('simulate',False)
    filelist=program.GetFileList(prog,self.inp.GetAll())
    self.xdname_old, self.xdname_new = [], []
    if consistent is None:
      consistent = self.CheckMTZ(filelist)
    if consistent and consistency_check_sufficient:
      self.out.AddNew('datafile', filelist[0], 'mtz')
      return
    outmtz=kwargs.pop('outmtz',None)
    if not outmtz:
      if prog.nick!='program':
        outmtz=os.path.join(prog.rundir,prog.nick+'_inpmerge.mtz')
      else:
        outmtz=os.path.join(self.rundir,'inpmerge.mtz')
    if not simulate:
      if not use_cad:
        self.lbl,self.lbl_uniq=program.GetUniqueLabels(prog,objects,filelist)
        sft = self.AddProg('sftools')
        sft.runname=sft.name+'_merge'
        filelist.reverse()
        for f in filelist:
          for l in self.lbl[f]:
            sft.SetKey( 'read', ('"'+f+'"', 'col', "\"{0}\"".format(l)) )
          for l_ren in [l  for l in zip(self.lbl[f],self.lbl_uniq[f])  if l[0]!=l[1]]:
            sft.SetKey('set label col',l_ren[0])
            sft.SetKey(l_ren[1])
        if not keep_nodata:
          sft.SetKey('purge', 'nodata')
        sft.SetKey('write',outmtz)
        sft.SetKey('Y')
        sft.SetKey('quit\nY')
        sft.Run()
        # assuming outmtz was passed as full path here
        sft.out.AddNew('datafile', outmtz, 'mtz')
        self.programs.remove(sft)
      else:
        cad = self.AddProg('cad')
        cad.runname+='_'+self.runname
        cad.outfilename['mtz']=outmtz
        if not keep_nodata:
          cad.SetKey('VALM', ('NaN', 'NOOUTPUT'))
        cad.Run()
        self.programs.remove(cad)
        self.lbl,self.lbl_uniq=cad.lbl,cad.lbl_uniq


  def CheckMTZ(self,filelist):
    """Check the consistency between mtz dname/xname and object dname/xname
       MTZ d/x-names are set to object d/x-names unless object d/x-name is default_unknown
    """
    consistent=True
    sft = self.AddProg('sftools')
    sft_run=None
    for i,fname in enumerate(filelist):
      self.xdrename[fname], self.xdrename_multiple[fname] = {}, {}
      for o in self.inp.GetAll():
        if o.GetFileName('mtz')==fname: # and (o.xname!=o.default_unknown or o.dname!=o.default_unknown):
          if sft_run!=i:
            sft.ClearAnyParams()
            sft.SetKey('read', '"'+fname+'"')
            sft.SetKey('Y')  # for the rare case of interactive question, typically irrelevant about xplor free ref...
            sft.SetKey('list id')
            sft.SetKey('quit\nY')
            sft.runname=sft.name+'_check'+str(i)
            sft.Run()
            sft_run=i
            if sft.GetStat('illegal_cell',accept_none=True):
              common.Error('Illegal cell parameters in the mtz. Have a look at {0} for more details.'.format(sft.GetLogFileName()))
          for l in [o.GetLabel(c) for c in o.col_list  if o.GetLabel(c)]:
            try:
              xname=sft.GetStat('xname',l)
            except Exception as e:
              if hasattr(sys,'exc_clear'): sys.exc_clear()
              common.Error('Error reading mtz file {0} to determine crystal name for label {1}, check the input data'.format(fname,l))
            dname=sft.GetStat('dname',l)
            if ((xname,dname),(o.xname,o.dname)) not in self.xdrename[fname].items():
              if (xname,dname) in self.xdrename[fname]:
                self.xdrename_multiple[fname].setdefault((xname,dname),[self.xdrename[fname][(xname,dname)],]).append((o.xname,o.dname))
              self.xdrename[fname][(xname,dname)] = (o.xname,o.dname)
            if xname!=o.xname or dname!=o.dname:
              consistent=False
    self.programs.remove(sft)
    return consistent
