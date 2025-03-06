from __future__ import with_statement
from future.utils import raise_
import os,sys,copy,shutil
import common,manager,parse


def CallCrankFromCCP4i2(ccp4i2crank, xmlfile=None, inpfile=None, defaults=False, rvapi_style=None):
  from core import CCP4PluginScript, CCP4ErrorHandling
  # we need to return to the original cwd when leaving otherwise i2 gets confused
  cwd_saved = os.getcwd()
  os.chdir(ccp4i2crank.workDirectory)
  parser=parse.parse()
  parser.pars_arg.dirout = ccp4i2crank.workDirectory
  if xmlfile:
    f = open(xmlfile,'r')
    parser.pars_arg.xmlin = f
  elif inpfile:
    f = open(inpfile,'r')
    parser.pars_arg.keyin = f
  else:
    os.chdir(cwd_saved)
    common.Error('Crank process could not be initialized.')
  global crank_logfile, stdout_save
  with open(os.path.join(ccp4i2crank.workDirectory,'log.txt') ,'w', 1) as crank_logfile:
    stdout_save = sys.stdout
    sys.stdout = crank_logfile
    error=0
    try:
      crank = parser.ParseAndRun(ccp4i2crank,defaults,rvapi_style=rvapi_style)
    except Exception as e:
      error=e
    else:
      crank.ccp4i2job = ccp4i2crank
    sys.stdout = stdout_save
    f.close()
    if error:
      os.chdir(cwd_saved)
      # simple raise used to lose the trace...  raise_ is needed to make a python2/3 compatible specific error raise with trace
      #raise_(error,None,sys.exc_info()[2])
      try: #python2
        raise
      except RuntimeError: #python3
        raise error
  # register output objects from the last step... (or from previous step if not present in last)
  if not defaults:
    if hasattr(crank.processes[-1],'ccp4i2job'):
      last_job = crank.processes[-1].ccp4i2job
      for outd_name in last_job.container.outputData._dataOrder:
        RegisterSubOutputAsMain(ccp4i2crank,crank,last_job,outd_name)
      if not ccp4i2crank.container.outputData.XYZOUT_SUBSTR.isSet() and len(crank.processes)>=2 and \
         hasattr(crank.processes[-2].ccp4i2job.container.outputData,'XYZOUT_SUBSTR'):
        RegisterSubOutputAsMain(ccp4i2crank,crank,crank.processes[-2].ccp4i2job,'XYZOUT_SUBSTR')
      if not ccp4i2crank.container.outputData.FPHOUT_DIFFANOM.isSet() and len(crank.processes)>=2 and \
         hasattr(crank.processes[-2].ccp4i2job.container.outputData,'FPHOUT_DIFFANOM'):
        RegisterSubOutputAsMain(ccp4i2crank,crank,crank.processes[-2].ccp4i2job,'FPHOUT_DIFFANOM')
      handp=crank.GetProcess('handdet') or crank.GetProcess('phdmmb')
      if not ccp4i2crank.container.outputData.F_SIGF_OUT.isSet() and handp and \
         hasattr(handp.ccp4i2job.container.outputData,'F_SIGF_OUT'):
        RegisterSubOutputAsMain(ccp4i2crank,crank,handp.ccp4i2job,'F_SIGF_OUT')
      for i in range(4):
        label = 'F_SIGFanom_OUT' if not i  else 'F_SIGFanom_OUT'+str(i+1)
        if not getattr(ccp4i2crank.container.outputData,label).isSet() and handp and \
           hasattr(handp.ccp4i2job.container.outputData, label):
          RegisterSubOutputAsMain(ccp4i2crank,crank,handp.ccp4i2job, label)
      free = crank.processes[0].out.Get('exclude',typ='freeR')
      if free and hasattr(ccp4i2crank.container.outputData,'FREEROUT'):
        filepath=OutFilesDirMatch(free,crank,filetype='mtz')
        getattr(ccp4i2crank.container.outputData,'FREEROUT').setOutputPath(projectId=ccp4i2crank.projectId(), relPath=ccp4i2crank.relPath())
        error = ccp4i2crank.splitHklout(['FREEROUT',], [free.GetLabel('free'),], infile=filepath)
    #crank.ccp4i2.reportStatus(CCP4PluginScript.CPluginScript.SUCCEEDED)
  #else:
  #  shutil.rmtree(ccp4i2crank.workDirectory)
  os.chdir(cwd_saved)
  return crank

def RegisterSubOutputAsMain(i2crank,crank,i2subjob,outd_name):
  if outd_name in i2crank.container.outputData._dataOrder:
    i2_subjob_obj = getattr(i2subjob.container.outputData, outd_name)
    setattr( i2crank.container.outputData, outd_name, i2_subjob_obj )
    if outd_name not in('PERFORMANCE',):
      filepath=OutFilesDirMatch(str(i2_subjob_obj),crank)
      if filepath:
        from core import CCP4ErrorHandling
        try:
          getattr(i2crank.container.outputData, outd_name).setFullPath(filepath)
          getattr(i2crank.container.outputData, outd_name).annotation.set( i2_subjob_obj.annotation )
        except (CCP4ErrorHandling.CException, AttributeError):  # if eg copied as string
          getattr(i2crank.container.outputData, outd_name).set(filepath)

def RegisterProcessToCCP4i2(ccp4i2crank, process):
  from core import CCP4PluginScript, CCP4ErrorHandling, CCP4XtalData
  sys.stdout = stdout_save
  subjob = ccp4i2crank.makePluginObject(pluginName='crank2_'+process.nick, pluginTitle=process.short_name[0].upper()+process.short_name[1:])
  if subjob is None:
    subjob = ccp4i2crank.makePluginObject(pluginName=process.short_name.capitalize(), dummy=True)
  process.SetRunDir(subjob.workDirectory)
  process.ccp4i2job = subjob
  # copy the i2 crank2 input into i2 subjob input containers
  for ic in ('inputData','controlParameters'):
    if hasattr(subjob.container,ic): 
      i2cont, subcont = getattr(ccp4i2crank.container,ic), getattr(subjob.container,ic)
      if not getattr(subjob.container,ic):
        for dc in i2cont._dataOrder:
          subcont.addObject( getattr(i2cont,dc), reparent=True )
      else:
        subcont.copyData( otherContainer=i2cont )
  # define derived (not inputted by user) input objects with i2 (for subprocesses rerunning)
  inp_obj=process.inp.Get(filetype='pdb',typ='substr')
  if inp_obj:
    filepath=OutFilesDirMatch(inp_obj,process,filetype='pdb')
    subjob.container.inputData.XYZIN_SUB.setFullPath(filepath)
    if process.nick=='phdmmb' and len(filepath)>4 and inp_obj.GetFileName('res'):
      subjob.container.inputData.XYZIN_SUB_RES.set(inp_obj.GetFileName('res'))
  inp_obj=process.inp.Get(filetype='pdb',typ=('partial+substr','partial'))
  if inp_obj:
    filepath=OutFilesDirMatch(inp_obj,process)
    subjob.container.inputData.XYZIN.setFullPath(filepath)
  inp_obj=process.inp.Get('mapcoef',filetype='mtz',typ=('best','combined'),inp_cont=process.inp.Get('fsigf',typ='average'))
  if inp_obj:
    filepath,filepath_out = inp_obj.GetFileName('mtz'),'FPHIN_HL.mtz'
    i2obj = subjob.container.inputData.FPHIN_HL
    if inp_obj.GetLabel('hla'):
      i2obj.contentFlag.set(CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL)
      status = subjob.splitMtz(filepath, [[filepath_out, \
        inp_obj.GetLabel('hla')+','+inp_obj.GetLabel('hlb')+','+inp_obj.GetLabel('hlc')+','+inp_obj.GetLabel('hld'),i2obj.columnNames(True)],])
    else:
      i2obj.contentFlag.set(CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM)
      status = subjob.splitMtz(filepath, [[filepath_out, inp_obj.GetLabel('ph')+','+inp_obj.GetLabel('fom'),i2obj.columnNames(True)],])
    if status != CCP4PluginScript.CPluginScript.SUCCEEDED:
      subjob.reportStatus(CCP4PluginScript.CPluginScript.FAILED)
    i2obj.setFullPath(filepath_out)
    subjob.container.inputData.INPUT_PHASES.set(True)
  # this is useful if sca/hkl was inputted and we will want to rerun mbref, ref etc
  inp_obj=process.inp.Get('fsigf',filetype='mtz',typ='average',col='sigf',custom=['native'],try_convert=False)
  if not inp_obj:
    inp_obj=process.inp.Get('fsigf',filetype='mtz',typ='average',col='sigf',xname='native',try_convert=False)
  if inp_obj:
    subjob.container.inputData.SAVED_FAVFILE_NATIVE = inp_obj.GetFileName('mtz')
    subjob.container.inputData.SAVED_FAVER_NATIVE = inp_obj.GetLabel('f')
    subjob.container.inputData.SAVED_SIGFAVER_NATIVE = inp_obj.GetLabel('sigf')
  # previous steps incurred parameters
  if process.nick=='refatompick' and process.IsParam('num_iter') and subjob.container.controlParameters.REFATOMPICK_NUM_ITER and \
     int(process.GetParam('num_iter')) != int(subjob.container.controlParameters.REFATOMPICK_NUM_ITER):
    subjob.container.controlParameters.REFATOMPICK_NUM_ITER.set(int(process.GetParam('num_iter')))
  if process.nick=='phdmmb' and process.IsParam('thorough_build') and subjob.container.controlParameters.PHDMMB_THOROUGH_BUILD and \
     bool(process.GetParam('thorough_build')) != bool(subjob.container.controlParameters.PHDMMB_THOROUGH_BUILD):
    subjob.container.controlParameters.PHDMMB_THOROUGH_BUILD.set(process.GetParam('thorough_build'))
  sys.stdout = crank_logfile
  # "dummy" update of the program.xml file (not needed if the i2 rvapi converter is used)
  #if not ccp4i2crank.rvapi_converter:
  shutil.copy(os.path.join(subjob.workDirectory,'..',process.nick+'.xml'), os.path.join(subjob.workDirectory,'..','program.xml'))
  shutil.copy(os.path.join(subjob.workDirectory,'..',process.nick+'.xml'), os.path.join(subjob.workDirectory,'program.xml'))

# annotations for the ccp4i2 output containers (does not seem possible to assign in xml?)
annotations = { 'FPHOUT_HAND2':    '"Best" density - other hand (not chosen)', \
                'FPHOUT_HL_HAND2': '"Best" phases - other hand (not chosen)', \
                'FPHOUT_2FOFC':    '2mFo-DFc weighted density map coef.', \
                'FPHOUT_DIFF':     'Difference map coef.', \
                'FPHOUT_DIFFANOM': 'Difference anomalous diff. map coef.', \
                'FPHOUT_HL':       '"Best" phases', \
                'FPHOUT':          '"Best" electron density map coefficients', \
                'XYZOUT_HAND2':    'Model coord. - other hand (not chosen)', \
                'XYZOUT_SUBSTR':   'Anomalous substructure coordinates', \
                'XYZOUT':          'Model coordinates', \
                'F_SIGFanom_OUT':  'Anom. data  (changed polar spacegroup)', \
                'F_SIGFanom_OUT2':  'Anom. data 2 (changed polar spacegroup)', \
                'F_SIGFanom_OUT3':  'Anom. data 3 (changed polar spacegroup)', \
                'F_SIGFanom_OUT4':  'Anom. data 4 (changed polar spacegroup)', \
                'F_SIGF_OUT':  'Native data  (changed polar spacegroup)', \
              }

def RegisterOutputToCCP4i2(process,error,nosuccess=False):
  from core import CCP4PluginScript, CCP4ErrorHandling, CCP4XtalData
  sys.stdout = stdout_save
  i2job=process.ccp4i2job
  if error:
    i2job.reportStatus(CCP4PluginScript.CPluginScript.FAILED)
    raise_(error,None,sys.exc_info()[2])
  else:
    for outd_name in i2job.container.outputData._dataOrder:
     if not hasattr(i2job,'out_params') or outd_name in i2job.out_params:
      outd_obj = getattr(i2job.container.outputData, outd_name)
      if outd_name.startswith('XYZOUT'):
        ft='pdb'
        if outd_name=='XYZOUT':
          out_obj=process.out.Get(filetype=ft,typ=('partial+substr','partial',))
        elif outd_name=='XYZOUT_SUBSTR':
          out_obj=process.out.Get(filetype=ft,typ='substr')
        elif outd_name=='XYZOUT_SUB_RES':
          out_obj,ft=process.out.Get(filetype='res',typ='substr'),'res'
        elif outd_name=='XYZOUT_HAND2':
          out_obj=process.other_phas.out.Get(filetype=ft)
        if out_obj is not None:
          filepath=OutFilesDirMatch(out_obj,process,filetype=ft)
          s=outd_obj.set(filepath)  if outd_name=='XYZOUT_SUB_RES'  else outd_obj.setFullPath(filepath)
      elif outd_name.startswith('FPHOUT'):
        if outd_name=='FPHOUT_HAND2':
          out_obj=process.other_phas.out.Get('mapcoef',filetype='mtz',typ=('best','combined'),inp_cont=process.inp.Get('fsigf',filetype='mtz'))
        elif outd_name=='FPHOUT_2FOFC':
          out_obj=process.out.Get('mapcoef',filetype='mtz',typ='weighted')
        elif outd_name=='FPHOUT_DIFFANOM':
          out_obj=process.out.Get('mapcoef',filetype='mtz',typ='anom-diff')
        elif outd_name=='FPHOUT_DIFF':
          out_obj=process.out.Get('mapcoef',filetype='mtz',typ='diff')
        else:
          out_obj=process.out.Get('mapcoef',filetype='mtz',typ=('best','combined'),inp_cont=process.inp.Get('fsigf',typ='average'))
        if out_obj is not None:
          filepath=OutFilesDirMatch(out_obj,process,filetype='mtz')
          outd_obj.setOutputPath(projectId=i2job.projectId(), relPath=i2job.relPath())
          if outd_name.startswith('FPHOUT_HL'):
            if out_obj.GetLabel('hla'):
              outd_obj.contentFlag.set(CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL)
              error = i2job.splitHklout([outd_name,], [out_obj.GetLabel('hla')+','+out_obj.GetLabel('hlb')+','+out_obj.GetLabel('hlc')+','+out_obj.GetLabel('hld'),], infile=filepath)
            else:
              outd_obj.contentFlag.set(CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM)
              error = i2job.splitHklout([outd_name,], [out_obj.GetLabel('ph')+','+out_obj.GetLabel('fom'),], infile=filepath)
          else:
            error = i2job.splitHklout([outd_name,], [out_obj.GetLabel('f')+','+out_obj.GetLabel('ph'),], infile=filepath)
          if error.maxSeverity()>CCP4ErrorHandling.SEVERITY_WARNING:
            i2job.reportStatus(CCP4PluginScript.CPluginScript.FAILED)
      elif outd_name.startswith('F_SIGF'):
        out_obj=None
        if outd_name.startswith('F_SIGFanom_OUT') and (outd_name[-1]=='T' or int(outd_name[-1])<=len(process.out.GetDataList())):
          ind = 0 if outd_name[-1]=='T' else int(outd_name[-1])-1
          out_obj=process.out.Get('fsigf',filetype='mtz',typ='plus',custom='otherhand',dname=process.out.GetDataList()[ind])
          out_obj2=process.out.Get('fsigf',filetype='mtz',typ='minus',custom='otherhand',dname=process.out.GetDataList()[ind])
          if len(process.out.GetDataList())>1 and out_obj and process.out.GetDataList()[ind]!=out_obj.default_unknown:
            annotations[outd_name] = annotations[outd_name].replace(' '+str(ind+1) if ind else '  ',' '+process.out.GetDataList()[ind])
        if outd_name=='F_SIGF_OUT':
          out_obj=process.out.Get('fsigf',filetype='mtz',typ='average',custom='otherhand',is_native=True)
        if out_obj is not None:
          fi = 'f' if out_obj.GetLabel('f') else 'i'
          if outd_name=='F_SIGF_OUT':
            i2flag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN if fi=='f' else CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN
          else:
            i2flag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_FPAIR if fi=='f' else CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR
          filepath=OutFilesDirMatch(out_obj,process,filetype='mtz')
          outd_obj.setOutputPath(projectId=i2job.projectId(), relPath=i2job.relPath())
          outd_obj.contentFlag.set(i2flag)
          if outd_name=='F_SIGF_OUT':
            error = i2job.splitHklout([outd_name,], [out_obj.GetLabel(fi)+','+out_obj.GetLabel('sig'+fi)], infile=filepath)
          else:
            error = i2job.splitHklout([outd_name,], [out_obj.GetLabel(fi)+','+out_obj.GetLabel('sig'+fi)+','+out_obj2.GetLabel(fi)+','+out_obj2.GetLabel('sig'+fi),], infile=filepath)
      if outd_name in annotations:
        outd_obj.annotation.set(annotations[outd_name])
        if outd_name=='FPHOUT_DIFFANOM' and process.GetParam('target')=='SAD':
          outd_obj.annotation.set("Anomalous gradient LLG (difference) map")
    if 'PERFORMANCE' in i2job.container.outputData._dataOrder:
      for perf_name in i2job.container.outputData.PERFORMANCE.CONTENTS_ORDER:
        if hasattr(i2job,'perform') and perf_name in i2job.perform:
          if perf_name=='RFactor' and hasattr(process,'report_R'):
            i2job.container.outputData.PERFORMANCE.RFactor.set(process.report_R)
          elif perf_name=='RFree' and hasattr(process,'report_Rfree'):
            i2job.container.outputData.PERFORMANCE.RFree.set(process.report_Rfree)
          elif perf_name=='FOM' and hasattr(process,'report_fom'):
            i2job.container.outputData.PERFORMANCE.FOM.set(process.report_fom)
          elif perf_name=='CFOM' and hasattr(process,'score') and process.score>=0.0:
            i2job.container.outputData.PERFORMANCE.CFOM.set(process.score)
          elif perf_name in ('Hand1Score','Hand2Score') and hasattr(process,'score'):
            i2job.container.outputData.PERFORMANCE.Hand1Score.set(process.score[0])
            i2job.container.outputData.PERFORMANCE.Hand2Score.set(process.score[1])
          elif perf_name=='CC' and hasattr(process,'best_cc'):
            i2job.container.outputData.PERFORMANCE.CC.set(process.best_cc)
    i2job.reportStatus(CCP4PluginScript.CPluginScript.SUCCEEDED)
    sys.stdout = crank_logfile
  if nosuccess:
    #i2job.reportStatus(CCP4PluginScript.CPluginScript.UNSATISFACTORY)
    raise

def OutFilesDirMatch(out_obj,process,filetype=None):
  # makes sure the output file in the output directory expected by ccp4i2
  try:
    f=out_obj.GetFileName(filetype=filetype)
  except AttributeError:
    f=out_obj
  if not f:
    return None
  prefix, f2 = '', f
  while os.path.isfile(f2):
    f2 = os.path.join(process.ccp4i2job.workDirectory, prefix+os.path.basename(f))
    prefix=prefix+'n_'
  if os.path.isfile(f):
    shutil.copy(f,f2)
  return f2
