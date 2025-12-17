import copy
import os
import re
import subprocess
import sys
import threading
import time
from collections.abc import Iterable
from distutils import spawn
from queue import Queue
from xml.etree import ElementTree as ET

import common
import data
import inout


class program(object):
  """base program class"""
  binary=None
  name=None
  stat={}
  references=()
  # modifier (eg '-', '--') to be added to the front of each program argument
  modif=''
  # True if the program uses spaces between argument and its value, False otherwise (eg -s 1 vs -s1)
  arg_divided=True
  # prefix added to (the front of) output labels of the program
  labelout_prefix=''
  # log filename extension
  log_suffix='.log'
  # if True then new mtz file is created on "merge" even if we only have a single mtz file
  # the difference is that only specified labels will be included (useful eg for clipper progs
  # which seem to work incorrectly if a label is to be created that already exists in input mtz)
  always_merge=False
  # merging can be skipped for programs that support multiple input mtz
  never_merge=False
  # merging behaviour for reflections with no data in the merged file: True=keeping, False=purging
  keep_nodata=False
  # run directory is assigned immediately as the program instance is added to its process parent.  Use with care
  always_rundir=False
  # ccp4 parsing - if yes then setting/comparison is done on the first 4 capitalized chars keys and capitalized args
  ccp4_parsing=False
  # ccp4 keys (only) parsing - if yes then setting/comparison is done on the first 4 capitalized chars keys (args unchanged)
  ccp4_keys_parsing=False
  # if True then continuous output parsing (using Interact_output() member function) is enabled
  # should be only used when needed as it comes with some performance drawbacks
  interact_output = False
  # if no_rewrite is True then output files assigned using (Set/Add)FileToChild are checked 
  # for existence and adjusted to a non-existing filename
  no_rewrite=False
  # list of supported (virtual) parameters (disabled for programs as of now)
  supported_params = {}
  supported_params['binary'] = common.parameter( desc='Specifies the binary to be used for the program', typ=str )
  supported_params['no_rewrite'] = common.parameter( desc='Output files names will be adjusted not to rewrite existing files', typ=bool )
  supported_params['ignore_labelout_prefix'] = common.parameter( desc='Default program prefix for output mtz labels will be ignored', typ=bool )
  supported_params['allow_basiclib'] = common.parameter( desc='Whether the quick but rudimentary basiclib is allowed', typ=str )

  def __init__(self, parent_process, xmlelem=None, inpline=None, dummy=False):
    self.nick=self.__class__.__name__
    self.process = parent_process
    self.rundir=None
    if self.name:
      self.runname=re.sub(r'\W+', '', self.name)
    else:
      self.runname="prog"
    self.inp = inout.input_output(is_output=False,parent=self)
    self.out = inout.input_output(is_output=True,parent=self)
    # dictionary of 'virtual' parameters
    self.param = {}
    # dictionary of (actual) program keywords (passed to the program as input)
    self.key = {}
    # dictionary of (actual) program arguments (passed to the program from the command line)
    self.arg = {}
    # list of all keys and args - the args and keys will be used in this order in the run of the program
    self.arg_list = []
    self.key_list = []
    # program specific environment may be set here
    self.env={}
    self.init_xmlelem = None
    if xmlelem is not None:
      self.XMLElemInit(xmlelem,dummy)
    elif inpline is not None:
      self.InputElemInit(inpline,dummy)
    if not dummy:
      self.Init()

  def Init(self):
    """Custom initialization of the program can be placed here"""
    pass

  def CheckBinary(self, silent=False, try_run=False):
    """Checks existence of binary of the program.
       If 'silent' then False is returned on non-existence otherwise error message is printed (False by default)
       If 'try_run' then an attempt is made to run the binary with no args and errors reported if not 'silent'
    """
    if not spawn.find_executable(self.binary):
      if silent:
        return False
      else:
        common.Error('Binary {0} for program {1} does not exist.'.format(self.binary, self.name))
    if try_run:
      self.SetRunDir('checkbinary')
      self.runname='{0}_version'.format(self.name)
      bin_issue=self.name+' binary problem: '
      try:
        self.ExternalRun( [self.binary,], [] )
      except Exception as e:
        bin_issue += str(e)
        if 'Permission denied' in str(e):  bin_issue+='. Please make the binary executable.'
        if hasattr(sys,'exc_clear'): sys.exc_clear()
        if silent:
          return False
        else:
          common.Error(bin_issue)
    return True

  @classmethod
  def from_xml(cls, xml, parent_process):
    """Create program instance from XML file or XML tree object"""
    if common.is_file(xml) or common.is_string(xml):
      xml=common.ReadXML(xml)
    inst = cls( xml, parent_process )
    return inst

  @classmethod
  def from_name(cls, prognick, parent, xmlelem=None, inpline=None, dummy=False):
    """Create specific (derived) program instance from the (nick)name of the program"""
    try:
      # import the specific program/process
      setattr(sys.modules[__name__], 'programs', __import__('programs.'+prognick))
    except (AttributeError,ImportError):
      common.Error('Program {0} not supported (check the spelling)'.format(prognick))
    try:
      # create an instance of the program
      inst = getattr(sys.modules['programs.'+prognick], prognick)(parent,xmlelem,inpline,dummy=dummy)
    except (AttributeError,KeyError):
      common.Error('Error when creating program {0} instance.'.format(prognick))
    else:
      return inst

  def XMLElemInit(self, xelem, dummy=False):
    """Read input, parameters definitions from XML into Crank object"""
    tags=["inp","out","param","key","arg"]
    inpnum,outnum,stri,stro=0,0,'',''
    for xchild in xelem:
      if xchild.tag in tags:
        if xchild.tag=="inp":
          if not hasattr(self,'inp'+stri):
            setattr(self,'inp'+stri, inout.input_output(is_output=False,parent=self) )
          getattr(self,'inp'+stri).XMLElemInit(xchild,dummy=dummy)
          inpnum+=1
          stri=str(inpnum+1)
        elif xchild.tag=="out":
          if not hasattr(self,'out'+stro):
            setattr(self,'out'+stro, inout.input_output(is_output=True,parent=self) )
          getattr(self,'out'+stro).XMLElemInit(xchild,dummy=dummy)
          outnum+=1
          stro=str(outnum+1)
        elif xchild.tag in ("param","key","arg"):
          for xparam in list(xchild):
            self.SetParam( xparam.tag, common.AutoConvert(xparam.text), is_key=xchild.tag=="key", is_arg=xchild.tag=="arg" )
      else:
        common.Error('Wrong tag when parsing XML for program {0}: {1}'.format(self.name,xchild.tag))
    self.init_xmlelem = xelem

  def InputElemInit(self, input_line, dummy=False):
    """Read input, parameters definitions from the input line into the Crank object
       Consumes those tokens that were used, stops at first token that cannot be used
    """
    while input_line:
      token=input_line[0]
      input_line.remove(token)
      if token in [ sc.__name__ for sc in data.data_container.__subclasses__() ]:
        self.inp.AddCopy( data.data_container.from_name(token, inpline=input_line), propagate=not dummy )
      # this is a workaround for ccp4i1 - an alternative way to pass keywords to the program
      elif token.startswith('keywords=') or token.startswith('keyword='):
        trash,sep,keyval = token.partition('=')
        key,sep,val = keyval.partition(' ')
        key,val = key.strip(),val.strip()
        if val=='':
          val=True
        if key and token.startswith('keyword'):
          self.AddToKey( key, common.AutoConvert(val) )
      elif '::' in token:
        key,sep,val = token.partition('::')
        self.SetParam( key, common.AutoConvert(val) )
      elif ';' in token:
        key,sep,val = token.partition(';')
        self.AddToArg( key, common.AutoConvert(val) )
      elif ':' in token:
        key,sep,val = token.partition(':')
        self.AddToKey( key, common.AutoConvert(val) )
      else:
        input_line.insert(0,token)
        break

  def Data2XML(self, ET_element):
    """Convert the actual program object into XML"""
    xroot = ET.SubElement(ET_element,'program')
    xroot.text=self.__class__.__name__
    for t in ('param','key','arg'):
      if getattr(self,t):
        xt = ET.SubElement(xroot,t)
        if hasattr(self,t+'_list'):
          num={}
          for key in getattr(self,t):
            num[key]=0
          for key in getattr(self,t+'_list'):
            ET.SubElement(xt, key).text = str( getattr(self,t)[key].value[num[key]] )
            num[key]+=1
        else:
          for key,par in self.param.items():
            ET.SubElement(xt, key).text = str(par.value[0])
    hand2io = [getattr(self,io+'2') for io in ['inp','out'] if hasattr(self,io+'2')]
    for obj in [self.inp, self.out] + hand2io:
      obj.Data2XML(xroot)

  def GetFileList(self, objects, filetype='mtz'):
    # returns list of (unique) names of files of given type from the inputted objects
    filelist = []
    for o in objects:
      f=o.GetFile(filetype)
      if f  and  f.name is not None  and  f.name not in filelist:
        filelist.append(f.name)
    return filelist

  def GetUniqueLabels(self, objects, filelist=None, allow_multiple=False):
    # returns unique labels for inputted objects
    if filelist is None:
      filelist=self.GetFileList(objects)
    # get the list of all labels for each mtz file
    lbl, lbl_uniq = dict([(f,[]) for f in filelist]), {}
    for f in filelist:
      # first creating a dictionary with filenames as keys and list of labels as values
      for o in objects:
        if o.GetFile('mtz') and o.GetFile('mtz').name==f:
          lbl[f].extend( [o.GetLabel(c) for c in o.col_list  if o.GetLabel(c)] )
    all_labels=[l for f in filelist for l in lbl[f]]
    for f in filelist:
      # now assigning unique labels
      if not allow_multiple:
        lbl[f] = [l  for i,l in enumerate(lbl[f])  if lbl[f].index(l)==i]
      prev_uniq = [l  for f2 in filelist  if f2 in lbl_uniq  for l in lbl_uniq[f2]]
      lbl_uniq[f] = []
      for l in lbl[f]:
        suf=0
        l_new=l
        while l_new in prev_uniq or l_new in lbl_uniq[f] or (suf>0 and l_new in all_labels):
          suf+=1
          l_new=l+'_'+str(suf)
        lbl_uniq[f].append(l_new)
    return lbl,lbl_uniq

  def MergeMTZ(self, *objects, **options):
    """Determines whether the mtz file(s) of inputted objects require adjustment or merging.
       Consistency of dname/xname in mtz files and parent object is checked before the actual merge.
       In case of inconsistency and/or multiple input mtz files the adjustment and/or merging is done
       (this behaviour can be modified by 'force_merge' and 'no_check' options)
       If a merging/adjustment is done then the merged/adjusted mtz file is returned.
       If no merging/adjustment done then the input mtz file(s) or None (no input mtz files) is returned.
       The input objects are adjusted to point to the new merged/adjusted mtz file.
       Any label conflicts are resolved by numerical suffix addition.
    """
    force_merge=options.pop('force_merge',False)
    no_check=options.pop('no_check',False)
    keep_nodata=options.pop('keep_nodata',self.keep_nodata)
    filelist=self.GetFileList(objects)
    # determine whether the xname/dname consistency check is needed
    consistent=None
    if no_check or self.never_merge or \
       all(o.xname==o.default_unknown and o.dname==o.default_unknown  for o in objects):
        # consistency check not needed
        consistent=True
    cons_check_suff=False
    self.files_before_merge=[]
    if self.never_merge and not force_merge:
      # merging is never done and multiple files are possible
      return filelist
    elif len(filelist)==1 and not self.always_merge and not force_merge:
      if consistent:
        # no mtz merging/adjusting is necessary
        return filelist[0]
      else:
        # if consistency test passes then no mtz merging/adjusting will be necessary
        cons_check_suff=True
    if len(filelist)>=1:
      # mtz merging/adjusting
      mergemtz_proc=self.process.AddProcess('mergemtz',propagate_inp=False,propagate_out=False)
      for o in objects:
        mergemtz_proc.inp.Add(o)
      mergemtz_proc.Run( consistent=consistent, cons_check_suff=cons_check_suff, \
                         keep_nodata=keep_nodata, prog=self )
      outmtz=mergemtz_proc.out.Get(filetype='mtz').GetFileName('mtz')
      # adjust the mtz files and labels to point to the merged file for all the objects
      if outmtz!=filelist[0] or len(filelist)>1:
        for o in objects:
          f=o.GetFile('mtz')
          if f and f.name is not None and f.name in mergemtz_proc.lbl:
            labels_changed=False
            for i,l in enumerate(mergemtz_proc.lbl[f.name]):
              for k,l2 in o.label.items():
                if l2==l and l!=mergemtz_proc.lbl_uniq[f.name][i]:
                  o.SetLabel(k,mergemtz_proc.lbl_uniq[f.name][i],ignore_prefix=True)
                  labels_changed=True
            self.files_before_merge.append((o,f.name,labels_changed))
            f.name = outmtz
      self.process.processes.remove(mergemtz_proc)
      return outmtz
    else:
      return None

  def BackupAnyPars(self):
    # backups all actual parameters into private tmp arrays
    self.partmp, self.argtmp, self.keytmp = copy.deepcopy(self.param), copy.deepcopy(self.arg), copy.deepcopy(self.key)
    self.arglsttmp, self.keylsttmp = copy.deepcopy(self.arg_list), copy.deepcopy(self.key_list)
    return (copy.deepcopy(self.param), copy.deepcopy(self.arg), copy.deepcopy(self.key), \
           copy.deepcopy(self.arg_list), copy.deepcopy(self.key_list))

  class RestoreParamsError(Exception):
    pass

  def RestoreAnyPars(self,rest_list=None):
    # restores the previously backuped parameters from private tmp arrays
    if rest_list is not None:
      self.param, self.arg, self.key = copy.deepcopy(rest_list[0]), copy.deepcopy(rest_list[1]), copy.deepcopy(rest_list[2])
      self.arg_list, self.key_list = copy.deepcopy(rest_list[3]), copy.deepcopy(rest_list[4])
    else:
      try:
        self.param, self.arg, self.key = copy.deepcopy(self.partmp), copy.deepcopy(self.argtmp), copy.deepcopy(self.keytmp)
        self.arg_list, self.key_list = copy.deepcopy(self.arglsttmp), copy.deepcopy(self.keylsttmp)
      except AttributeError:
        raise self.RestoreParamsError('Restoring backup parameters failed, backup does not exist.')

  def SetRunDir(self,rundir=None,error_on_failure=True,change_cwd=False):
    if rundir is not None:
      self.rundir=rundir
    else:
      rd_obj=self.process
      while self.rundir is None and rd_obj is not None:
        self.rundir=rd_obj.rundir
        if self.rundir is None:
          rd_obj=rd_obj.parent_process
        elif len(self.process.programs+self.process.processes)>1 or rd_obj!=self.process:
          self.rundir=os.path.join(self.rundir,self.nick)
      if self.rundir is None:
        if error_on_failure:
          common.Error('Run directory for program {0} could not be set as it is not set for its parent process {1}'.
              format(self.name,self.process.name))
        else:
          return
    if not os.path.isdir(self.rundir):
      os.mkdir(self.rundir)
    if not os.path.isabs(self.rundir):
      self.rundir=os.path.abspath(self.rundir)
    if change_cwd:
      self.previous_cwd=os.getcwd()
      os.chdir(self.rundir)

  def CatchMTZObjects(self):
    # catch MTZ objects by performing the actual input and param treatment and keeping track of used MTZ
    for o in self.inp.GetAll(try_convert=False):
      o.mtzfile_asked_before[self.name]=o.mtzfile_asked
      o.mtzfile_asked=False
    #suppressing warnings, to prevent them being printed twice as this routine is typically followed by TreatParams() again
    with open(os.devnull, 'w') as fnull:
      sys_stdout = sys.stdout
      sys.stdout = fnull
      self.TreatParams()
      sys.stdout = sys_stdout
    self.TreatInput()
    self.RestoreAnyPars()
    result = [o for o in self.inp.GetAll(try_convert=False) if o.mtzfile_asked]
    for o in self.inp.GetAll(try_convert=False):
      o.mtzfile_asked=o.mtzfile_asked_before[self.name]
    self.caught_mtz_objs=result
    return result

  def Run(self,rundir=None,restore=True,clear_out=True,check_bin=False,lock=None,return_to_orig_filenames=True,restore_env=True):
    """run the program"""
    if check_bin:
     self.CheckBinary()
    self.SetRunDir(rundir,change_cwd=True)
    self.BackupAnyPars()
    self.env_backup=dict(self.env)
    if clear_out:
      self.out.ClearAll(propagate=False)
    # merge and treat parameters and input
    self.MergeMTZ( *self.CatchMTZObjects() )
    self.TreatParams()
    self.TreatInput()
    #define and treat output
    self.DefineOutput()
    self.TreatOutput()
    if self.GetCrankParent() and self.GetCrankParent().logfilehandle:
      self.GetCrankParent(prep=True).Info("{0}Going to run program {1}, script {2}.".format( \
              '  '*self.GetNumParents(), self.name, self.GetScrFileName(relative=self.GetCrankParent())), stdout=False)
    self.CheckParamAccess()
    # create the actual script
    args, scr_lines = self.GetRunLists()
    # run!
    log = self.ExternalRun(args,scr_lines,lock=lock)
    # return to unprocessed parameters for case we will want to rerun
    if restore_env:
      self.env=self.env_backup
    if restore:
      self.RestoreAnyPars()
    # return to original filenames before the merge
    # it would be a bug in case the labels are changed! the original input filenames can be found in files_before_merge.
    for o,fname,label_change in self.files_before_merge:
      if not label_change and return_to_orig_filenames:
        o.GetFile('mtz').name=fname
    # return to the previous working directory
    os.chdir(self.previous_cwd)
    return log

  def GetStat(self, stat, param=[], accept_none=False, param_escape=True, multiple=False):
    """Returns the requested statistics; assumes the regexp/xpath has been specified in 
       program class attribute dictionary stat"""
    result = self.GetStatXML(stat,multiple=multiple)
    if result is not None: 
      return result
    result = self.GetStatGrep(stat,param,param_escape=param_escape,multiple=multiple)
    if result is not None or accept_none or self.stat[stat].accept_none:
      return result
    else:
      common.Error('The requested statistics {0} could not be obtained for program {1}'.format(stat,self.name))

  def GetStatGrep(self, stat, param=[], from_str=None, param_escape=True, multiple=False):
    if from_str is None:
      logfile=self.GetLogFileName()
      if not logfile:
        common.Warning("{0} logfile of program {1} could not be found.".format(logfile,self.name))
        return None
      if os.path.getsize(logfile)>100000000:
        common.Warning("{0} logfile of program {1} too big for grepping.".format(logfile,self.name))
        return None
      if stat not in self.stat or self.stat[stat].regexp is None:
        common.Warning("No definition for grepping stat {0} of program {1}.".format(stat,self.name))
        return None
      f = open(logfile,"r")
      from_str = f.read()
    if param:
      if common.is_string(param) or not isinstance(param, Iterable):
        param = [param,]
      if param_escape:
        param = [re.escape(str(p)) for p in param]
    value=[]
    if common.is_string(self.stat[stat].regexp):
      self.stat[stat].regexp = [self.stat[stat].regexp,]
    for gs in self.stat[stat].regexp:
      value.extend( re.findall(gs.format(*param), from_str) )
    if value:
      if self.stat[stat].multiple or multiple:
        value=[self.stat[stat].ConvertIfAsked(v) for v in value]
      else:
        value=self.stat[stat].ConvertIfAsked(value[-1])
    else:
      if not self.stat[stat].multiple or multiple:
        value=None
      #common.Warning("Grepping stat {0} of program {1} not succesfull.".format(stat,self.name))
    return value

  def GetStatXML(self, stat, multiple=False):
    value=None
    xmlobj=self.out.Get(filetype='xml')
    if not xmlobj or stat not in self.stat or self.stat[stat].xpath is None:
      return None
    else:
      xmlfile=xmlobj.GetFileName('xml')
    try:
      xtree = ET.ElementTree( file=xmlfile )
      value = xtree.findall(self.stat[stat].xpath)
    except Exception:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
    if not value:
      try:
        tag_split = self.stat[stat].xpath.split('/',1)
        elems = xtree.findall(tag_split[1])
      except Exception:
        if hasattr(sys,'exc_clear'): sys.exc_clear()
    if not value:
      if not self.stat[stat].multiple or multiple:
        return None
      return []
    else:
      attr=self.stat[stat].attrib
      if self.stat[stat].multiple or multiple:
        return [self.stat[stat].ConvertIfAsked(v.attrib[attr] if attr else v.text) for v in value]
      else:
        return self.stat[stat].ConvertIfAsked(value[-1].attrib[attr] if attr else value[-1].text)

  def Clean(self, *files):
    if not self.debug:
      for f in files:
        if os.path.isfile(f):
          os.remove(f)

  def TreatInput(self):
    """Translates actual program input into keyworded program input;
       (mtz) file names and labels must retrieved using container's GetFileName(),GetLabel() methods
       in order for automatic merging to work properly
    """
    pass

  def TreatParams(self):
    """Translates 'virtual' parameters into program keywords"""
    if 'binary' in self.param:
      self.binary=self.GetParam('binary')

  def TreatOutput(self):
    """Translates actual program output definition into keyworded program output"""
    pass

  def DefineOutput(self):
    """Defines program output (as a function of input and parameters)"""
    pass

  def CheckOutput(self):
    pass
    # checks whether the defined output files do not exist already
    #for of in [ f  for o in self.out.GetAll() for f in o.GetFiles() ]:
    #  while os.path.isfile(of.name):
      #if os.path.isfile(of.name):
        # if the file exists already then check whether it is used by any other object of this step
        #greatest_parent=self.process
        #while greatest_parent.parent_process is not None:
        #  greatest_parent=greatest_parent.parent_process
        #loop over a list of all objects anywhere in the process/program tree of this step
        #for p in greatest_parent.GetFlatSubTree():
        #  #print 'debug',p.nick
        #  if p!=self:
        #    for o in p.inp.GetAll()+p.out.GetAll():
        #      #print 'debug2',o.nick
        #      if of.name in [f.name for f in o.GetFiles()] and o not in self.out.GetAll():
        #        common.Warning('File {0} already exists'.format(of.name))
        #        break


  def GetFlatSubTree(self,lst=None):
    if lst is None:
      lst=[]
    lst.append(self)
    return lst

  def CheckParamAccess(self):
    """Checks whether all the inputted 'virtual' parameters have been evaluated"""
    for parname,parinst in self.param.items():
      if not parinst.is_key and not parinst.is_arg and not parinst.accessed:
        common.Warning('Virtual parameter {0} has not been used by program {1}.'.format(parname,self.name))

  def AddArgKeyToList(self, par, is_key, is_arg, ins=None):
    # Adding the arg/key to the list of args/keys
    if is_arg:
      if ins is None:
        self.arg_list.append(par)
      else:
        self.arg_list.insert(ins,par)
    elif is_key:
      if ins is None:
        self.key_list.append(par)
      else:
        self.key_list.insert(ins,par)

  def RemoveArgKeyFromList(self, par, is_key, is_arg):
    # Removes the arg/key from the list of args/keys (if set multiple times then removes all the settings)
    if is_arg or is_key:
      if is_arg:
        self.arg_list = [p  for p in self.arg_list  if p!=par]
      else:
        self.key_list = [p  for p in self.key_list  if p!=par]

  def GetParParsingStyle(self, par, is_key, is_arg):
    if (self.ccp4_parsing or self.ccp4_keys_parsing) and (is_key or is_arg):
      return self.GetParCCP4(par, is_key, self.ccp4_parsing)
    return par

  def GetParCCP4(self, par, is_key, full_ccp4_parsing=True):
    # Return the inputted parameter (key/arg) in the CCP4 style - first 4 capitalized characters
    # unmodified key/arg returned if it contains spaces
    if ' ' in par:
      return par
    elif is_key:
      return par[:4].upper()
    elif full_ccp4_parsing:
      return par.upper()
    return par

  def GetParArray(self, is_key=None, is_arg=None):
    # Returns the key/arg/virtpar dictionary based on whether we deal with a key,arg or virtpar
    array=self.param
    if is_key:
      array=self.key
    elif is_arg:
      array=self.arg
    return array

  def GetSupportedParam(self,par):
    # specific virtual parameters for programs disabled as of now except of the general program 
    # subclass ones; we can enable when they are needed.
    if par in program.supported_params:
      self.supported_params[par] = program.supported_params[par]
      return par
    else:
      common.Error('Cannot set parameter {0} to program {1}: (virtual) parameters can be only set for a process.'.format(
                  par,self.name))

  def ConvRelPath(self, value):
    # called by SetParam() to convert path to relative path to self.rundir
    try:
      iter(value)
    except TypeError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
    else:
      if common.is_string(value):
        value = [value,]
      else:
        value = list(value)
      for i,v in enumerate(value):
        if common.is_string(v):
          v=v.strip('"\'')
          if len(os.path.commonprefix([self.rundir,v]))>5 and len(os.path.relpath(v, self.rundir))<len(v):
            value[i] = os.path.relpath(v, self.rundir)
      if len(value)==1:
        value=value[0]
    return value

  def SetParam(self, par, value=True, is_key=None, is_arg=None, append=True, ind=None,  \
               relpath=False, ignore_modif=False, insert=None, inputted=False):
    # Set program parameter 'par' to inputted value
    array=self.GetParArray(is_key, is_arg)
    par=self.GetParParsingStyle(par, is_key, is_arg)
    if not is_key and not is_arg:
      append=False
      self.GetSupportedParam(par)
    if par not in array:
      array[par] = common.parameter()
      append=True
    # conversion to relative path
    if relpath and self.rundir:
      value = self.ConvRelPath(value)
    # set in the parameter dictionary
    try:
      array[par].Set(value,is_key,is_arg,append=append,ind=ind,ignore_modif=ignore_modif)
    except Exception as e:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      common.Error('Could not set \'{0}\' for program {1}: {2}'.format(par,self.name,str(e)))
    # set in the list of keywords/arguments if needed
    if append:
      self.AddArgKeyToList(par,is_key,is_arg,ins=insert)
    elif ind is None:
      self.RemoveArgKeyFromList(par,is_key,is_arg)
      self.AddArgKeyToList(par,is_key,is_arg,ins=insert)

  def UnsetParam(self, par, is_key=None, is_arg=None):
    # unset the parameter.  experimental.
    array=self.GetParArray(is_key, is_arg)
    par=self.GetParParsingStyle(par, is_key, is_arg)
    if not is_key and not is_arg:
      self.GetSupportedParam(par)
    array[par].Unset(is_key,is_arg)
    self.RemoveArgKeyFromList(par, is_key, is_arg)

  def GetParam(self, par, is_key=None, is_arg=None, ind=-1, as_list=False, capitalized=False, allval=False):
    # Returns the value of parameter 'par', None if not defined
    array=self.GetParArray(is_key, is_arg)
    par=self.GetParParsingStyle(par, is_key, is_arg)
    if par in array:
      array[par].accessed=1
      if allval:
        return array[par].GetAll(capitalized=capitalized)
      else:
        return array[par].Get(ind=ind,as_list=as_list,capitalized=capitalized)
    else:
      if as_list or allval:
        return []
      return None

  def AddToParam(self, par, value=True, is_key=None, is_arg=None, ind=-1):
    #Adds inputted value to program parameter 'par' by appending/expending the previous value list
    array=self.GetParArray(is_key, is_arg)
    par=self.GetParParsingStyle(par, is_key, is_arg)
    if par not in array:
      new_value = value
    else:
      prev_value = self.GetParam(par,is_key,is_arg,ind=ind,as_list=True)
      new_value = copy.copy(prev_value)
      try:
        assert not common.is_string(value)
        new_value.extend(value)
      except (AssertionError,TypeError):
        if hasattr(sys,'exc_clear'): sys.exc_clear()
        new_value.append(value)
    self.SetParam(par,new_value,is_key,is_arg,append=False,ind=ind)

  def AddToArg(self, arg, value=True):
    """Adds inputted 'value' to program argument 'arg' by appending/expending the previous value list
       If the argument's value before addition was not list (str,int etc) then it is converted to a list
       If the added 'value' is a scalar then it is appended, if it is a list then extension takes place
       If there was no such argument set before then it is set
       In command line input, multiple argument values are printed after the keyword separated by space
    """
    self.AddToParam(arg, value, is_key=False, is_arg=True)

  def AddToKey(self, key, value=True, ind=-1):
    """Adds inputted 'value' to program argument 'key' by appending/expending the previous value list
       If the argument's value before addition was not list (str,int etc) then it is converted to a list
       If the added 'value' is a scalar then it is appended, if it is a list then extension takes place
       If there was no such keyword set before then it is set
       If there were multiple such keywords set before then 'ind' specifies the index for addition
       Multiple keyword values mean single line with keyword followed by values separated by space
    """
    self.AddToParam(key, value, is_key=True, is_arg=False, ind=ind)

  def SetArg(self, arg, value=True, keep_previous=True, relpath=True, ignore_modif=False, insert_index=None):
    """Set inputted arg as program argument"""
    assert keep_previous in (False,True), \
      'Setting par {0} to {1}: {2} might need to be included in a list?'.format(arg,value,keep_previous)
    self.SetParam(arg,value,is_key=0,is_arg=1,append=keep_previous,relpath=relpath,ignore_modif=ignore_modif,insert=insert_index)

  def GetArg(self, arg, ind=-1, capitalized=False, allval=False):
    """Returns the value of argument, None if no such argument defined
       In case of multiple argument values, 'ind' specifies the index of value returned (-1 by default)
       and if 'allval' is True then a list of all values is returned (False by default)
    """
    return self.GetParam(arg,is_key=0,is_arg=1,capitalized=capitalized,ind=ind,allval=allval)

  def SetKey(self, key, value=True, keep_previous=True, relpath=True, insert_index=None):
    """Set inputted key as program keyword"""
    assert keep_previous in (False,True), \
      'Setting par {0} to {1}: {2} might need to be included in a list?'.format(key,value,keep_previous)
    self.SetParam(key,value,is_key=1,is_arg=0,append=keep_previous,relpath=relpath,insert=insert_index)

  def GetKey(self, key, ind=-1, capitalized=False, allval=False, as_list=False):
    """Returns the value of keyword, None if no such program keyword defined

       Parameters:
       ind - specifies the index of value returned in case of multiple keyword values (default: the last set)
       capitalized - return capitalized value(s) (False by default)
       allval - a list of all values is returned (False by default)
       as_list - return as a list (ie if the last setting was a single value, a list with this value is returned)
    """
    return self.GetParam(key,is_key=1,is_arg=0,ind=ind,capitalized=capitalized,allval=allval,as_list=as_list)

  def SetVirtPar(self, par, value=True):
    """Set inputted key as program virtual parameter"""
    self.SetParam(par,value,is_key=0,is_arg=0)

  def GetVirtPar(self, par):
    """Returns the value of virtual parameters, None if no such parameter defined"""
    return self.GetParam(par,is_key=0,is_arg=0)

  def GetAnyParam(self, par, value=True):
    """Searches all program parameters - virtuals, arguments, keywords in this order
       are returns a list of parameters found, empty if no such parameter 'par' found
    """
    lst=[]
    if self.GetParam(par,is_key=0,is_arg=0) is not None:
      lst.append( self.GetParam(par,is_key=0,is_arg=0) )
    if self.GetParam(par,is_key=0,is_arg=1) is not None:
      lst.append( self.GetParam(par,is_key=0,is_arg=1) )
    if self.GetParam(par,is_key=1,is_arg=0) is not None:
      lst.append( self.GetParam(par,is_key=1,is_arg=0) )
    return lst

  def CapitalizeArg(self, arg):
    """Capitalize value of the inputted arg inplace (if existing and if possible)"""
    self.CapitalizeParam(arg,value,is_key=0,is_arg=1)

  def CapitalizeKey(self, key):
    """Capitalize value of the inputted key inplace (if existing and if possible)"""
    self.CapitalizeParam(key,value,is_key=1,is_arg=0)

  def CapitalizeVirtPar(self, par):
    """Capitalize value of the inputted virtual parameter inplace (if existing and if possible)"""
    self.CapitalizeParam(par,value,is_key=0,is_arg=0)

  def CapitalizeParam(self, par, is_key=0, is_arg=0):
    # Capitalize value of the inputted parameter inplace (if existing and if possible)
    self.SetParam( par, self.GetParam(par,is_key=is_key,is_arg=is_arg,capitalized=True), 
                   is_key=is_key, is_arg=is_arg )

  def CapitalizeAnyParam(self, par):
    """Searches all program parameters - virtuals, arguments, keywords in this order
       are capitalizes a list of parameters found, empty if no such parameter 'par' found
    """
    if self.GetParam(par,is_key=0,is_arg=0) is not None:
      self.CapitalizeParam(par,is_key=0,is_arg=0)
    if self.GetParam(par,is_key=1,is_arg=0) is not None:
      self.CapitalizeParam(par,is_key=1,is_arg=0)
    if self.GetParam(par,is_key=0,is_arg=1) is not None:
      self.CapitalizeParam(par,is_key=0,is_arg=1)

  def ClearAllKeys(self):
    """Deletes all the currently defined program keywords"""
    self.key={}
    self.key_list=[]

  def ClearAllArgs(self):
    """Deletes all the currently defined program arguments"""
    self.arg={}
    self.arg_list=[]

  def ClearAllVirtPars(self):
    """Deletes all the currently defined program virtual parameters"""
    self.param={}

  def ClearAnyParams(self):
    """Deletes all the currently defined program keywrods, arguments and virtual parameters"""
    self.ClearAllKeys()
    self.ClearAllArgs()
    self.ClearAllVirtPars()

  def InitRunLists(self):
    for arg in self.arg:
      if arg not in self.arg_list:
        self.arg_list.append(arg)
    for key in self.key:
      if key not in self.key_list:
        self.key_list.append(key)

  def IsKey(self,key):
    return self.IsParam(key,is_key=True)

  def IsArg(self,arg):
    return self.IsParam(arg,is_arg=True)

  def IsKeyOrArg(self,arg):
    return self.IsParam(arg,is_arg=True) or self.IsParam(arg,is_key=True)

  def IsVirtPar(self,par):
    return self.IsParam(par)

  def IsParam(self,par,is_key=False,is_arg=False):
    array=self.GetParArray(is_key, is_arg)
    par=self.GetParParsingStyle(par, is_key, is_arg)
    if par in array and array[par] is not [] and array[par] is not None:
      return True
    else:
      return False

  def IsNonFalseKey(self,key):
    return self.IsNonFalseParam(key,is_key=True)

  def IsNonFalseArg(self,arg):
    return self.IsNonFalseParam(arg,is_arg=True)

  def IsNonFalseVirtPar(self,par):
    return self.IsNonFalseParam(par)

  def IsNonFalseParam(self,par,is_key=False,is_arg=False):
    array=self.GetParArray(is_key, is_arg)
    par=self.GetParParsingStyle(par, is_key, is_arg)
    if self.IsParam(par) and array[par] is not False:
      return True
    else:
      return False

  # was param inputted by user?
  def IsInputtedParam(self,par,test_parents=False):
    if self.IsParam(par) and self.param[par].inputted:
      return True
    else:
      if not test_parents or not self.process or not self.process.IsParam(par):
        return False
      else:
        return self.process.IsInputtedParam(par,test_parents=True)

  def Add2Args(self, keyarg, val, args, is_key, add_keyarg=True, multiple=True, ignore_modif=False):
    # adds the parameter 'keyarg' with value 'val' to the arguments array 'args'
    array=self.arg
    if is_key:
      array=self.key
    if multiple:
      array[keyarg].multuse += 1
      if array[keyarg].multuse >= len(array[keyarg].value):
        array[keyarg].multuse = 0
    if val is not False and val is not None:
      try:
        # if multiple values then call recursively
        if common.is_string(val):
          raise TypeError
        for i,v in enumerate(val):
          add_keyarg=True
          if i>0:
            add_keyarg=False
          self.Add2Args(keyarg,v,args,is_key,add_keyarg=add_keyarg,multiple=False)
      except TypeError:
        if add_keyarg:
          modif=''  if is_key or ignore_modif  else self.modif
          value=''  if self.arg_divided or is_key or val is True  else str(val)
          args.append(modif+keyarg+value)
        if val is not True and (self.arg_divided or is_key):
          args.append(str(val))
        if hasattr(sys,'exc_clear'): sys.exc_clear()

  def GetNumParents(self):
    # returns the number of parents of the process (the length of the shortest path to the top of the tree)
    numpar=0
    proc=self.process
    while proc:
      proc=proc.parent_process
      numpar+=1
    return numpar

  def GetCrankParent(self,prep=None):
    # returns the (main) crank process
    proc=self.process
    while proc:
      proc_save=proc
      proc=proc.parent_process
    if proc_save.nick=='crank':
      if prep and hasattr(proc_save,'prep'):
        return proc_save.prep
      else:
        return proc_save
    return None

  def ReplaceInpOut(self, inp1='inp', inp2='inp2', recursive=True):
    # typically used when for hand2 input replacing the current input;  the 'recursive' parameter has no function but keep compatibility with the process analoge
    if hasattr(self,inp1) and hasattr(self,inp2):
      setattr(self,inp1,copy.copy(getattr(self,inp2)))

  def RemoveInpOut(self, inp='inp2', recursive=True):
    # typically used when removing hand2 input
    if hasattr(self,inp):
      delattr(self,inp)

  def AddToRunLists(self, arg_list, key_list, args, scr_lines):
    """Adds specified args/keywords to the args tuple and keyword lines tuple"""
    temp=[]
    for key in key_list:
      if self.IsKey(key):
        self.Add2Args(key, self.GetKey(key,ind=self.key[key].multuse), temp, is_key=1)
        if temp:
          scr_lines.append(' '.join(temp))
          temp=[]
    for arg in arg_list:
      if self.IsArg(arg):
        self.Add2Args(arg, self.GetArg(arg,ind=self.arg[arg].multuse), args, is_key=0, \
          ignore_modif=self.arg[arg].GetModifIgnorance(ind=self.arg[arg].multuse) )

  def GetRunLists(self):
    """Returns a pair consisting of the actual args tuple and the actual keyword lines tuple,
       ready to be inputted to run of the program
    """
    if not self.binary:
      common.Error('No binary defined for program {0}'.format(self.nick))
    args, scr_lines = [self.binary,], []
    self.InitRunLists()
    self.AddToRunLists(self.arg_list, self.key_list, args, scr_lines)
    return (args,scr_lines)

  def GetLogFileName(self):
    return os.path.join(self.rundir, self.runname+self.log_suffix)

  def GetScrFileName(self,relative=None):
    path = os.path.join(self.rundir, self.runname+'_input')
    if relative:
      return os.path.relpath( path, os.path.commonprefix([relative.rundir, path]) )
    else:
      return path

  class ProgramRunError(Exception):
    pass

  def queued_output(self,out,queue):
    for line in iter(out.readline, ''):
      queue.put(line)
    out.close()


  # the lock for interactive processing
  lock_int=threading.Lock()

  # unified routine for external program calling
  def ExternalRun(self, args, inp_scr_lines, clean=0, lock=None):
    #clean = 0 if self.debug>1  else 1
    self.prun=prog_run()
    #clean=0
    if self.rundir is None:
      self.rundir=os.getcwd()
    for e in os.environ:
      if e not in self.env:
        self.env[e]=os.environ[e]
    envs = [ "{0}={1}".format(e,self.env[e])  for e in self.env  if e not in os.environ or os.environ[e]!=self.env[e] ]
    #for i,a in enumerate(args[1:]):
    #  if pwd and len(os.path.commonprefix([pwd,a]))>5 and len(os.path.relpath(a, pwd))<len(a):
    #    print os.path.relpath(a, pwd), a, pwd
    #    args[i+1] = os.path.relpath(a, pwd)
    inp_scr = '\n'.join(inp_scr_lines)
    if inp_scr:  inp_scr+='\n'
    # an obscure fix to an obscure problem of Windows crashes "IOError: [Errno 22] Invalid argument"
    if os.name == 'nt':
      if inp_scr:  inp_scr+=' '
      if len(inp_scr)%2==1:  inp_scr+=' '
    arg_scr = ' '.join(args)
    env_scr = ' '.join(envs)+' '
    if not clean:
      if self.rundir and not os.path.isdir(self.rundir):
        os.mkdir(self.rundir)
      with open(self.GetScrFileName(), 'w') as scrin_f:
        scrin_f.write(env_scr)
        scrin_f.write(arg_scr)
        if inp_scr:
          scrin_f.write(' << END\n')
          scrin_f.write(inp_scr)
          scrin_f.write("END\n")
    if self.log_suffix is not None:
      self.prun.out_f = open(self.GetLogFileName(), 'w', 1)
    # fix for Windows to prevent terminal window being created
    startupinfo=None
    if os.name == 'nt':
      startupinfo = subprocess.STARTUPINFO()
      if sys.version_info == (2, 7):
        startupinfo.dwFlags |= subprocess._subprocess.STARTF_USESHOWWINDOW
      else:
        startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    # interactive output - use with care
    # runs the program in a separate thread, allowing to process its output simultanously
    # the output lines are passed one by one to the Interact_output method of the program
    if self.interact_output:
      self.prun.values=[]
      self.q=Queue()
      self.prun.popen = subprocess.Popen( args, stdin=subprocess.PIPE, stderr=subprocess.PIPE, \
              stdout=subprocess.PIPE, env=self.env, cwd=self.rundir, bufsize=1, startupinfo=startupinfo, universal_newlines=True )
      self.prun.popen.stdin.write(inp_scr)
      self.prun.popen.stdin.close()
      # releasing lock, thus the variables used below should not be local (or another lock would be needed)
      if lock and lock.locked():
        lock.release()
      self.t=threading.Thread( target=self.queued_output, args=(self.prun.popen.stdout,self.q) )
      self.t.daemon=True
      self.t.start()
      while not self.prun.popen.stdout.closed or not self.q.empty():
        if self.q.empty():
          time.sleep(0.1)
        else:
          with program.lock_int:
            self.prun.line = self.q.get_nowait()
            try:
              self.Interact_output(self.prun.line,self.prun.popen,empty=self.q.empty())
              if self.log_suffix is not None:
                self.prun.out_f.write(self.prun.line)
              else:
                self.prun.values.append(self.prun.line)
            # this makes sure that the subprocess dies together with the main program.
            except Exception as e:
              self.prun.popen.terminate()
              raise
      self.prun.values='\n'.join(self.prun.values)
      self.prun.popen.wait()
      self.prun.err=self.prun.popen.stderr.read()
      self.prun.popen.stdout.close()
      self.prun.popen.stderr.close()
    # non-interactive output
    else:
      if self.log_suffix is not None:
        stdo=self.prun.out_f
      else:
        stdo=subprocess.PIPE
      self.prun.popen = subprocess.Popen( args, stdin=subprocess.PIPE, stderr=subprocess.PIPE, \
              stdout=stdo, env=self.env, cwd=self.rundir, startupinfo=startupinfo, universal_newlines=True )
      if lock and lock.locked():
        lock.release()
      self.prun.values,self.prun.err = self.prun.popen.communicate(input=inp_scr)
    self.prun.out_f.close()
    if self.prun.popen.returncode != 0:
      if self.prun.err:
        raise self.ProgramRunError( "{0} failed with error message: {1}".format(self.runname, self.prun.err) )
      else:
        raise self.ProgramRunError( "{0} failed with error code {2} (look at {1} for more information)".format(self.runname, self.GetLogFileName(), self.prun.popen.returncode) )
    if self.log_suffix is None:
      return self.prun.values
    # cleaning prun due to mutiprocessing as popen, out_f are not pickable
    self.prun=prog_run()

class prog_run(object):
  # dummy object collecting the run meta-attributes
  pass
