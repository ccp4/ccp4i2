import os,sys,copy,shutil
import common,data,inout
from xml.etree import ElementTree as ET
from program import program
try:
  import crvapi
except:
  if hasattr(sys,'exc_clear'): sys.exc_clear()
  crvapi = False

class process(object):
  """ Base class for processes such as model building, density modification etc"""
  name="base class for process"
  references=()
  # if no_reporting then no information about the process will be printed in log/loggraph
  no_reporting=False
  # debug=0: no intermediate files are kept
  #       1: intermediate working files are kept
  #       2: all intermediate files are kept, incl. all the temporary scripts
  # THE DEBUG OPTION IS NOT IMPLEMENTED YET!
  debug=0
  # list of supported programs (names)
  supported_progs = []
  # list of supported processes (names)
  supported_procs = []
  # list of supported (virtual) parameters
  supported_params = {}
  supported_params['no_output_to_next_step'] = common.parameter( desc='Disables propagation of output objects to the next step', typ=bool )
  supported_params['no_rewrite'] = common.parameter( desc='Output files names will be adjusted not to rewrite existing files', typ=bool, share=True )
  # log filename extension
  log_suffix='.log'
  # filehandle of the logfile - self.Info() will send output to this filehandle if not None
  logfilehandle=None
  # has *this process* opened logfilehandle or not (even if not opened_log, logfilehandle may be still opened eg by a parent process)
  opened_log=False
  opened_loggraph=False
  # rvapi report (for crank) or section (for other processes)
  rv_report=None

  def __init__(self, xmlelem=None, inpline=None, parent_process=None, rundir=None, 
                     dummy=False, no_support_check=False, extra_attrs=None):
    self.nick = self.__class__.__name__
    self.inp = inout.input_output(is_output=False,parent=self)
    self.out = inout.input_output(is_output=True,parent=self)
    self.rundir=rundir
    self.runname=self.nick
    # actual dictionary of all parameters
    self.param = {}
    # parent process using this process
    self.parent_process = parent_process
    # list of program objects used by the process
    self.programs = []
    # list of process objects used by the process
    self.processes = []
    self.dummy=dummy
    # extra attributes
    if extra_attrs:
      for attr,val in extra_attrs.items():
        setattr(self,attr,val)
    self.init_xmlelem = None
    # initialize by elementree
    if xmlelem is not None:
      self.XMLElemInit(xmlelem,dummy,no_support_check)
    elif inpline is not None:
      self.InputElemInit(inpline,dummy)
    if not dummy:
      self.Init()

  def Init(self):
    """Custom initialization of the process can be placed here"""
    pass

  def PrintParams(self):
    """Prints parameters of the process and their description."""
    if self.supported_params:
      common.Info('The following parameters are defined for process {0}:'.format(self.nick))
      for key,par in self.supported_params.items():
        common.Info('  {0:20} {1}  (type {2})'.format(key, par.description, par.typ[0].__name__))
    else:
      common.Info('No parameters exist for the process.')

  @classmethod
  def from_xml(cls, xml, parent_process=None, rundir=None, dummy=False, no_support_check=False, extra_attrs=None):
    """Create process instance from XML file or XML tree object"""
    if common.is_file(xml) or common.is_string(xml):
      xml = common.ReadXML(xml)
    inst = cls( xml, None, parent_process, rundir, dummy=dummy, no_support_check=no_support_check, extra_attrs=extra_attrs )
    return inst

  @classmethod
  def from_name(cls, procnick, parent, xmlelem=None, inpline=None, dummy=False, no_support_check=False):
    """Create specific (derived) process instance from the (nick)name of the process
       Parameters:
       procnick - the (nickname) process whose instance will be created
       parent - the parent process
       xmlelem - if specified then initialization is performed using this xml element
       inpline - if specified then initialization is performed using this input line (preprocessed list assumed)
    """
    if procnick=='crank':
      from manager import crank
      return crank(xmlelem,inpline,parent,dummy=dummy)
    try:
      # import the specific process
      setattr(sys.modules[__name__], 'processes', __import__('processes.'+procnick))
    except (AttributeError,ImportError):
      common.Error('Process {0} not supported (check the spelling)'.format(procnick))
    try:
      # create an instance of the process
      inst = getattr(sys.modules['processes.'+procnick], procnick)(xmlelem,inpline,parent,dummy=dummy,no_support_check=no_support_check)
    except (AttributeError,KeyError):
      common.Error('Error when creating process {0} instance.'.format(procnick))
    else:
      return inst


  def XMLElemInit(self, xelem, dummy=False, no_support_check=False):
    """Initialize this process from the inputted XML element"""
    param_set_by_crank=[]
    tags=["inp","out","param","program","process","param_set_by_crank"]
    inpnum,outnum,stri,stro=0,0,'',''
    for xchild in xelem:
      if xchild.tag in tags:
        if xchild.tag=="param_set_by_crank":
          param_set_by_crank = common.AutoConvert(xchild.text)
        elif xchild.tag=="inp":
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
        elif xchild.tag=="param":
          for xparam in list(xchild):
            input_by_user = not xparam.tag in param_set_by_crank
            self.SetParam( xparam.tag, common.AutoConvert(xparam.text), inputted=input_by_user )
        elif xchild.tag=="program":
          pr_name=xchild.text.strip()
          if not pr_name in self.supported_progs and not no_support_check:
            common.Error('The program {0} is not supported by {1}.'.format(pr_name,self.name))
          else:
            self.AddProg(pr_name,xchild,dummy=dummy)
        elif xchild.tag=="process":
          pr_name=xchild.text.strip()
          if not pr_name in self.supported_procs and not no_support_check:
            common.Error('The process {0} is not supported by {1}.'.format(pr_name,self.name))
          else:
            self.AddProcess(pr_name,xchild,dummy=dummy,no_support_check=no_support_check)
      else:
        common.Error('Wrong tag when parsing XML for process {0}: {1}'.format(self.name,xchild.tag))
    self.init_xmlelem = xelem

  def InputElemInit(self, input_line, dummy=False):
    """Initialize this process by parsing the input line
       Initializes until a token that cannot be used; consumes the tokens that were used
    """
    first_proc = None
    while input_line:
      token=input_line[0]
      input_line.remove(token)
      if token in self.supported_procs:
        # special treatment needed to disable possibility of multiple crank subprocesses defined in one line by mistake
        if not first_proc:
          first_proc=token
        elif self.nick=='crank':
          common.Warning('Process {0} not supported by {1}'.format(token,first_proc))
          input_line.insert(0,token)
          break
        pr=self.AddProcess(token, inpline=input_line, dummy=dummy)
      elif token in self.supported_progs:
        self.AddProg(token, inpline=input_line, dummy=dummy)
      elif token in [ sc.__name__ for sc in data.data_container.__subclasses__() ]:
        cont=self.inp.AddCopy( data.data_container.from_name(token, inpline=input_line), propagate=not dummy )
        self.CheckInpContFile(cont)
      elif '::' in token:
        key,sep,val = token.partition('::')
        if val:
          self.SetParam( key, common.AutoConvert(val), inputted=True )
      else:
        input_line.insert(0,token)
        break

  def CheckInpContFile(self,cont):
    #global hklin,seqin,xyzin
    # if file was not supplied then the previous mtz or hklin,seqin,xyzin are assumed
    if not cont.GetFile() and cont.GetAllLabels():
      if hasattr(self,'hklin'):
        cont.AddFile(self.hklin, self.hklin.split('.')[-1], adjust_filetype=True)
      else:
        common.Error("No data file associated with object {0} of type {1}".format(cont.nick,cont.typ))
    if not cont.GetFile() and cont.nick=='sequence':
      if hasattr(self,'seqin'):
        cont.AddFile(self.seqin,self.seqin.split('.')[-1])
    if not cont.GetFile() and cont.nick=='model' and cont.typ=='substr' and hasattr(self,'xyzin'):
      cont.AddFile(self.xyzin,self.xyzin.split('.')[-1])
    # save hklin file
    if cont.GetFile() and cont.GetAllLabels():
      self.hklin=cont.GetFileName()
    # copy cell/spgr from plus to minus
    if isinstance(cont,data.fsigf) and cont.GetType()=='minus' and \
       self.inp.Get('fsigf',typ='plus',xname=cont.xname,dname=cont.dname):
      cont.spgr=self.inp.Get(typ='plus',xname=cont.xname,dname=cont.dname).spgr
      cont.cell=self.inp.Get(typ='plus',xname=cont.xname,dname=cont.dname).cell
    # create minus by copy if fsigf "anom" type was used
    if hasattr(cont,'minus_needs_to_be_added'):
      self.inp.AddCopy(cont,allow_duplicates=True).SetType('minus')
      del cont.minus_needs_to_be_added

  def Data2XML(self, ET_element=None):
    """Convert the actual process object into XML"""
    if ET_element is None:
      xroot = ET.Element('process')
    else:
      xroot = ET.SubElement(ET_element,'process')
    xroot.text=self.__class__.__name__
    if self.param:
      param_set_by_crank = [ p for p in self.param if self.IsParam(p) and not self.IsInputtedParam(p) ]
      if param_set_by_crank:
        ET.SubElement(xroot,'param_set_by_crank').text = str(param_set_by_crank)
      xparam = ET.SubElement(xroot,'param')
      for key,par in self.param.items():
        ET.SubElement(xparam, key).text = str(par.value[0])
    hand2io = [getattr(self,io+'2') for io in ['inp','out'] if hasattr(self,io+'2')]
    for obj in [self.inp, self.out] + hand2io + self.processes + self.programs:
      obj.Data2XML(xroot)
    return xroot

  def GetCCP4Header(self):
    # the banner functionality is not (yet?) provided either by a ccp4 python module or
    #  a standalone binary afaik, so this is a reimplementation
    # (calling static C libccp from python would be too much trouble for this single purpose)
    import datetime,getpass
    version=self.GetCrankParent().GetVersion()
    mdate=datetime.datetime.fromtimestamp( os.path.getmtime(sys.modules[self.__class__.__module__].__file__) )
    mdate=str(mdate.day)+'/'+str(mdate.month)+'/'+str(mdate.year)[-2:]
    name='crank2.'+self.nick  if self.nick!='crank'  else 'crank2'
    rundate=datetime.datetime.today()
    runtime=str(rundate.hour)+':'+str(rundate.minute)+':'+str(rundate.second)
    rundate=str(rundate.day)+'/'+str(rundate.month)+'/'+str(rundate.year)
    bars=' ###############################################################\n'
    proc_info=' ### CCP4 {0:3}: {1:21} version {2:7}: {3:8}##\n'.format(
              self.GetCrankParent().GetCCP4Version(), name, version, mdate)
    run_info=' User: {0}  Run date:  {1} Run time: {2}\n'.format(getpass.getuser(),rundate,runtime)
    header=bars+bars+bars+proc_info+bars+run_info+'\n'
    return header

  def Info(self, message, eol=True, stdout=True):
    common.Info(message, self.logfilehandle, eol, also_stdout=stdout)

  def LGInfo(self, message, eol=True, stdout=False):
    if self.GetCrankParent():
      common.Info(message, self.GetLogGraphHandle(), eol, also_stdout=stdout)
    else:
      common.Warning('Could not print loggraph: crank parent does not exist.')

  def GetSupportedParam(self,par,ignore_error=False):
    # checks whether inputted par is amongst supported and completes it if necessary and possible
    # if ignore_error is True then gracefully returns None on issue, otherwise prints the error
    if par in self.supported_params:
      return par
    if par in process.supported_params:
      self.supported_params[par] = process.supported_params[par]
      return par
    # complete or return error message
    sp_lst = [ sp  for sp in self.supported_params  if sp.startswith(par) ]
    if ignore_error and (not sp_lst or len(sp_lst)>1):
      return None
    if not sp_lst:
      common.Error('No such parameter {0} supported by process of {1}'.format(par,self.name))
    if len(sp_lst)>1:
      common.Error('Multiple choices for parameter {0} by process of {1}: {2}'.format(
        par,self.name,', '.join(sp_lst)))
    return sp_lst[0]

  def IsInstance(self,value,typ):
    # somewhat "corrected" isinstance() - bool!=int
    if typ==bool and type(value)==int:
      return False
    if typ==int and type(value)==bool:
      return False
    return isinstance(value,typ)

  def CheckType(self,par,value,accept_none=True):
    # checks type of the parameter to be set - called by SetParam()
    if not [ typ  for typ in self.supported_params[par].typ  if self.IsInstance(value,typ) ]:
      # we do not want to allow all conversions. just those specified here. 
      # for example, I consider str an incorrect input for bool or bool an incorrect input for int.
      if bool in self.supported_params[par].typ and self.IsInstance(value,int):
        value=bool(value)
      elif float in self.supported_params[par].typ and self.IsInstance(value,int):
        value=float(value)
      elif value is not None or not accept_none:
        common.Error( 'Wrong type {0} of value {1} for parameter {2}: must be {3}'.format(
          type(value).__name__, value, par, ' or '.join(t.__name__ for t in self.supported_params[par].typ)) )

  def SetParam(self, par, value=True, no_type_check=False, inputted=False):
    # the same as self.SetVirtPar()
    par = self.GetSupportedParam(par)
    if not no_type_check:
      self.CheckType(par,value)
    if par not in self.param:
      self.param[par] = copy.deepcopy(self.supported_params[par])
    self.param[par].Set(value,is_key=False,is_arg=False,append=False,inputted=inputted)
    # shares the parameter par with children of this process having the same parameter supported
    # if the par has a shared_with_children attribute (recursive)
    if self.param[par].shared_with_children:
      for p in self.GetProcesses() + self.GetProgs():
        if (par in p.supported_params or par in process.supported_params) and not p.IsInputtedParam(par):
          p.SetParam(par, self.GetParam(par), inputted=inputted)

  def SetVirtPar(self, par, value=True):
    """Set program parameter 'par' to inputted value (only 'virtual' parameters allowed for process)"""
    self.SetParam(par,value)

  def GetParam(self, par, capitalized=None, ignore_error=False):
    # the same as self.GetVirtPar()
    par = self.GetSupportedParam(par,ignore_error)
    if par and par in self.param:
      self.param[par].accessed=1
      return self.param[par].Get(capitalized=capitalized)
    else:
      return None

  def GetVirtPar(self, par):
    """Returns the value of parameter 'par', None if not defined"""
    return self.GetParam(par)

  def IsVirtPar(self,par):
    return self.IsParam(par)

  # is param defined?
  def IsParam(self,par):
    if par in self.param and self.GetParam(par) is not [] and self.GetParam(par) is not None:
      return True
    else:
      return False

  # was param inputted by user?
  def IsInputtedParam(self,par,test_parents=False):
    if self.IsParam(par) and self.param[par].inputted:
      return True
    else:
      if not test_parents or not self.parent_process or not self.parent_process.IsParam(par):
        return False
      else:
        return self.parent_process.IsInputtedParam(par,test_parents=True)

  def IsNonFalseVirtPar(self,par):
    return self.IsNonFalseParam(par)

  # is param defined (incl. non-None) and non-False?
  def IsNonFalseParam(self,par):
    if self.IsParam(par) and self.param[par] is not False:
      return True
    else:
      return False

  def IsTrueOrNoneVirtPar(self,par):
    return self.IsTrueOrNoneParam(par)

  # is param not defined or None or True?  (usually used to test whether program defaults should be redefined - if False then orig.prog.def. are forced)
  def IsTrueOrNoneParam(self,par):
    if par not in self.param or self.GetParam(par) is True or self.GetParam(par) is None:
      return True
    else:
      return False

  def CapitalizeVirtPar(self, par):
    """Capitalize value of the inputted virtual parameter (if existing and if possible)"""
    self.CapitalizeParam(par)

  def CapitalizeParam(self, par):
    # Capitalize value of the inputted parameter (if existing and if possible)
    par = self.GetSupportedParam(par)
    self.SetParam( par, self.GetParam(par,capitalized=True) )

  def TreatInOutPar(self,set_all_par=False):
    """processes the actual process input,output and/or parameters"""
    if set_all_par:
      for par in self.supported_params:
        if par not in self.param:
          self.SetParam(par, None)
    # accessing the parameter - the warning is not meant to be printed for this parameter.
    self.GetParam('no_output_to_next_step')


  def CheckParamAccess(self):
    """checks whether all the inputted 'virtual' parameters have been evaluated"""
    for parname,parinst in self.param.items():
      if not parinst.accessed:
        common.Warning('Virtual parameter {0} has not been used by process of {1}.'.format(parname,self.name))

  def GetProgs(self,name=None,supported=False):
    """returns the list of programs (of this process) with specified (nick)name (empty list if there are none)"""
    return [ p  for p in self.programs  if (p.nick==name or name is None) and (not supported or p.nick in self.supported_progs) ]

  def GetProg(self,name=None,supported=False):
    """returns the first program (of this process) with specified name (None if there is none)"""
    try:
      return self.GetProgs(name,supported)[0]
    except IndexError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      return None

  def GetOrAddProg(self,name,*arg,**kwarg):
    """returns the first program object with the specified name or creates a new one if none found"""
    prog = self.GetProg(name)
    if prog is None:
      prog = self.AddProg(name,*arg,**kwarg)
    return prog

  def GetProcesses(self,name=None,supported=False):
    """returns the list of processes (of this process) with specified (nick)name(s) (empty list if there are none)"""
    if common.is_string(name) or name is None:
      name=[name,]
    return [ p  for p in self.processes  for n in name if (p.nick==n or n is None) and (not supported or p.nick in self.supported_procs) ]

  def GetProcess(self,name=None,supported=False):
    """returns the first process (of this process) with specified name (None if there is none)"""
    try:
      return self.GetProcesses(name,supported)[0]
    except IndexError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      return None

  def GetOrAddProcess(self,name,*arg,**kwarg):
    """returns the first subprocess with the specified name or creates a new one if none found"""
    proc = self.GetProcess(name)
    if proc is None:
      proc = self.AddProcess(name,*arg,**kwarg)
    return proc

  def GetFlatSubTree(self,lst=None):
    if lst is None:
      lst=[]
    lst.append(self)
    for prog in self.programs:
      prog.GetFlatSubTree(lst)
    for proc in self.processes:
      proc.GetFlatSubTree(lst)
    return lst

  def AddProg(self, name, xmlelem=None, inpline=None, ind=None, propagate_inp=True, propagate_out=True, dummy=False):
    """Adds a new program to this process.

       Parameters:
       name - (nick)name of the program to be added
       xmlelem - xml element used for initialization of the added program object (None by default)
       inpline - input line (preprocessed list) used for initialization of the program (None by default)
       ind - add to this index in the programs list of this process (after the last by default)
       propagate_inp - if True then the added program's input is propagated from the self process (default)
                       if False then nothing is propagated to the input
       propagate_out - if True then the output will be propagated (default)
                       if False then the output is flagged as never_propagate 
       dummy - create a 'dummy' object - without its specific initialization and turns propagate to False
    """
    if dummy:
      propagate_inp=False
      propagate_out=False
    if ind is None:
      ind=len(self.programs)
    self.programs.insert(ind, program.from_name(name,self,xmlelem,inpline,dummy))
    if propagate_inp:
      for container in self.inp.GetAll():
        # propagating to the beginning of the list so that the program's own objects have priority
        self.programs[ind].inp.Add(container,ind=0)
      # propagate inp2 etc
      inpind=2
      while hasattr(self,'inp'+str(inpind)):
        for container in getattr(self,'inp'+str(inpind)).GetAll():
          if not hasattr(self.programs[ind],'inp'+str(inpind)):
            setattr(self.programs[ind],'inp'+str(inpind), inout.input_output(is_output=False,parent=self.programs[ind]))
          getattr(self.programs[ind],'inp'+str(inpind)).Add(container,ind=0)
        inpind+=1
      inpind=2
      # propagate inp also to inp2 etc if self has a single inp, but the child has multiple inps
      while not hasattr(self,'inp'+str(inpind)) and hasattr(self.programs[ind],'inp'+str(inpind)):
        for container in getattr(self,'inp').GetAll():
          getattr(self.programs[ind],'inp'+str(inpind)).Add(container,ind=0)
        inpind+=1
    if not propagate_out:
      self.programs[ind].out.never_propagate=True
    if self.programs[ind].always_rundir:
      self.programs[ind].SetRunDir(error_on_failure=False)
    # share the value of the parameters with .shared_with_children attribute with the added program child
    for par in self.programs[ind].supported_params:
      if self.IsParam(par) and self.param[par].shared_with_children and not self.programs[ind].IsInputtedParam(par):
        self.programs[ind].SetParam(par, self.GetParam(par))
    # disable output mtz label prefixes if requested
    if (self.GetCrankParent() and hasattr(self.GetCrankParent(),'disable_mtz_label_prefix')) or \
        self.programs[ind].GetParam('ignore_labelout_prefix'):
      self.programs[ind].labelout_prefix=''
    return self.programs[ind]

  def AddProcess(self, name, xmlelem=None, inpline=None, ind=None, propagate_inp=True, propagate_out=True, 
                       dummy=False, no_support_check=False):
    """Adds a new subprocess to this process. The same parameters as AddProg()."""
    if dummy:
      propagate_inp=False
      propagate_out=False
    if ind is None:
      ind=len(self.processes)
    self.processes.insert(ind, process.from_name(name,self,xmlelem,inpline,dummy,no_support_check))
    if propagate_inp:
      for container in self.inp.GetAll():
        # propagating to the beginning of the list so that the process's own objects have priority
        self.processes[ind].inp.Add(container,ind=0)
      # propagate inp2 etc
      inpind=2
      while hasattr(self,'inp'+str(inpind)):
        for container in getattr(self,'inp'+str(inpind)).GetAll():
          if not hasattr(self.processes[ind],'inp'+str(inpind)):
            setattr(self.processes[ind],'inp'+str(inpind), inout.input_output(is_output=False,parent=self.processes[ind]))
          getattr(self.processes[ind],'inp'+str(inpind)).Add(container,ind=0)
        inpind+=1
      inpind=2
      # propagate inp also to inp2 etc if self has a single inp, but the child has multiple inps
      while not hasattr(self,'inp'+str(inpind)) and hasattr(self.processes[ind],'inp'+str(inpind)):
        for container in getattr(self,'inp').GetAll():
          getattr(self.processes[ind],'inp'+str(inpind)).Add(container,ind=0)
        inpind+=1
    if not propagate_out:
      self.processes[ind].out.never_propagate=True
    # share the value of the parameters with .shared_with_children attribute with the added child
    for par in list(self.processes[ind].supported_params.keys())+list(process.supported_params.keys()):
      if self.IsParam(par) and self.param[par].shared_with_children and not self.processes[ind].IsInputtedParam(par):
        self.processes[ind].SetParam(par, self.GetParam(par))
    return self.processes[ind]

  def AddProgCopy(self, prog, ind=None, propagate_inp=False, deeper_copy=False):
    """Adds a copy of the inputted prog program object to this process. 
       No input objects are propagated to the added program from the parent process by default.
       The output propagation is the same as that of the program that was copied.

       Parameters:
       prog - the program object a copy of which will be added to this process
       ind - add to this index in the programs list of this process (after the last by default)
       propagate_inp - if True then the added program's input is propagated from the self process
                       if False then nothing is propagated to the input (default)
       deeper_copy: if True then the input/output and args/keys are copied 
                    if False then the original proc's i/o,args/keys objects are kept (ie shared with prog)
                    (False by default)
    """
    assert isinstance(prog, program)
    if ind is None:
      ind=len(self.programs)
    self.programs.insert( ind, copy.copy(prog) )
    self.programs[ind].process=self
    if deeper_copy:
      self.programs[ind].inp, self.programs[ind].out = copy.copy(prog.inp), copy.copy(prog.out)
      self.programs[ind].inp.ClearAll(propagate=False),  self.programs[ind].out.ClearAll(propagate=False)
      for container in prog.inp.GetAll(stored_order=True):
        self.programs[ind].inp.Add(container,propagate=False)
      for container in prog.out.GetAll(stored_order=True):
        self.programs[ind].out.Add(container,propagate=False)
      self.programs[ind].inp.parent = self.programs[ind].out.parent = self.programs[ind]
      self.programs[ind].param = copy.copy(prog.param)
      self.programs[ind].key, self.programs[ind].arg = copy.deepcopy(prog.key), copy.deepcopy(prog.arg)
      self.programs[ind].key_list = copy.deepcopy(prog.key_list)
      self.programs[ind].arg_list = copy.deepcopy(prog.arg_list)
    if propagate_inp:
      for container in self.inp.GetAll():
        # propagating to the beginning of the list so that the process's own objects have priority
        self.programs[ind].inp.Add(container,ind=0)
    return self.programs[ind]

  def AddProcessCopy(self, proc, ind=None, propagate_inp=False, deeper_copy=False):
    """Adds a copy of the proc process object to this process. The same parameters as AddProgCopy() +
       deeper_copy: if True then all the subprocesses/subprograms,input/output,param are copied 
                    if False then the original proc's subprocesses/subprograms,i/o,par are kept (ie shared with proc)
                    (False by default)
    """
    assert isinstance(proc, process)
    if ind is None:
      ind=len(self.processes)
    self.processes.insert( ind, copy.copy(proc) )
    self.processes[ind].parent_process=self
    if deeper_copy:
      self.processes[ind].inp = copy.copy(proc.inp)
      self.processes[ind].out = copy.copy(proc.out)
      self.processes[ind].inp.ClearAll(propagate=False),  self.processes[ind].out.ClearAll(propagate=False)
      for container in proc.inp.GetAll(stored_order=True):
        self.processes[ind].inp.Add(container,propagate=False)
      for container in proc.out.GetAll(stored_order=True):
        self.processes[ind].out.Add(container,propagate=False)
      self.processes[ind].inp.parent = self.processes[ind].out.parent = self.processes[ind]
      self.processes[ind].param = copy.deepcopy(proc.param)
      self.processes[ind].processes, self.processes[ind].programs = [], []
      for pr in proc.processes:
        prcopy = self.processes[ind].AddProcessCopy(pr, deeper_copy=True)
        # fixing attribute assignements to the nested processes:  however, only works for child processes assigned at the parent level
        [setattr(self.processes[ind],atr,prcopy)  for atr in proc.__dict__  if getattr(proc,atr) is pr]
      for pr in proc.programs:
        prcopy = self.processes[ind].AddProgCopy(pr, deeper_copy=True)
        [setattr(self.processes[ind],atr,prcopy)  for atr in proc.__dict__  if getattr(proc,atr) is pr]
    if propagate_inp:
      for container in self.inp.GetAll():
        # propagating to the beginning of the list so that the process's own objects have priority
        self.processes[ind].inp.Add(container,ind=0)
    return self.processes[ind]

  def GetProcessCopy(self):
    # returns a copy of this process by conversion to and from xml.
    # does not copy all the meta-attributes: use with care!
    proc_copy = self.from_xml(self.Data2XML(), no_support_check=True)
    proc_copy.rundir = self.rundir
    proc_copy.logfilehandle = self.logfilehandle
    if hasattr(self,'run'):
      proc_copy.run = self.run.GetProcessCopy()
    #proc_copy = copy.copy(self)
    #proc_copy.param = {}
    #for key,par in self.param.items():
    #  proc_copy.param[key] = par
    #proc_copy.processes=[]
    #for proc in self.processes:
    #  proc_copy.processes.append( proc.GetProcessCopy() )
    #proc_copy.programs=[]
    #for prog in self.programs:
    #  proc_copy.programs.append( prog.GetProgramCopy() )
    return proc_copy

  def SetRunDir(self, rundir=None, reset_subtree=False, change_cwd=False):
    if rundir is not None:
      self.rundir=rundir
    else:
      if self.rundir is None:
        if self.parent_process and self.parent_process.rundir:
          self.rundir=self.parent_process.rundir
          if len(self.parent_process.processes+self.parent_process.programs)>1:
            self.rundir=os.path.join(self.rundir,self.nick)
        else:
          self.rundir=os.getcwd()
    if not os.path.isdir(self.rundir):
      os.mkdir(self.rundir)
    if not os.path.isabs(self.rundir):
      self.rundir=os.path.abspath(self.rundir)
    if reset_subtree:
      for pr in self.processes:
        pr.rundir=None
        pr.SetRunDir(reset_subtree=True)
      for pr in self.programs:
        pr.rundir=None
        pr.SetRunDir()
    if change_cwd:
      self.previous_cwd=os.getcwd()
      os.chdir(self.rundir)

  def GetLogFileName(self):
    verb = ''  if self.nick!='crank'  else '_verb'
    return os.path.join(self.rundir, self.runname + verb + self.log_suffix)

  def GetLogGraphFileName(self, actual=False):
    """Returns the loggraph filename for this process.
       For crank process, this is the 'merge' loggraph filename; otherwise the 'actual' loggraph filename.
       Returns None if crank parent does not exist or loggraph disabled by user using --lgout None
       If the parameter 'actual' is True then the actual filename in crank directory is retuned (used to 
         be default with crank2<2.0.10 and still used for gui2 as symlink)
    """
    crank = self.GetCrankParent()
    if crank:
      if actual:
        lgdir, lgname, lgsuffix = crank.rundir, crank.runname, '.loggraph'
      else:
        lgdir, lgname, lgsuffix = self.rundir, self.runname, '.loggraph'
      if hasattr(self,'lgout') or (hasattr(crank,'lgout') and actual):
        if crank.lgout == 'None':
          if self.nick=='crank':
            common.Info('Loggraph output disabled.'.format(lgdir))
          return None
        lgdir=os.path.dirname(crank.lgout)
        lgname,lgsuffix=os.path.splitext(os.path.basename(crank.lgout))
      actual_str='_act'
      if self.nick=='crank' or not actual:
        actual_str=''
      return os.path.join(lgdir, lgname + actual_str + lgsuffix)
    else:
      return None

  def GetLogGraphHandle(self):
    """Gets the loggraph handle of this process.
       For crank process, this is the 'merge' loggraphfile handle; otherwise the 'actual' handle.
    """
    if self.nick=='crank':
      return self.loggraphfilehandle
    elif self.GetCrankParent():
      return self.GetCrankParent().loggraph_actfilehandle
    else:
      return None

  def SetLogGraphHandle(self, handle):
    """Sets the loggraph handle of this process to the inputted handle.
       For crank process, sets the 'merge' loggraphfile handle; otherwise the 'actual' handle
    """
    if self.nick=='crank':
      self.loggraphfilehandle = handle
    else:
      self.GetCrankParent().loggraph_actfilehandle = handle

  def GetLogGraphReferences(self, crank_only=True, recursive=True, rvapi=False, bib_to_file=False):
    """Returns list of references strings for the current process if rvapi
       Returns a loggraph formatted string of references if not rvapi
       Returns empty list/string if no references are defined.
       Arguments:
       crank_only - if True then returns empty list/string for any other processes than crank (default True)
       recursive - if true then searches all the subprocesses recursively for references as well
    """
    if crank_only and self.nick!='crank':
      return '' if not rvapi  else []
    try:
      instance=self.emulate.run
    except:
      instance=self
    pros = [instance,] + instance.programs
    if recursive:
      pros = instance.GetAllSubPrograms() + instance.GetAllSubProcesses()
    # the SHELX reference is currently the one of SHELXE
    if self.GetProcess('phdmmb') and not recursive:
      pros.insert(0,program.from_name('shelxe',None))
      #section.Text('&emsp;'+program.from_name('shelxe',None).references[0])
    # the Crank2 reference for MR-SAD
    if self.GetProcess('refatompick') and not recursive:
      pros.insert(0,program.from_name('crank2_mrsad',None))
    refs = [ p.references[0]  for p in pros  if p.references ]
    # make unique
    refs = [ r  for i,r in enumerate(refs)  if refs.index(r)==i ]
    if refs and bib_to_file and self.GetCrankParent():
      ref_dir = os.path.join(os.path.dirname(__file__),'programs','references')
      done=[]
      with open(os.path.join(self.GetCrankParent().rundir,'crank2.bibtex.txt'),'w' if pros[0].nick=='crank' else 'a') as g1:
       with open(os.path.join(self.GetCrankParent().rundir,'crank2.medline.txt'),'w' if pros[0].nick=='crank' else 'a') as g2:
        for p in pros:
          if os.path.isfile(os.path.join(ref_dir,p.nick+'.bib')):
            with open(os.path.join(ref_dir,p.nick+'.bib')) as f:
              ref_bib = f.read()
              if not ref_bib or ref_bib in done:
                continue
              done.append( ref_bib )
              g1.write( ref_bib+'\n' )
          if os.path.isfile(os.path.join(ref_dir,p.nick+'.nbib')):
            with open(os.path.join(ref_dir,p.nick+'.nbib')) as f:
              g2.write( f.read()+'\n' )
    if rvapi:
      return refs
    else:
      return '$TEXT:Reference: $$ Please reference $$\n{0}\n$$\n'.format('\n'.join(refs))

  def GetLogGraphPrograms(self,rvapi=False,binaries=False):
    """Returns list of program (name) strings used by the process if rvapi
       Returns loggraph formatted string containing information about programs used by the process.
       Returns list of program binaries rather than names if the binaries parameter is True"""
    try:
      programs=self.emulate.run.GetAllSubPrograms()
    except:
      programs=self.GetAllSubPrograms()
    # make list of unique program names
    programs = [p.name if not binaries else p.binary for p in programs]
    programs = [p for i,p in enumerate(programs) if programs.index(p)==i]
    if 'sftools' in programs:
      programs.remove('sftools')
    if rvapi:
      return programs
    else:
      if not programs:
        return ''
      return '$TEXT:Programs used: $$ Programs used in the run $$\n{0}\n$$\n'.format(', '.join(programs))

  def GetAllSubPrograms(self, programs=None):
    """Returns list of all (currently defined) programs of this process or its subprocesses (recursively)"""
    if programs is None:
      programs=[]
    programs.extend(self.programs)
    for proc in self.processes:
      proc.GetAllSubPrograms(programs)
    return programs

  def GetAllSubProcesses(self, processes=None):
    """Returns list of all (currently defined) subprocesses of this process or its subprocesses (recursively)"""
    if processes is None:
      processes=[]
    processes.extend(self.processes)
    for proc in self.processes:
      proc.GetAllSubProcesses(processes)
    return processes

  def ReplaceInpOut(self, inp1='inp', inp2='inp2', recursive=True):
    # typically used when for hand2 input replacing the current input
    if hasattr(self,inp1) and hasattr(self,inp2):
      keep_order_num = getattr(self,inp1).order_num
      setattr(self,inp1, copy.copy(getattr(self,inp2)))
      getattr(self,inp1).order_num=keep_order_num
    if recursive and hasattr(self,inp1):
      for io in getattr(self,inp1).FindPropagateObjects():
        io.ReplaceInpOut(inp1,inp2,recursive=True)

  def RemoveInpOut(self, inp='inp2', recursive=True):
    # typically used when removing hand2 input
    if hasattr(self,inp):
      if recursive:
        for io in getattr(self,inp).FindPropagateObjects():
          io.RemoveInpOut(inp,recursive=True)
      delattr(self,inp)

  def BackupAnyPars(self):
    # backups all actual parameters into private tmp arrays
    self.partmp = copy.deepcopy(self.param)
    return self.partmp

  class RestoreParamsError(Exception):
    pass

  def RestoreAnyPars(self,par=None):
    # restores the previously backuped parameters from private tmp arrays
    if par:
      self.param = copy.deepcopy(par)
    else:
      try:
        self.param = copy.deepcopy(self.partmp)
      except AttributeError:
        raise self.RestoreParamsError('Restoring backup parameters failed, backup does not exist.')

  def GetNumParents(self):
    # returns the number of parents of the process (the length of the shortest path to the top of the tree)
    numpar=0
    proc=self.parent_process
    while proc:
      proc=proc.parent_process
      numpar+=1
    return numpar

  def GetCrankParent(self,prep=None):
    # returns the (main) crank process
    proc=self
    while proc:
      proc_save=proc
      proc=proc.parent_process
    if proc_save.nick=='crank':
      if prep and hasattr(proc_save,'prep'):
        return proc_save.prep
      else:
        return proc_save
    return None

  def RunPreprocess(self, rundir=None, clear_out=True, save=True, *args, **kwargs):
    """Any preprocessing before the actual run of the process (including TreatInOutPar)"""
    setallpar = kwargs.pop('set_all_par', False)
    from_ccp4i2 = kwargs.pop('from_ccp4i2', False)
    emulate = kwargs.pop('emulate', '')
    info = kwargs.pop('info', '')
    global crvapi
    crvapi = kwargs.pop('crvapi', crvapi)
    self.SetRunDir(rundir, change_cwd=True)
    if self.GetCrankParent() and self.GetCrankParent().logfilehandle:
      self.GetCrankParent(prep=True).Info('{0}Starting process of {1} {2}'.format('  '*self.GetNumParents(), self.name, info), stdout=False)
    if not emulate and self.nick=='crank':
      self.CheckBinaries()
    if clear_out:
      self.out.ClearAll(propagate=False)
    if save:
      self.partmp=copy.copy(self.param)
    if self.nick!='crank':
      self.TreatInOutPar(set_all_par=setallpar)
    if not self.no_reporting:
      if self.logfilehandle is None or self.logfilehandle.closed:
        self.logfilehandle=open(self.GetLogFileName(),'ab',0)
        self.opened_log=True
        if from_ccp4i2 and os.path.isfile(self.GetLogFileName()):
          # create symlinks to the logfile - as of now required by ccp4i2 to register the logfile...
          common.SymLink(self.GetLogFileName(), os.path.join(self.rundir, 'log.txt'))
      if crvapi and self.GetCrankParent():
        section=None
        if self.nick=='crank':
          #self.rv_is_tree = self.run.rv_is_tree = crvapi.tree  if hasattr(crvapi,'tree')  else True
          #self.rv_separate_steps = self.run.rv_separate_steps = False
          viewer = self.rvapi_viewer  if hasattr(self,'rvapi_viewer')  else None
          rvdoc = self.rvapi_document  if hasattr(self,'rvapi_document')  else None
          # full layout for non i2;  nothing for i2 tree;  otherwise tabs in i2 (sections in grid do not work in rvapi)
          layout = 7 if not self.run.ccp4i2 else 0 if crvapi.tree else 4
          crank2doc = crvapi.Document( ID='crank2doc', outdir=self.rundir, header='CRANK2 job viewer', \
                                       viewer=viewer, rvapidoc=rvdoc, i2=self.run.ccp4i2, layout=layout)
          if self.run.ccp4i2 and crvapi.tree:
            self.rv_body = self.rv_report = crank2doc.Grid("body")
          else:
            self.rv_body = self.rv_report = crank2doc.Tab("results_tab", "Results")
            if hasattr(self,'logout'):
              self.rv_log = crank2doc.Tab("log_tab", "Logfile")
              crvapi.AppendContent(self.logout, "log_tab")
          title = "SHELX (run via CRANK2)"  if self.GetProcess('phdmmb')  else "CRANK2 {0}".format(self.GetVersion())
          section = self.rv_report.Section('crank2', title=title, opened=False)
          section.Text('<i><small>Please cite:</small></i>')
          for r in self.GetLogGraphReferences(recursive=False,rvapi=True,bib_to_file=self.run.ccp4i2):
            section.Text('&emsp;'+r)
          section.Text('<BR>')
          binaries={'programs_used':self.GetLogGraphPrograms(rvapi=True,binaries=True)}
          import json
          crvapi.PutMetaData( json.dumps(binaries), rvdoc )
          if crvapi.tree:
            self.rv_report = self.rv_report.Tree("Processes")
        elif self.parent_process and self.parent_process.nick=='crank' and self.GetCrankParent().rv_report \
             and not self.no_reporting:
          sec_name = self.short_name if crvapi.tree else self.name
          if crvapi.separate_steps:
            rv_proc = crvapi.Document( ID=self.nick+'doc', outdir=self.rundir, i2=True, layout=4 ).Tab(self.nick+"_tab", "Results")
            row = -1
          else:
            rv_proc = self.GetCrankParent().rv_report
            row = self.GetCrankParent().processes.index(self)+1
          #sec_name = self.short_name if crvapi.tree else self.name
          #section = self.rv_report = self.GetCrankParent().rv_report.Section(self.nick,sec_name[0].upper()+sec_name[1:])
          section = self.rv_report = rv_proc.Section(self.nick,sec_name[0].upper()+sec_name[1:])#,row=row)
        if section is not None:
          section.Text("<i><small>Programs used:</small></i>")
          section.Text('&emsp;'+', '.join(self.GetLogGraphPrograms(rvapi=True)))
          if self.nick=='crank':
            section.Text("<BR><i><small>Program references:</small></i>")
            for r in self.GetLogGraphReferences(rvapi=True,bib_to_file=self.run.ccp4i2):
              section.Text('&emsp;'+r)
            if not self.run.ccp4i2:
              section.Text("<BR><i><small>Input data:</small></i>")
              for fname in set([o.GetFileName() for o in self.inp.GetAll('fsigf')]):
                section.DataFile(fname, "hkl:hkl", "<small>Input data</small>")
          else:
            crvapi.Flush(ignore_timing_restr=True)
      # at the moment, there is one "previous steps' merge" loggraphfile and the actual step loggraphfile
      # both of them are opened from here; the actual step one is only opened if not opened yet by another 
      # process and if the main merge one is opened
      # currently, only the step processes (direct crank's children processes) have a loggraphfile
      if self.GetLogGraphFileName() and self.GetLogGraphHandle() is None and \
         (self.nick=='crank' or self.GetCrankParent().opened_loggraph):
        self.SetLogGraphHandle( open(self.GetLogGraphFileName(),'wb',0) )
        self.opened_loggraph=True
        # create symlink as _act to the main crank directory (used by both ccp4i1 and ccp4i2)
        if self.nick!='crank':
          common.SymLink(self.GetLogGraphFileName(), self.GetLogGraphFileName(actual=True))
        self.LGInfo(self.GetCCP4Header())
        self.LGInfo(self.GetLogGraphReferences(recursive=False))
        self.LGInfo(self.GetLogGraphPrograms())
        self.LGInfo(self.GetLogGraphReferences())
      if self.parent_process is self.GetCrankParent():
        common.Info('\n*** Starting process of {0} ***\n'.format(self.name))

  def RunBody(self,*args,**kwargs):
    """Actual run of the process, by default just calling the associated program(s, in sequence)"""
    prev_prog=None
    for prog in self.programs:
      if prev_prog:
        for container in prev_prog.out.GetAll():
          prog.inp.AddCopy(container)
      prog.Run()
      prev_prog=prog

  def RunPostprocess(self,restore=True,*args,**kwargs):
    """Any postprocessing after the actual run of the process"""
    from_ccp4i2 = kwargs.pop('from_ccp4i2', False)
    info = kwargs.pop('info', '')
    if restore:
      self.param=copy.copy(self.partmp)
    self.CheckParamAccess()
    os.chdir(self.previous_cwd)
    if self.GetCrankParent():
      self.GetCrankParent(prep=True).Info('{0}Process of {1} {2} has finished'.format('  '*self.GetNumParents(),self.name,info), stdout=False)
    if self.opened_log:
      if self.parent_process is self.GetCrankParent():
        common.Info('\n*** Process of {0} has finished. ***\n\n'.format(self.name))
      self.logfilehandle.close()
      self.logfilehandle=None
      self.opened_log=False
      if os.path.getsize(self.GetLogFileName())==0:
        os.remove(self.GetLogFileName())
    if self.opened_loggraph:
      self.GetLogGraphHandle().close()
      self.SetLogGraphHandle(None)
      self.opened_loggraph=False
      # the content of the "actual" loggraphfile is copied to the main crank loggraphfile
      if self.nick!='crank':
        if os.path.getsize(self.GetLogGraphFileName())>0:
          with open(self.GetLogGraphFileName(),'rb') as act_loggraph:
            shutil.copyfileobj(act_loggraph, self.GetCrankParent().loggraphfilehandle)
        # fix for Windows: otherewise file not found and crank crashes
        fpath = self.GetLogGraphFileName(actual=True)
        if os.path.isfile(fpath):
          os.remove(fpath)
        if not self.no_reporting and crvapi and not crvapi.tree and not crvapi.separate_steps and self.rv_report:
          # backwards compatibility - as Eugene explained, the fallback should not be used, thus should be removed once distributed jsrview rev>=104 can be assumed
          if not self.rv_report.SetState(close=True):
            self.GetCrankParent().rv_report.Section(self.nick,"",row=20,opened=False)
      elif not self.no_reporting and crvapi:
        build = [p for p in self.run.processes if p.nick in ('comb_phdmmb','mbref')]
        handchanger = [p for p in self.run.processes if p.nick in ('phdmmb','handdet','comb_phdmmb','mbref') and hasattr(p,'spacegroup_change') and p.spacegroup_change]
        if build:
          summary = self.rv_body.Section("summary","Job Summary",row=0)
          summary.Text('<i>'+build[-1].result_str+'</i>')
          summary.Text('Number of residues in the built model: {0}'.format(build[-1].mb_res_all[-1][1]))
          if hasattr(self.run.processes[-1],'report_R'):
            summary.Text('Final R factor: {0}'.format(self.run.processes[-1].report_R))
          if hasattr(self.run.processes[-1],'report_Rfree') and self.run.processes[-1].report_Rfree:
            summary.Text('Final Rfree factor: {0}'.format(self.run.processes[-1].report_Rfree))
          if handchanger:
            summary.Text('<p style="color:DarkOrange">Spacegroup set to: {0}</p>'.format(handchanger[-1].spacegroup_change))
          if hasattr(self,'xyzout') and hasattr(self,'hklout') and not self.ccp4i2:
            # mainly for jsCoFe.  the requested files were prepared in manager PostProcess.
            summary.Text("<BR><i><small>Output:</small></i>")
            summary.rv_files = summary.DataFile(self.xyzout, "xyz", "<small>Built model and map</small>")
            summary.rv_files.DataFile(self.hklout, "hkl:map")
            #summary.rv_files.DataFile(self.hklout+'.map', "hkl:ccp4_map")
            #summary.rv_files.DataFile(self.hklout+'_diff.map', "hkl:ccp4_dmap")
          crvapi.Flush(ignore_timing_restr=True)

  def Run(self,rundir=None,*args,**kwargs):
    """Run the process!"""
    clearout=kwargs.pop('clear_out',True)
    restore=kwargs.pop('restore',False)
    self.RunPreprocess(rundir=rundir,clear_out=clearout,save=restore,*args,**kwargs)
    self.RunBody(*args,**kwargs)
    self.RunPostprocess(restore,*args,**kwargs)
