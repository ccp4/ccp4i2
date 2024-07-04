import os,sys,copy
from xml.etree import ElementTree as ET
import common



class files(object):
  # mtz labels handling should be moved here (?)

  def __init__(self,fname,ftyp):
    self.name=fname
    self.typ=ftyp

  def Data2XML(self, ET_element):
    """Convert the actual file object into XML"""
    xroot = ET.SubElement(ET_element, 'file')
    xroot.text = self.name
    xroot.attrib['typ']=self.typ

  def Convert(self,newtyp,par_proc,inp_cont,conv_opts=[]):
    if all(c.GetFile(newtyp) for c in inp_cont):
      return inp_cont[0].GetFile(newtyp)
    from_to=(self.typ,newtyp)
    # in case we are called from program, find its parent process
    try:
      par_proc=par_proc.process
    except AttributeError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
    # remember objects with mtzfile_asked initially, to be restored after conversion
    mtzfile_asked_objs = [o for o in par_proc.inp.GetAll(try_convert=False) if o.mtzfile_asked]
    # create the convert process and run it
    convert_proc=par_proc.AddProcess('convert',propagate_inp=False,propagate_out=False)
    convert_proc.no_conversion = par_proc.no_conversion if hasattr(par_proc,'no_conversion') else False
    if from_to not in convert_proc.setup:
      par_proc.processes.remove(convert_proc)
      return None
    # we have to addcopy rather than add as the conversion can fail
    cont_dict={}
    for cont in inp_cont:
      cont_dict[cont] = convert_proc.inp.AddCopy(cont)
    convert_proc.opts = conv_opts[:]
    if 'no_rewrite' in conv_opts:
      convert_proc.SetParam('no_rewrite')
    convert_proc.SetRunDir()
    converted_file=None
    if not convert_proc.no_conversion:
      #suppressing any warnings (this could be made optional)
      with open(os.devnull, 'w') as fnull:
        sys_stdout, sys_stderr = sys.stdout, sys.stderr
        sys.stdout = sys.stderr = fnull
        try:
          convert_proc.setup[from_to](convert_proc,*from_to)
          if not convert_proc.no_conversion:
            convert_proc.Run(clear_out=False,info='from '+from_to[0]+' to '+from_to[1])
        except Exception as e:
          sys.stdout, sys.stderr = sys_stdout, sys_stderr
          #common.Warning('Conversion from {0} to {1} of process {2} failed with message: {3}'.format(self.typ,newtyp,par_proc.nick,e))
          if hasattr(sys,'exc_clear'): sys.exc_clear()
        else:
          sys.stdout, sys.stderr = sys_stdout, sys_stderr
          # success - copy the objects (inplace conversion assumed, ie the converted format must be returned in the input container(s))
          for cont in cont_dict:
            for k,v in cont_dict[cont].__dict__.items():
              setattr(cont,k,copy.deepcopy(v))
    # clean up and return
    par_proc.processes.remove(convert_proc)
    for o in par_proc.inp.GetAll(try_convert=False):
      if o.mtzfile_asked and o not in mtzfile_asked_objs:
        o.mtzfile_asked=False
    return converted_file



class data_container(object):
  """base container for data (observed and calculated) storage classes"""
  # the internal booking changes - cell, spacegroup etc might be treated as property of crystal etc
  # mtz column labels should be hidden on the files level etc
  default_unknown='Undefined'
  deffiletype=None
  supported_filetypes=[]
  xname=default_unknown
  dname=default_unknown
  cell=None
  spgr=None
  spgr_num=None
  resol=None
  resol_low=None
  wavel=None
  col_list=[]
  label_default=None
  labelout_prefix=''
  custom=[]
  typ=''
  solvent_content=None
  # number of (equal) monomers in the asym unit.  This gets tricky if sequence is inputted with
  # (equal) monomers in it:  monomers_asym actually considers the entire inputted sequence as one monomer.
  # therefore, another parameter seq_monomers is needed, defining the number of monomers within the sequence
  monomers_asym=None
  seq_monomers=None
  residues_mon=None
  # list of supported arguments set by user (as opposed to set by crank)
  # this might get tricky when new objects are created from objects with user inputted args...
  set_by_user = []

  def __init__(self, filename=None, filetype=None, typ=None, lab_dict={}, custom=None, \
                     xname=None, dname=None, cell=None, spgr=None, resol=None, xmlelem=None, inpline=None):
    self.nick=self.__class__.__name__
    self.mtzfile_asked=False
    self.mtzfile_asked_before={}
    if typ is not None:
      self.typ=typ
    self.files=[]
    if filename is not None:
      self.AddFile(filename, filetype)
    if self.col_list:
      self.InitLabels(lab_dict)
    self.custom=[]
    if custom:
      self.custom.append(custom)
    if xname:
      self.xname=str(xname)
    if dname:
      self.dname=str(dname)
    if cell:
      self.cell=cell
    if spgr:
      self.spgr=spgr.replace(' ','')
    if resol:
      self.resol=resol
    self.delayed_init={}
    if xmlelem is not None:
      self.XMLElemInit(xmlelem)
    elif inpline is not None:
      self.InputElemInit(inpline)

  @classmethod
  def from_name(cls, container_name, xmlelem=None, inpline=None):
    """Create specific (derived) container instance from the (nick)name of the container"""
    try:
      inst = getattr(sys.modules[__name__], container_name)(xmlelem=xmlelem,inpline=inpline)
    except (AttributeError,KeyError):
      common.Error('Error when creating container {0} instance.'.format(container_name))
    else:
      return inst

  def DelayedInit(self, inp, numstep, crank):
    """Performs delayed initialization of this object: 
       replaces this object by the specified previously outputted object.
       Leaves gracefully if this object does not contain any delayed init. information.
       It is an error if this object contains delayed init. information but the delayed init. failed.

       Arguments:
       inp - the program/process input asking this object's delayed initialization
       numstep - order number of the step asking this object's delayed initialization
       crank - the main crank process
    """
    if self.delayed_init:
      step=self.delayed_init.pop('step')
      if step>=numstep:
        common.Error('Delayed init of {0} in step {1} failed: requested step {2} has not happened yet'.format(
          self.nick, numstep, step))
      if not crank.run.processes[step].out.Get(**self.delayed_init):
        common.Error('Delayed initialization of {0} in step {1} failed: no such object from step {2}.'.format(
          self.nick, numstep, step))
      inp.AddCopy(crank.run.processes[step].out.Get(**self.delayed_init), propagate=False)
      self.delayed_init={}
      inp.Delete(self)

  def XMLElemInit(self,xelem):
    """Initialize from inputted XML element"""
    for atr,val in xelem.attrib.items():
      try:
        getattr(self,atr)
      except AttributeError:
        common.Error('Unknown attribute {0} for data container {1}'.format(atr,self.nick))
      else:
        setattr(self,atr,val)
    attr_set_by_crank=[]
    tags = list(self.supported_attributes) + [ "delayed_init","attr_set_by_crank" ] + self.col_list
    for xchild in xelem:
      if xchild.tag in tags:
        if xchild.tag=='attr_set_by_crank':
          attr_set_by_crank = common.AutoConvert(xchild.text)
        elif xchild.tag=='file':
          filetyp=None
          if 'typ' in xchild.attrib:
            filetyp=xchild.attrib['typ']
          self.AddFile(xchild.text,filetyp)
        elif xchild.tag=='atomtype':
          fp,fpp,d_name,guessed_dtype=None,None,None,None
          for xatomchild in xchild:
            if xatomchild.tag=='fp':
              fp=common.AutoConvert(xatomchild.text)
            elif xatomchild.tag=='fpp':
              fpp=common.AutoConvert(xatomchild.text)
            elif xatomchild.tag=='d_name':
              d_name=xatomchild.text
            elif xatomchild.tag=='guessed_dtype':
              guessed_dtype=xatomchild.text
            else:
              common.Error('Wrong subtag {0} of atomtype when initiating data container {1} from XML'.format(
                xatomchild.tag,self.nick))
          self.SetAtomType(xchild.text.strip(),fp,fpp,d_name,guessed_dtype)
        elif xchild.tag=='xname':
          self.SetCrystalName(xchild.text)
        elif xchild.tag=='dname':
          self.SetDataName(xchild.text)
        elif xchild.tag in ('solvent_content','monomers_asym','seq_monomers','residues_mon','exp_num_atoms','custom','wavel','spgr'):
          setattr(self,xchild.tag,common.AutoConvert(xchild.text))
        elif xchild.tag=='cell':
          self.SetCell(common.AutoConvert(xchild.text))
        elif xchild.tag=='resol':
          self.resol=common.AutoConvert(xchild.text)
        elif xchild.tag=='delayed_init':
          for xdelinichild in xchild:
            self.delayed_init[xdelinichild.tag] = common.AutoConvert(xdelinichild.text)
        elif xchild.tag in self.col_list:
          self.SetLabel(xchild.tag,xchild.text,error_on_dupl=True)
        if xchild.tag in self.supported_attributes and xchild.tag not in attr_set_by_crank and \
           xchild.tag not in self.set_by_user:
          self.set_by_user.append(xchild.tag)
      else:
        common.Error('Wrong tag {0} when initializing data container {1} from XML'.format(xchild.tag,self.nick))
    if 'typ' not in attr_set_by_crank and 'typ' not in self.set_by_user:
      self.set_by_user.append('typ')


  def InputElemInit(self,input_line):
    """Initialize from input line"""
    ad={'atomtype':None,'fp':None,'fpp':None,'d_name':self.dname}
    keys = list(self.supported_attributes) + [ "obj_from", ] + self.col_list
    while input_line:
      token=input_line.pop(0)
      if token in self.types or token=='anom':
        self.SetType(token)
        if not 'typ' in self.set_by_user:
          self.set_by_user.append('typ')
      else:
        key,sep,val = token.partition('=')
        if key not in keys:
          input_line.insert(0,token)
          break
        if key=='typ':
          self.SetType(val)
        elif key =='xname':
          self.SetCrystalName(val)
        elif key =='dname':
          self.SetDataName(val)
        elif key in ('solvent_content','monomers_asym','seq_monomers','residues_mon','exp_num_atoms','wavel','spgr'):
          setattr(self,key,common.AutoConvert(val))
        elif key =='custom':
          if isinstance(common.AutoConvert(val),(list,tuple)):
            self.custom.extend(common.AutoConvert(val))
          else:
            self.custom.append(common.AutoConvert(val))
        elif key =='cell':
          self.SetCell(val)
        elif key =='resol':
          self.resol=val
        # filetype is always determined from suffix as of now!
        elif key=='file':
          filetyp=val.split('.')[-1].lower()
          if not os.path.isabs(val):
            val = os.path.join(os.getcwd(), val)
          self.AddFile(val,filetyp,adjust_filetype=True)
        elif key in ('atomtype','fp','fpp','d_name'):
          ad[key]=val
          if ad['atomtype']:
            self.SetAtomType(ad['atomtype'],ad['fp'],ad['fpp'],ad['d_name'])
            if ad['fp'] and ad['fpp']:
              ad['fp'], ad['fpp'] = None, None
        elif key in self.col_list:
          self.SetLabel(key,val,error_on_dupl=True)
        elif key=='obj_from':
          args = val.split(',')
          if not args or type(common.AutoConvert(args[0]))!=int:
            common.Error("Input keyword {0} needs an integer specifying the step".format(key))
          self.delayed_init['step']=common.AutoConvert(args[0])
          for cond in args[1:]:
            k,s,v=cond.partition('=')
            if not v:
              common.Error("Condition {0} of input keyword {1} malformed.".format(cond,key))
            import inout
            try:
              if k not in getattr(inout.input_output,'GetAll').func_code.co_varnames:
                common.Error("Wrong subkeyword {0} of input keyword {1}.".format(k,key))
            except AttributeError:
              # python 3
              if k not in getattr(inout.input_output,'GetAll').__code__.co_varnames:
                common.Error("Wrong subkeyword {0} of input keyword {1}.".format(k,key))
            self.delayed_init[k]=common.AutoConvert(v)
        if not key in self.set_by_user:
          self.set_by_user.append(key)

  def Data2XML(self, ET_element):
    """Convert the actual container object into XML"""
    xroot = ET.SubElement(ET_element, self.__class__.__name__)
    xroot.attrib['typ']=self.typ
    attr_set_by_crank = []
    for attr in self.supported_attributes:
      if hasattr(self,attr) and getattr(self,attr) and attr not in self.set_by_user:
        attr_set_by_crank.append(attr)
    if attr_set_by_crank:
      ET.SubElement(xroot,"attr_set_by_crank").text = str(attr_set_by_crank)
    if self.dname != self.default_unknown:
      ET.SubElement(xroot,'dname').text = self.dname
    if self.xname != self.default_unknown:
      ET.SubElement(xroot,'xname').text = self.xname
    for attr in ('solvent_content','monomers_asym','seq_monomers','residues_mon','custom','exp_num_atoms','wavel','spgr','cell','resol'):
      if hasattr(self,attr) and getattr(self,attr):
        ET.SubElement(xroot,attr).text = str(getattr(self,attr))
    for f in self.files:
      f.Data2XML(xroot)
    for col,lab in self.GetAllLabels():
      ET.SubElement(xroot,col).text = lab
    if 'atomtypes' in self.__dict__:
      atomtypes = [self.atomtype1,] + [ at  for at in self.atomtypes  if at!=self.atomtype1 ]
      for at in [at for at in atomtypes if at is not None]:
        for v in self.atomtypes[at]:
          if v:
            xat = ET.SubElement(xroot,'atomtype')
            xat.text = at
            ET.SubElement(xat,'fp').text=str(v[0])
            ET.SubElement(xat,'fpp').text=str(v[1])
            ET.SubElement(xat,'d_name').text=str(v[2])
            ET.SubElement(xat,'guessed_dtype').text=str(v[3])
    if self.delayed_init:
      xdelini = ET.SubElement(xroot,'delayed_init')
      for k,v in self.delayed_init.items():
        ET.SubElement(xdelini,k).text=str(v)

  def IsInputtedAttr(self, attr):
    """Returns True if the data object attribute has been inputted by user"""
    if attr in self.set_by_user:
      return True
    return False

  def GetFiles(self, filetype=None, conv_parent=None, inp_cont=None, conv_opts=None):
    """Returns list of file objects of the specified file type associated with this container

       Arguments:
       filetype - type of the files to be returned; if None or omitted then all files are returned; 
                  a vector/list of string filetypes is also accepted meaning logical OR
       conv_parent - specifies that conversion of file should be tried and parent process or program that asks it 
                     upon a first successful conversion the converted file is returned
       inp_cont - list of input containers (other than this container) for conversion (single container object also accepted)
    """
    if common.is_string(filetype):
      filetype=[filetype,]
    if filetype is None:
      result = self.files
    else:
      result = [ f  for ft in reversed(filetype)  for f in self.files   if ft is None or f.typ==ft ]
    if not result and conv_parent:
      inp_cont_loc=copy.copy(inp_cont)
      if inp_cont_loc is None:
        inp_cont_loc=[]
      try:
        inp_cont_loc.append(self)
      except AttributeError:
        inp_cont_loc=[inp_cont,self]
        if hasattr(sys,'exc_clear'): sys.exc_clear()
      inp_cont_loc=list(filter(None,inp_cont_loc))
      if conv_opts is None:
        conv_opts=[]
      converted=False
      for f in self.files:
        for ftyp in filetype:
          if f.Convert(ftyp,conv_parent,inp_cont_loc,conv_opts):
            converted=True
            break
        if converted:  break
      result = [ f  for f in self.files  if f.typ in filetype ]
    return result

  def GetFile(self, filetype=None, conv_parent=None, inp_cont=None, conv_opts=None, ind=-1):
    """Returns container's file object with the specified file type; None if there is none;
       in case of multiple files of the file type, 'ind' specifies the index (last by default)
       The same arguments as GetFiles() + the ind parameter
    """
    if inp_cont is None:
      inp_cont=[]
    if conv_opts is None:
      conv_opts=[]
    try:
      return self.GetFiles(filetype,conv_parent,inp_cont,conv_opts)[ind]
    except IndexError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      return None

  def GetFileTypes(self):
    """Returns the actually existing filetypes of the object"""
    return [f.typ for f in self.GetFiles()]

  def GetFileName(self, filetype=None, conv_parent=None, inp_cont=None, conv_opts=None, ind=-1, trim_path=False):
    """Returns filename of a file object of the specified file type; None if there is none
       In case of multiple files of the file type, 'ind' specifies the index (last by default)
       The same arguments as GetFiles() + the ind parameter
    """
    if conv_opts is None:
      conv_opts=[]
    try:
      f=self.GetFiles(filetype, conv_parent, inp_cont, conv_opts)[ind]
    except IndexError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      return None
    else:
      if f.typ=='mtz':
        self.mtzfile_asked=True
      fname = f.name
      if trim_path:
        fname = os.path.basename(fname)
      return fname

  def SetLabel(self, col, lbl=None, ignore_prefix=False, bad_lbl_obj=[], error_on_dupl=False):
    """Set (mtz) column(s) to label(s)

       Arguments:
       col - column(s) whose labels should be set; can be either a string or a list of strings
       lbl - label(s) to be set, needs to be of the same type as col (and the same order assumed for list)
             if lbl is None or omitted then the default label is set (if defined, otherwise it is an error)
       bad_lbl_obj - if the default or lbl inputted label is used by an object from the excl_lbl_obj list,
                     the smallest positive integer will be appended to the label such that it is not in bad_lbl_obj
       error_on_dupl - if true then it is an error if there are duplicate labels set for the data object
    """
    # if string inputted then convert to list
    if common.is_string(col):
      if lbl and not common.is_string(lbl):
        common.Error('Wrong argument types when setting label for {0}'.format(self.nick))
      col=[col,]
      if lbl:
        lbl=[lbl,]
    # check the inputted col
    for c in col:
      if c not in self.col_list:
        common.Error('No column list or no column {0} defined for {1}.'.format(c,self.nick))
    # assign defaults if asked for
    if not lbl:
      if not self.label_default or self.typ not in self.label_default:
        common.Error('No default label definition for {0} of type {1}'.format(self.nick, self.typ))
      lbl = [ self.label_default[self.typ][self.col_list.index(c)]  for c in col ]
    # do the actual setting
    for c,l in zip(col,lbl):
      if not ignore_prefix:
        l=self.labelout_prefix+l
      # check whether the inputted/default lbl do not collide with existing labels
      if l in [self.label[cl] for cl in self.col_list if cl!=c]:
        issue = "Label {0} to be assigned for column {1} already exists for another column in {2}".format(l,c,self.GetFileName())
        if error_on_dupl:
          common.Error(issue)
        else:
          common.Warning(issue)
      self.label[c]=l
    # make sure the labels are unique if requested and reset them
    if bad_lbl_obj:
      from program import program
      lbl_o,lbl_u = program(None).GetUniqueLabels( list(bad_lbl_obj)+[self,], allow_multiple=True )
      lbl = [lbl_u[self.GetFileName('mtz')][::-1][lbl_o[self.GetFileName('mtz')][::-1].index(l)]  for l in lbl]
      for c,l in zip(col,lbl):
        self.label[c]=l

  def GetLabel(self, col, mtz_not_asked=False, ignore_undef_col=False):
    """Get label of the inputted column"""
    if col not in self.col_list: 
      if not ignore_undef_col:
        common.Error('Column {0} not defined for container {1}'.format(col,self.nick))
      else:
        return None
    if not mtz_not_asked:
      self.mtzfile_asked=True
    return self.label[col]

  def InitLabels(self, lab_dict={}):
    """Initialize (mtz) label dictionary with keys from internal col_list and (optionally) inputted dictionary values"""
    self.label={}
    for col in self.col_list:
      self.label[col]=None
    for col in lab_dict:
      self.SetLabel(col,lab_dict[col])

  def GetAllLabels(self, labels_only=False, mtz_not_asked=False):
    """Returns all defined (non None) labels in a list of pairs (column,label)
       Returns just a list of (non None) labels if labels_only is True (False by default)
    """
    if not labels_only:
      return [ (c,self.GetLabel(c,mtz_not_asked))  for c in self.col_list  if self.GetLabel(c) is not None ]
    else:
      return [ self.GetLabel(c,mtz_not_asked)  for c in self.col_list  if self.GetLabel(c) is not None ]

  def GetFullLabel(self, *cols, **kwargs):
    """Return full (mtz) label path for specified column(s) (ie 'crystal/data/[labels]')"""
    objs = [self,] + list(kwargs.get('other_objs', []))
    if not self.col_list or not self.GetFile('mtz'):
      common.Error('No column list or no mtz file for {0}.'.format(self.nick))
    crystal_name = self.xname if self.xname!=self.default_unknown  else '*'
    data_name = self.dname if self.dname!=self.default_unknown  else '*'
    labels=[]
    for obj in objs:
      if obj.xname!=self.xname:
        common.Warning('xname mismatch when making full label: {0} vs {1}'.format(obj.xname,self.xname))
      if obj.dname!=self.dname:
        common.Warning('dname mismatch when making full label: {0} vs {1}'.format(obj.dname,self.dname))
      for col in cols:
        labels.append(obj.GetLabel(col))
    return '/'.join((crystal_name,data_name,'['+','.join(labels)+']'))

  def SetType(self, typ):
    """Set type of the data container"""
    if typ in self.types:
      self.typ=typ
    elif typ=='anom' and self.nick=='fsigf':
      # allowing fsigf 'anom' type = plus+minus, simplifying non-mtz anomalous inputs (eg shelx)
      self.typ='plus'
      self.minus_needs_to_be_added=True
    else:
      common.Error('Type {0} not allowed for container {1} (the following types are allowed: {2})'.format(
        typ,self.nick,','.join(self.types)))

  def GetType(self):
    """Get type of the data container"""
    return self.typ

  def AddFile(self, filename, filetype=None, adjust_filetype=False):
    """Add/associate file with filename of filetype to the data container
       The filename must be supplied as absolute path otherwise there is an error
    """
    assert os.path.isabs(filename), 'File {0} is not specified as absolute path.'.format(filename)
    if not filetype:
      filetype=self.deffiletype
    elif adjust_filetype:
      # detecting hkl vs. HKL since it can be easily switched, esp. with hkl on/from Windows
      if filetype=='hkl' and os.path.isfile(filename):
        with open(filename) as f:
          if 'FORMAT=XDS_ASCII' in f.readline():  filetype='HKL'
          else:                                   filetype='hkl'
    if filetype not in self.supported_filetypes:
      common.Error('Wrong filetype {0} for container {1}'.format(filetype,self.nick))
    self.files.append(files(filename,filetype))

  def SetFile(self, filename, filetype=None):
    """The same as AddFile() but any files of the same filetype will be removed before adding the file"""
    if filetype:
      filetype=filetype
    try:
      self.DetachFile(filetype=filetype)
    except self.NoDetachFileError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
    self.AddFile(filename,filetype)

  class NoDetachFileError(Exception):
    pass

  def DetachFile(self, filename=None, filetype=None):
    """Removes association of file with filename and/or file(s) of filetype with the container
       It is an error if there is no such file
    """
    if not filetype:
      filetype=self.deffiletype
    removed=False
    for f in self.files:
      if f.typ==filetype or f.name==filename:
        self.files.remove(f)
        removed=True
    if not removed:
      self.NoDetachFileError('No files with name {0} and/or type {1} could be detached from container {2}'.format(
                   filename, filetype, self.nick) )

  def SetCrystalName(self, name):
    if name:
      self.xname=str(name)

  def GetCrystalName(self):
    return self.xname

  def SetDataName(self, name):
    if name:
      self.dname=str(name)

  def GetDataName(self):
    return self.dname

  def SetWavel(self, wavel):
    if wavel:
      try:
        wavel = float(wavel)
      except ValueError:
        common.Warning('Wavelength must be a real number, "{0}" cannot be accepted.'.format(wavel))
      self.wavel=wavel
      return self.wavel

  def SetCell(self, cell):
    if cell:
      if common.is_string(cell):
        cell=cell.split(',')
      if len(cell)<3:
        common.Warning('At least 3 cell parameters required, only {0} supplied: {1}'.format(len(cell),*cell))
        cell=[]
      try:
        cell = [float(c) for c in cell]
      except ValueError:
        common.Warning('Cell parameters must be real numbers, {0} cannot be accepted.'.format(*cell))
        cell=[]
      self.cell=cell
      return self.cell

  def GetCellSpacegroupResol(self, pro, get='', accept_none=False, non_mtz=None):
    # pro is here to define the directory (associated to the program/process pro) in which sftools will be run
    # furthermore, if the information cannot be obtained from this object, it will attempt to get it 
    # from another object from pro.inp with the same xname, dname if possible
    info=None
    mtz_obj=None
    if not non_mtz:
      if self.GetFile('mtz'):
        mtz_obj = self
      else:
        mtz_obj = pro.inp.Get('fsigf',filetype='mtz',xname=self.GetCrystalName(),dname=self.GetDataName(),try_convert=False)
    # non mtz retrieval attempt
    if not mtz_obj:
      hkl_sca = self.GetFile(['HKL','sca'])
      if hkl_sca and get in ('cell','spg','spgn','wavel'):
        import re
        with open(hkl_sca.name,'r') as f:
          for i,line in enumerate(f):
            if hkl_sca.typ=='sca' and i==2:
              cell_space = re.match('\s*(\d+\.\d+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)',line)
              if cell_space and len(cell_space.groups())>6:
                if get=='cell':  self.SetCell( cell_space.groups()[:6] )
                if get=='spg':  self.spgr = cell_space.groups()[6]
              break
            if hkl_sca.typ=='HKL':
              if re.match('!END_OF_HEADER',line):
                break
              if re.match('\s*-?\d+\s+-?\d+\s+-?\d+\s+-?\d\.\d+\s+\d\.\d+',line):
                break
              if get=='cell':
                cell = re.match('!UNIT_CELL_CONSTANTS=\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)\s*(\S+)',line)
                if cell:
                  if len(cell.groups())>5:
                    self.SetCell( cell.groups()[:6] )
                  break
              if get=='spgn':
                spg_num = re.match('!SPACE_GROUP_NUMBER=\s*(\d+)',line)
                if spg_num:
                  self.spgr_num=spg_num.group(1)
                  break
              if get=='wavel':
                wavel = re.match('!X-RAY_WAVELENGTH=\s*(\d*\.\d+)',line)
                if wavel:
                  self.SetWavel( wavel.group(1) )
                  break
        if get=='cell': info=self.cell
        if get=='spg':  info=self.spgr
        if get=='spgn': info=self.spgr_num
        if get=='wavel': info=self.wavel
      if hkl_sca and get=='res':
        # this is a bit of a hack.  A better solution should be found...
        pro.BackupAnyPars()
        pro.SetParam('target', 'RICE')
        shelxc=pro.AddProg('shelxc',propagate_out=False)
        shelxc.inp.Set(self)
        shelxc.SetRunDir('check')
        shelxc.Run()
        pro.RestoreAnyPars()
        resol_out=shelxc.GetStat('resol',accept_none=accept_none)
        if resol_out: 
          self.resol=resol_out[-1]
          info=self.resol
        pro.programs.remove(shelxc)
      if info is None:
        if accept_none:
          return None
        else:
          getname,hklname=get,'any'
          if get=='spg': getname='spacegroup'
          if get=='spgn': getname='spacegroup number'
          if hkl_sca:  hklname=hkl_sca.name
          common.Error('Could not retrieve {0} from {1} file'.format(getname,hklname))
    # mtz retrieval
    else:
      if get=='res' and not mtz_obj.GetLabel('f',ignore_undef_col=True) and not mtz_obj.GetLabel('i',ignore_undef_col=True):
        mtz_obj = pro.inp.Get(filetype='mtz',xname=self.GetCrystalName(),dname=self.GetDataName(),label='f')
        if not mtz_obj:
          mtz_obj = pro.inp.Get(filetype='mtz',xname=self.GetCrystalName(),dname=self.GetDataName(),label='i')
        if not mtz_obj:
          if accept_none:
            return None
          common.Error('Resolution could not be retrieved.')
      try:
        pro=pro.process
      except AttributeError:
        if hasattr(sys,'exc_clear'): sys.exc_clear()
      while pro.rundir is None:
        pro=pro.parent_process
      sft = pro.AddProg('sftools')
      sft.runname='sftools_getinfo'
      sft.SetKey('read', '"'+mtz_obj.GetFileName('mtz')+'"')
      sft.SetKey('Y')  # for the rare case of interactive question, typically irrelevant about xplor free ref...
      if get=='res':
        if mtz_obj.GetLabel('f'):
          sft.SetKey('select',['col',mtz_obj.GetLabel('f')])
        else:
          sft.SetKey('select',['col',mtz_obj.GetLabel('i')])
        sft.SetKey('checkhkl')
      elif get=='wavel':
        sft.SetKey('list dwave')
      else:
        sft.SetKey('list')
        if self.xname!=self.default_unknown:
          sft.SetKey('list dcell')
      sft.SetKey('quit','\nY')
      try:
        sft.Run()
      except:
        pro.programs.remove(sft)
        raise
      if get=='cell':
        self.SetCell( ( sft.GetStat('a',self.xname), sft.GetStat('b',self.xname), \
                        sft.GetStat('c',self.xname), sft.GetStat('alpha',self.xname), \
                        sft.GetStat('beta',self.xname), sft.GetStat('gamma',self.xname) ) )
        info = self.cell
      elif get=='spg':
        self.spgr=sft.GetStat('spacegroup').strip()
        info = self.spgr
      elif get=='spgn':
        self.spgr_num=sft.GetStat('spacegroup_num',accept_none=accept_none)
        info = self.spgr_num
      elif get=='res':
        if self.GetLabel('f'):
          self.resol=sft.GetStat('resolution',self.GetLabel('f')[:12],accept_none=accept_none)
          self.resol_low=sft.GetStat('low_resol',self.GetLabel('f')[:12],accept_none=accept_none)
        else:
          self.resol=sft.GetStat('resolution',self.GetLabel('i')[:12],accept_none=accept_none)
          self.resol_low=sft.GetStat('low_resol',self.GetLabel('i')[:12],accept_none=accept_none)
        info = self.resol
      elif get=='wavel':
        xname = self.xname # if self.xname!=self.default_unknown  else r'\S+'
        dname = self.dname # if self.dname!=self.default_unknown  else r'\S+'
        self.wavel = sft.GetStat( 'wavelength', (xname,dname), accept_none=True, param_escape=False )
        info = self.wavel
      else:
        common.Error('Internal error in GetCellSpacegroup.')
      pro.programs.remove(sft)
      if not info and non_mtz is None:
        info=self.GetCellSpacegroupResol(pro, get, accept_none, non_mtz=True)
    return info

  def GetWavelength(self,pro,accept_none=False):
    if self.wavel:
      return self.wavel
    else:
      return self.GetCellSpacegroupResol(pro,get='wavel',accept_none=accept_none)

  def GetCell(self,pro,accept_none=False):
    if self.cell:
      return self.cell
    else:
      return self.GetCellSpacegroupResol(pro,get='cell',accept_none=accept_none)

  def GetSpacegroup(self,pro,accept_none=False):
    if self.spgr:
      return self.spgr
    else:
      return self.GetCellSpacegroupResol(pro,get='spg',accept_none=accept_none)

  def GetSpacegroupNumber(self,pro):
    if self.spgr_num:
      return self.spgr_num
    else:
      return self.GetCellSpacegroupResol(pro,get='spgn')

  def GetResolution(self,pro,accept_none=False,lowres=False):
    if lowres:
      if self.resol_low:
        return self.resol_low
      else:
        return self.GetCellSpacegroupResol(pro,get='res',accept_none=accept_none)
    else:
      if self.resol:
        return self.resol
      else:
        return self.GetCellSpacegroupResol(pro,get='res',accept_none=accept_none)

  def AddCustomTag(self,tag):
    self.custom.append(tag)


class fsigf(data_container):
  """Observed amplitudes/intensities container class"""
  description="observed amplitudes/intensities"
  supported_filetypes=('mtz','hkl','HKL','sca','drear','cif')
  deffiletype='mtz'
  col_list=['f','sigf','i','sigi','e','sige','alpha']
  types=('average','plus','minus','fa','delta-anom')
  typ='average'
  label_default = { 'average': ['F','SIGF','I','SIGI','E','SIGE',None], \
                    'plus': ['F+','SIGF+','I+','SIGI+','E+','SIGE+',None], \
                    'minus': ['F-','SIGF-','I-','SIGI-','E-','SIGE-',None], \
                    'fa': ['FA','SIGFA',None,None,'EA','SIGEA','ALPHA'], \
#                    'ea': ['EA','SIGEA',None,None,'EA','SIGEA','ALPHA'], \
                    'delta-anom': ['DA','SIGDA','IDA','SIGIDA','EDA','SIGEDA',None,], \
                    'average-derived': ['F','SIGF','I','SIGI','E','SIGE',None], \
                  }
  supported_attributes=('typ','file','xname','dname','custom','wavel','spgr','cell','resol')


class mapcoef(data_container):
  """Map coefficients (incl. phase information) container class"""
  description="map coefficients (incl. phase information)"
  supported_filetypes=('mtz','map','phs')
  deffiletype='mtz'
  col_list=['f', 'ph', 'fom', 'hla', 'hlb', 'hlc', 'hld',]
  types=('best','weighted','combined','densmod','anomalous','anom-diff','diff','mask')
  typ='best'
  label_default = { 'best': ['FB','PHIB','FOM','HLA','HLB','HLC','HLD'], \
                    'weighted': ['FWT','PHWT',None,None,None,None,None], \
                    'combined': ['FCOMB','PHCOMB','FOMCOMB','HLACOMB','HLBCOMB','HLCCOMB','HLDCOMB'], \
                    'densmod': ['FMOD','PHMOD','FOMMOD',None,None,None,None], \
                    'anomalous': ['FAN','PHAN',None,None,None,None,None], \
                    'anom-diff': ['DELFAN','PHDELAN',None,None,None,None,None], \
                    'diff': ['DELFWT','PHDELWT',None,None,None,None,None], \
                    'mask': ['FMSK','PHMSK',None,None,None,None,None], \
                  }
  supported_attributes=('typ','file','xname','dname','custom','spgr','cell','resol')


class sequence(data_container):
  """Sequence container class"""
  description="protein sequence"
  supported_filetypes=('pir','fasta','fas','seq')
  deffiletype='pir'
  types=('protein',)
  typ='protein'
  supported_attributes=('file','xname','solvent_content','monomers_asym','seq_monomers','residues_mon')
  # the sequence string, if known.
  seqstr=""

  def GetSequenceString(self):
    """Gets the sequence for the file for internal use/decision making.
       At the moment, there is no proper support for the different sequence formats implemented,
       just a simple heuristics is used to grab the sequence.
    """
    if self.seqstr:
      return self.seqstr
    self.seqstr=""
    if self.GetFileName() and os.path.isfile(self.GetFileName()):
      with open(self.GetFileName(),'r') as f:
        il,lines=0,0
        monom,im = ["",], 0
        for il,line in enumerate(f):
          line=line.strip().replace(' ','').upper()
          if not line.startswith('>') and not line.startswith('#'):
            #exclude RNA/DNA
            if len(line)<4 or (line.count('C')+line.count('G')+line.count('A')+line.count('U')+line.count('T'))/len(line)<0.9:
              self.seqstr+=line;  monom[im]+=line;  lines+=1
          else:
            monom.append("");  im+=1
          if il==1 and self.GetFile().typ=='pir' and len(line)<80:
            line1, self.seqstr = line, ""
          if lines>1000:
            common.Error('More than 1000 lines in the sequence file {}'.format(self.GetFileName()))
        if not self.seqstr and il in (1,2) and self.GetFile().typ=='pir':
          self.seqstr+=line1
        # 100% sequence identity required as of now for multiple monomers in the input sequence
        self.seq_monomers=1
        for mon in monom:
          if len(mon)>10 and monom.count(mon)>self.seq_monomers:
            self.seq_monomers = monom.count(mon)
    return self.seqstr

class model(data_container):
  """Model (of crystal content) container class"""
  description="crystal model and its description"
  # xyz is the PMF substr. output; res is the SHELXD output substr.; hat is the SHELXE susbtr. output; frac is the York fractional format
  supported_filetypes=('pdb','xyz','res','hat','frac')
  deffiletype='pdb'
  types=('partial','substr','partial+substr','unknown','patterson')
  typ='unknown'
  exp_num_atoms=None
  supported_attributes=('typ','file','xname','atomtype','fp','fpp','solvent_content','monomers_asym','residues_mon','exp_num_atoms','custom','d_name')

  def __init__(self,*args,**kwargs):
    # dictionary of atom types and their f',f''
    self.atomtypes={}
    self.atomtype1=kwargs.pop('atomtype1',None)
    for at,fs in kwargs.pop('atomtypes',{}).items():
      for f in fs:
        self.SetAtomType(at,f[0],f[1],f[2])
    if self.atomtype1 and self.atomtype1 not in self.atomtypes:
      self.SetAtomType(self.atomtype1)
    super(model,self).__init__(*args,**kwargs)

  def SetAtomType(self, atomtype, fp=None, fpp=None, dname=None, guessed_dtype=None):
    """Set atom type (typically substructure) and (optionally) its f' and f''"""
    atomtype=atomtype.upper()
    if fp is not None and fpp is not None:
      try:
        fp=float(fp)
        fpp=float(fpp)
      except:
        common.Error("f' and f'' must be float ({0},{1} inputted for atom type {2})".format(fp,fpp,atomtype))
    if not self.atomtype1:
      self.atomtype1=atomtype
    if atomtype not in self.atomtypes:
      self.atomtypes[atomtype] = []
    # only one assignment makes sense per atomtype and dname
    # we adjust existing record if its dname matches or if it is completely empty,  otherwise a new one added
    ia = [i for i,v in enumerate(self.atomtypes[atomtype]) if v[2]==dname or (v[2]==self.default_unknown and not any((v[0],v[1],v[3])))] if atomtype in self.atomtypes else []
    if ia:
      self.atomtypes[atomtype][ia[0]] = (fp,fpp,dname,guessed_dtype)
    else:
      self.atomtypes[atomtype].append((fp,fpp,dname,guessed_dtype))

  def SetAtomTypes(self,atomtypes,atomtype1=None):
    """Takes the input dictionary of atom types (and f',f'' if supplied) and sets it to this container"""
    if type(atomtypes) != dict:
      common.Error('SetAtomTypes() requires atom types dictionary as input')
    self.atomtypes=atomtypes
    if atomtype1:
      self.atomtype1=atomtype1

  def GetAtomTypes(self, getlist=False):
    """Returns dictionary of atom types (typically substructure) with a list of (f',f'',dname,guessed_type) 
       as dict. values  by default  or a list of atom types if 'getlist' is True
    """
    if getlist:
      lst = list(self.atomtypes.keys())
      if self.atomtype1 and lst[0]!=self.atomtype1:
        lst.remove(self.atomtype1)
        lst.insert(0, self.atomtype1)
      return lst
    else:
      return self.atomtypes

  def Getfpfpp(self, atomtype=None, dname=None, return_guessed_dtype=False):
    """Returns (fp,fpp,dname,atomtype) for the specified atomtype and dname
       If atomtype or dname is not specified (None) then any (first) stored atomtype and/or dname is used.
       Returns (None,None,None,None) if fp,fpp for no such atomtype+dname combination is stored.
       This function only reports fp,fpp currently stored values, ie does not try to estimate them.
    """
    at_lst, at = [], atomtype
    if atomtype in self.atomtypes:
      at_lst = self.atomtypes[atomtype]
    elif self.atomtypes.values() and atomtype is None:
      at = list(self.atomtypes.keys())[0]
      at_lst = self.atomtypes[at]
    guessed_dtype = []
    if return_guessed_dtype:
      guessed_dtype = [None,]
    if at_lst:
      dn_match = None
      if dname is None:
        dn_match = at_lst
      else:
        dn_match = [v for v in at_lst  if v[2]==dname]
      if dn_match:
        if return_guessed_dtype:
          guessed_dtype = [dn_match[0][3],]
        return [dn_match[0][0],dn_match[0][1],dn_match[0][2],at]+guessed_dtype
    return [None,None,None,None]+guessed_dtype

  def GetAtomType(self):
    """Returns the first (main) atom type """
    return self.atomtype1

  def GuessNumSubstrAtomsFromSeq(self,proces):
    seq_obj = proces.inp.Get('sequence', filetype=sequence.supported_filetypes, try_convert=False)
    at = self.GetAtomType()
    # I think this could fail if multiple models can be inputted for the same crystal 
    # (in which case this might be fixed eg by looping through all such crystals?)
    if seq_obj and at:
      if not seq_obj.seqstr:
        seq_obj.GetSequenceString()
      if at=='SE':
        self.exp_num_atoms = seq_obj.seqstr.count('M') + seq_obj.seqstr.count('U')
      elif at=='S':
        self.exp_num_atoms = seq_obj.seqstr.count('M') + seq_obj.seqstr.count('C')
      elif at=='I' or at=='BR':
        self.exp_num_atoms = max( len(seq_obj.seqstr)/20, 2)
      else:
        self.exp_num_atoms = max( len(seq_obj.seqstr)/200, 2)  # a wild guess only.
      monom_asym = self.monomers_asym if self.monomers_asym  else seq_obj.monomers_asym
      if monom_asym is None:
        matthews=proces.AddProcess('matthews', propagate_out=False)
        #matthews.no_reporting=True
        matthews.Run()
        proces.processes.remove(matthews)
        monom_asym = seq_obj.monomers_asym
      self.exp_num_atoms = self.exp_num_atoms * monom_asym


class exclude(data_container):
  """Excluded reflections container class"""
  description="'free' reflections to be excluded"
  supported_filetypes=('mtz',)
  deffiletype='mtz'
  types=('freeR','freebias')
  typ='freeR'
  col_list=['free',]
  label_default = { 'freeR': ['FREER'],  \
                    'freebias': ['BIASFREE']  }
  supported_attributes=('typ','file','xname','dname')


class datafile(data_container):
  """Container class for various other, typically program specific, data storage"""
  # type is given by filetype:
  # ins = SHELX input parameters
  # xml = .xml file
  # mtz = .mtz file (if we are not interested in its crystallographic information, eg for mtz merging)
  # dluz = Refmac/Multicomb Luzzati D parameters file
  # par = ARP/wARP .par file
  description="other data - defined by file type"
  supported_filetypes=('ins','xml','mtz','dluz','par')
  deffiletype='ins'
  types=supported_filetypes
  typ=deffiletype
  supported_attributes=('typ','file','custom')

  def __init__(self, filename=None, filetype=None, typ=None, *args, **kwargs):
    if filetype is None and typ:
      filetype=typ
    if typ is None and filetype:
      typ=filetype
    super(datafile,self).__init__(filename=filename,filetype=filetype,typ=typ,*args,**kwargs)
