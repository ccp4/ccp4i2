import os,sys,copy,re
from xml.etree import ElementTree as ET
from . import data, common



class input_output(object):
  """program or process input and output base class"""
  never_propagate=False

  def __init__(self, is_output, parent):
    # create empty list for each known data container class
    for data_class in data.data_container.__subclasses__():
      setattr(self, data_class.__name__, [])
    # program or process this input/output belongs to
    self.parent=parent
    self.is_output=is_output
    self.order_num=0
    # check whether we have multiple inputs/outputs
    if parent:
      numstr=''
      while hasattr(parent, 'out'+numstr if is_output else 'inp'+numstr):
        self.order_num+=1
        numstr=str(self.order_num+1)

  def XMLElemInit(self, xelem, dummy=False):
    """Initialize from inputted XML element"""
    tags = [ sc.__name__ for sc in data.data_container.__subclasses__()]
    for xchild in xelem:
      if xchild.tag in tags:
        self.AddCopy( data.data_container.from_name(xchild.tag,xmlelem=xchild), propagate=not dummy )
      else:
        common.Error('Wrong tag when parsing XML for input/output of {0}: {1}'.format(self.parent.nick,xchild.tag))

  def Data2XML(self, ET_element):
    """Convert the actual input/output object into XML"""
    if self.Get():
      if self.is_output:
        xroot = ET.SubElement(ET_element,'out')
      else:
        xroot = ET.SubElement(ET_element,'inp')
    for obj in self.GetAll(stored_order=True):
      obj.Data2XML(xroot)

  def DelayedInit(self, numstep, crank, recur=True):
    """Asks for delayed initialization of any objects of this input

       Arguments:
       numstep - order number of the step asking the delayed initialization
       crank - the main crank process object
       recur - optional; if True then also recursively asks delayed init for inputs of all 
               progs/processes of children of this input's parent (True by default)
    """
    from . import process
    for o in self.GetAll():
      o.DelayedInit(self, numstep, crank)
      if recur and isinstance(self.parent, process.process):
        for p in self.parent.processes:
          p.inp.DelayedInit(numstep, crank)
        for p in self.parent.programs:
          p.inp.DelayedInit(numstep, crank)

  def AddNew(self, container_name, *args, **kwargs):
    """Creates new data container object and adds it to program input/output.

       Arguments:
         container_name : data container name (eg 'fsigf','mapcoef','model') - compulsory
         filename : name of file to be associated with the new object
         filetype : type of file to be associated with the new object
         typ :      type of the new object
         lab_dict : dictionary of mtz column labels
         propagate: the container will be also added to all child processes/programs (if self 
                    is output class) or parent processes (if self is input class) (True by default)

       Returns the created container object.
    """
    # remove propagate and checkinpfile from kwargs and sets it True by default
    propagate = kwargs.pop('propagate', True)
    checkinpfile = kwargs.pop('checkinpfile', True)
    # append the object to the current list of objects
    if container_name in self.__dict__:
      # call addfiletochild in order to make sure we have absolute path for filename
      filename, filetype = kwargs.pop('filename', None), kwargs.pop('filetype', None)
      if args:
        args=list(args)
        filename=args.pop(0)
      if args:
        filetype=args.pop(0)
      getattr(self, container_name).append( \
        getattr(data,container_name)(*args, **kwargs) )
      if filename:
        self.AddFileToChild(getattr(self, container_name)[-1],filename,filetype)
    else:
      common.Error('Wrong data class name: {0} - object cannot be created.'.format(container_name))
    if self.is_output:
      self.SetContainerPrefix(getattr(self, container_name)[-1])
    if propagate:
      self.Propagate(getattr(self, container_name)[-1])
    if checkinpfile:
      self.CheckInpFile(getattr(self, container_name)[-1])
    return getattr(self, container_name)[-1]

  def AddCopy(self, container, propagate=True, allow_duplicates=False, checkinpfile=True):
    """Adds (deep)copy of the inputted data container object to program input/output, returns the copy
       propagate can be True (resursive propagation - default), False (no propagation) or 'shallow' (one level propagation only)
    """
    try:
      while container in getattr(self, container.nick) and not allow_duplicates:
        getattr(self, container.nick).remove(container)
      getattr(self, container.nick).append( copy.deepcopy(container) )
    except AttributeError:
      common.Error('Wrong type of container {0} to be copied for input/output of {1}.'.format(type(container),self.parent.nick))
    if self.is_output:
      self.SetContainerPrefix(container)
    if propagate:
      self.Propagate(getattr(self, container.nick)[-1], recursive=(propagate!='shallow'))
    if checkinpfile:
      self.CheckInpFile(getattr(self, container.nick)[-1])
    return getattr(self, container.nick)[-1]

  def Add(self, container, propagate=True, allow_duplicates=False, ind=None):
    """Adds the inputted data container object to program input/output

       Arguments:
       container - the container object to be added
       propagate - if True then the container is propagated up (for output) or down (for input) the tree
       allow_duplicates - if True then the same object is allowed in the input/output multiple times
                          if False then if the same object existed in the input/output it is deleted
    """
    try:
      while container in getattr(self, container.nick) and not allow_duplicates:
        getattr(self, container.nick).remove(container)
      if ind is None or len(getattr(self, container.nick))<=ind:
        getattr(self, container.nick).append(container)
      else:
        getattr(self, container.nick).insert(ind,container)
    except AttributeError:
      common.Error('Wrong type of container: {0} to be added for input/output of {1}.'.format(type(container),self.parent.nick))
    if self.is_output:
      self.SetContainerPrefix(container)
    if propagate:
      self.Propagate(container,ind)
    return container
    #return getattr(self, container.nick)[-1]

  def CheckInpFile(self,container):
    if not self.is_output:
      for f in container.GetFiles():
        if not os.path.isfile(f.name):
          common.Error('Input file {0} could not be found!'.format(f.name))

  def AddFileToChild(self, child, filename, filetype=None, no_rewrite=None):
    """Add file to the specified container belonging to this input/output.
       If the file name is supplied as relative path then the program/process's rundir is prepended

       Arguments:
       child - container object belonging to this input/output
       filename - name of the file
       filetype - type of the file (optional - default filetype of container used otherwise)
    """
    assert hasattr(child,'nick'), \
      'Cannot add file to object {0}: wrong object.'.format(child)
    assert child in getattr(self,child.nick), \
      'Cannot add file to object {0}: not child of {1}'.format(child,self)
    if not os.path.isabs(filename):
      if self.parent.rundir is not None:
        filename=os.path.join(self.parent.rundir, filename)
      else:
        filename=os.path.join(os.getcwd(), filename)
    if no_rewrite or (no_rewrite is None and self.parent.GetParam('no_rewrite')):
      while os.path.isfile(filename):
        base,dot,suffix=filename.rpartition('.')
        num=0
        re_res=re.search(r'_(\d+)$',base)
        if re_res:
          base=base[:len(base)-len(re_res.group(0))]
          num=int(re_res.group(1))
        filename=base+'_'+str(num+1)+dot+suffix
    child.AddFile(filename,filetype)

  def SetFileToChild(self, child, filename, filetype=None, no_rewrite=False):
    """The same as AddFileToChild() but removing any filenames of the same type before adding"""
    child.DetachFile(filetype=filetype)
    self.AddFileToChild(child,filename,filetype,no_rewrite=no_rewrite)

  def FindPropagateObjects(self):
    objects=[]
    if self.is_output:
      # find parent processes of the program/process this output belongs to
      try:
        objects.append(self.parent.process) 
      except AttributeError:
        if hasattr(sys,'exc_clear'): sys.exc_clear()
      try:
        if self.parent.parent_process is not None:
          objects.append(self.parent.parent_process)
      except AttributeError:
        if hasattr(sys,'exc_clear'): sys.exc_clear()
    else:
      # if this input belongs to a process, find its children (programs have no children)
      try:
        objects.extend(self.parent.programs) 
        objects.extend(self.parent.processes)
      except AttributeError:
        if hasattr(sys,'exc_clear'): sys.exc_clear()
    return objects

  def Propagate(self,container,ind=None,recursive=True):
    # propages the container to parent processes (for output) or child processes/programs (for input)
    if self.never_propagate:
      return
    objects=self.FindPropagateObjects()
    for o in objects:
      inpout = o.out if self.is_output else o.inp
      if hasattr(self,'order_num') and self.order_num:
        if self.is_output and not hasattr(o,'out'+str(self.order_num+1)):
          setattr(o,'out'+str(self.order_num+1), input_output(is_output=True,parent=o))
        elif not self.is_output and not hasattr(o,'inp'+str(self.order_num+1)):
          setattr(o,'inp'+str(self.order_num+1), input_output(is_output=False,parent=o))
        inpout = getattr(o,'out'+str(self.order_num+1)) if self.is_output else getattr(o,'inp'+str(self.order_num+1))
      inpout.Add(container,allow_duplicates=False,ind=ind,propagate=recursive)

  def SetContainerPrefix(self, container):
    # set the program output label prefix (if defined)
    if self.parent and hasattr(self.parent,'labelout_prefix'):
      container.labelout_prefix = self.parent.labelout_prefix

  def GetCopy(self, *args, **kwargs):
    """Filter that returns container with specified properties
       The same as Get() but returns the (deep)copy of the object
    """
    return copy.deepcopy( self.Get(*args, **kwargs) )

  def Get(self, *args, **kwargs):
    """Filter that returns container with specified properties
       if there is no such container then None is returned
       if there are multiple such containers then the last stored is returned by default unless 
       'ind' parameter is passed specifying the order number of the container to be passed
    """
    ind = kwargs.pop('ind', -1)
    getall=self.GetAll(stored_order=True,last=(ind==-1),*args,**kwargs)
    if getall:
      try:
        return getall[ind]
      except (TypeError,IndexError):
        if hasattr(sys,'exc_clear'): sys.exc_clear()
    return None
    #  return self.GetAll(stored_order=True,last=(ind==-1),*args,**kwargs)[ind]
    #except (TypeError,IndexError):
    #  if hasattr(sys,'exc_clear'): sys.exc_clear()
    #  return None

  def GetAll(self, container_name=None, typ=None, col=None, label=None, filetype=None, filename=None, custom=None, is_native=None,
                   has_atomtypes=None, has_solvent_content=None, has_monomers_asym=None, has_num_atoms=None, has_residues_mon=None,
                   has_wavel=None, has_fpp=None, xname=None, dname=None, stored_order=False, try_convert=True, inp_cont=None, conv_opts=None, 
                   last=None, exclude_cont=None):
    """Filter that returns a list of all containers of this input_output with specified properties
       Returns empty list if there are none
       The returned list is in reverse order (to the stored order, ie LIFO) unless 'stored_order' is True
       typ and filetype and custom can be iterables in which case disjunction filter (or) is assumed
       col and label can be iterable in which case conjuction filter (and) is assumed
    """
    # prevent indefinite convert recursion
    if try_convert and self.parent.nick in __import__('processes.convert').convert.convert.supported_progs:
      try_convert=False
    if try_convert:
      try_convert=self.parent
    if has_fpp:
      has_atomtypes=True
    # convert strings to lists for variables that are (in general) iterables
    if common.is_string(col):
      col=[col,]
    if common.is_string(label):
      label=[label,]
    if common.is_string(typ):
      typ=[typ,]
    if common.is_string(filetype):
      filetype=[filetype,]
    if common.is_string(custom):
      custom=[custom,]
    if inp_cont is None:
      inp_cont=[]
    if conv_opts is None:
      conv_opts=[]
    if isinstance(inp_cont,data.data_container):
      inp_cont=[inp_cont,]
    if isinstance(exclude_cont,data.data_container):
      exclude_cont=[exclude_cont,]
    # do the actual filtering
    filtered=[]
    finished=False
    if not container_name:
      container_names = [ dc.__name__ for dc in data.data_container.__subclasses__() ]
    else:
      container_names = [container_name,]
    for container_name in container_names:
      if container_name in self.__dict__:
        for c in reversed(getattr(self,container_name)):
          # input plus/minus partners together for conversion (needed due to cad crashes if separate...)
          if filetype and try_convert:
            inp_cont_orig=inp_cont[:]
            if c.GetType()=='plus' and not [ic for ic in inp_cont if ic.GetType()=='minus']:
              inp_cont.append( self.Get(typ='minus',filetype=c.GetFileTypes(),xname=c.xname,dname=c.dname,try_convert=False) )
            if c.GetType()=='minus' and not [ic for ic in inp_cont if ic.GetType()=='plus']:
              inp_cont.append( self.Get(typ='plus',filetype=c.GetFileTypes(),xname=c.xname,dname=c.dname,try_convert=False) )
          if ( (typ is None or any(c.GetType()==t for t in typ))                               and \
               (xname is None or c.GetCrystalName()==xname)                                    and \
               (dname is None or c.GetDataName()==dname)                                       and \
               (custom is None or any(cus in c.custom for cus in custom))                      and \
               (is_native is None or (is_native==('native' in c.custom or c.GetCrystalName()=='native')))  and \
               (has_atomtypes is None or (hasattr(c,'atomtypes') and c.atomtypes))             and \
               (has_fpp is None or [pp for (p,pp,d,t) in c.GetAtomTypes()[c.GetAtomType()] if pp]) and \
               (has_solvent_content is None or (hasattr(c,'solvent_content') and c.solvent_content)) and \
               (has_monomers_asym is None or (hasattr(c,'monomers_asym') and c.monomers_asym)) and \
               (has_residues_mon is None or (hasattr(c,'residues_mon') and c.residues_mon))    and \
               (has_num_atoms is None or (hasattr(c,'exp_num_atoms') and c.exp_num_atoms is not None)) and \
               (has_wavel is None or (hasattr(c,'wavel') and c.wavel))                         and \
               (filename is None or any(f.name==filename for f in c.files))                    and \
               (filetype is None or any(c.GetFile(t,try_convert,inp_cont,conv_opts) for t in filetype)) and \
               (col is None or all(c.GetLabel(l,True) for l in col))                           and \
               (label is None or all(l in c.GetAllLabels(labels_only=True) for l in label))    and \
               (exclude_cont is None or c not in exclude_cont)                                     \
             ):
            filtered.append(c)
            if last:
              break
          if filetype and try_convert:
            inp_cont=inp_cont_orig[:]
      else:
        common.Error('Wrong data class name: {0} - object cannot be filtered.'.format(container_name))
    if stored_order:   # reverse back - perhaps not the most efficient
      filtered.reverse()
    return filtered

  def ClearAll(self,propagate=True):
    """Removes any input/output containers currently stored, leaving empty lists storage"""
    for data_class in data.data_container.__subclasses__():
      self.Clear(data_class.__name__,propagate)

  def Clear(self,container_name,propagate=True):
    """Cleares the specified container type list, removing any objects stored in the list"""
    try:
      setattr(self,container_name,[])
    except AttributeError:
      common.Error('Wrong data class name: {0} - object cannot be cleared.'.format(container_name))
    if propagate:
      for o in self.FindPropagateObjects():
        if self.is_output:
          o.out.Clear(container_name,propagate)
        else:
          o.inp.Clear(container_name,propagate)

  def Set(self, container_list, propagate=True, checkinpfile=True, copy=False):
    """Sets the inputted container or container list as (the only) input/output of that container (nick)name
       If container list is inputted then the containers in it need to be of the same (nick)name
       Any previously stored containers of the container nick(name) are not stored anymore
       It is a warning if empty container list is inputted (there is Clear() member function for that)
    """
    if not container_list:
      inout_str='output'  if self.is_output  else 'input'
      common.Warning('No container to be set for {0} of {1}.'.format(inout_str,self.parent.nick))
      return
    try:
      container_list[0]
    except (TypeError,IndexError):
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      container_list=[container_list,]
    container_list=list(filter(None,container_list))
    try:
      cont_name = container_list[0].nick
    except AttributeError:
      common.Error('Attempting to set wrong or non-existing data class {0} to {1}.'.format(container_list[0], self.parent.nick))
    self.Clear(cont_name,propagate)
    for cont in container_list:
      if cont_name != cont.nick:
        common.Error('Setting container problem - inputted objects are of different type: {0}, {1}'.format(
               cont_name,cont.nick))
      else:
        if copy:
          self.AddCopy(cont,propagate,checkinpfile)
        else:
          self.Add(cont,propagate,checkinpfile)

  def SetCopy(self, container_list, propagate=True, checkinpfile=True):
    """The same as Set() but the inputted container(s) are copied (new instances created)"""
    self.Set(container_list, propagate, checkinpfile, copy=True)

  def Delete(self, container, propagate=True):
    """Deletes the inputted container from this input/output
       Leaves gracefully if no such container found in this input/output
    """
    for data_class in data.data_container.__subclasses__():
      if container in getattr(self, data_class.__name__):
        getattr(self, data_class.__name__).remove(container)
        if propagate:
          for o in self.FindPropagateObjects():
            if self.is_output:
              o.out.Delete(container,propagate)
            else:
              o.inp.Delete(container,propagate)

  def GetCrystalsList(self):
    """Returns a list of all the crystal names of (all the containers of) this input/output"""
    return list(set([ o.xname  for o in self.GetAll() ]))

  def GetDataList(self, crystal=None):
    """Returns a list of the dataset names of this input/output
       If 'crystal' is not None then only dataset names from that crystal are returned
    """
    return list(set([ o.dname  for o in self.GetAll(xname=crystal) ]))

  def GetTargetObjects(self, target, N=1, no_model=False, modtyp=None, fsftyp=None, inten=False, \
                       no_warn=False, accept_anom_nat=False, accept_anom_nosubstr=False):
    """Returns objects and N of this input/output suitable for inputted crystallographic llhood target,
       target=SAD   : f+, f-, mod (/substr)
       target=SIRAS : fn, f+, f-, modn (native-optional), mod (/substr, der.)
       target=MAD   : f1+, f1-, ... fN+, fN-, mod (/substr)
       target=MLHL  : f, hl, mod (optional)
       target=RICE  : f, mod
       (None,0) is returned if was not possible to find the objects above
       
       arguments:
       target   : one of {"SAD","SIRAS","MAD","MAD-2W","MAD-3W","MAD-4W","MLHL","RICE"}
       N        : number of anomalous datasets (default 2 for MAD unless '3W','4W' specified)
       no_model : if True then model is not required (default False)
       modtyp   : model needs to be of the specified filetype (default None)
       fsftyp   : reflection files need to be of the specified filetype (default None)
       inten    : if True then I's are searched instead of F's (default False); if None then both are accepted
       no_warn  : if True then no warning is printed about objects that could not be found (default False)
       accept_anom_nat : if True then data flagged native are accepted as anomalous for SAD/MAD/SIRAS
       accept_anom_nosubstr : if True then model without substr (partial) is accepted for SAD/MAD
    """
    if not target:
      return None,0
    N_strict = False
    if target.startswith('MAD') and N<2:
      N = 2
      for n,nw in enumerate(('2W','3W','4W')):
        if nw in target:
          N, N_strict = n, True
    if target.startswith('MAD'):
      target='MAD'
    fi='f'  if not inten  else 'i'
    if inten is None:
      fi=None
    missing_obj,fnat_xname,fnat_dname="",None,None
    anom,nonanom,hl={'f+':None,'f-':None},{'fn':None,'modn':None},{'hl':None}
    #if target=='SIRAS':
    fnat = self.Get('fsigf',typ='average',col=fi,filetype=fsftyp,is_native=True)
    if fnat:
      fnat_xname, fnat_dname = fnat.GetCrystalName(), fnat.GetDataName()
      nonanom = { 'fn':fnat, 'modn':self.Get('model',xname=fnat_xname,filetype=modtyp) }
    else:
      fnat_unusable = self.Get(is_native=True)
      if fnat_unusable:
        fnat_xname, fnat_dname = fnat_unusable.GetCrystalName(), fnat_unusable.GetDataName()
    if target in ('SAD','MAD','SIRAS'):
      for crys in self.GetCrystalsList():
        fplus = self.GetAll('fsigf',typ='plus',col=fi,xname=crys,filetype=fsftyp)
        fplus = [fp1  for i,fp1 in enumerate(fplus) if fp1.dname not in [fplus[j].dname for j in range(i)]]
        fminu = self.GetAll('fsigf',typ='minus',col=fi,xname=crys,filetype=fsftyp)
        fminu = [fp1  for i,fp1 in enumerate(fminu) if fp1.dname not in [fminu[j].dname for j in range(i)]]
        if not accept_anom_nat:
          fplus = [fp1  for i,fp1 in enumerate(fplus) if fp1.dname!=fnat_dname or fp1.xname!=fnat_xname]
          fminu = [fp1  for i,fp1 in enumerate(fminu) if fp1.dname!=fnat_dname or fp1.xname!=fnat_xname]
        if accept_anom_nosubstr:
          mdl = self.Get('model',xname=crys,filetype=modtyp)
        else:
          mdl = self.Get('model',typ=('substr','partial+substr'),xname=crys,filetype=modtyp)
        if len(fplus)>=N and len(fminu)>=N:
          if mdl or no_model:
            if not N_strict:
              N = min(4,len(fplus),len(fminu))
            anom = dict([ ('f'+str(n+1)+str(pm),f)  for n,p in enumerate(zip(fplus[:N],fminu[:N])) \
                                                    for pm,f in zip(('+','-'),p) ])
            anom['f+'],anom['f-']=anom['f1+'],anom['f1-']
            anom['mod']=mdl
          else:
            missing_obj='model'
          continue
        if not nonanom: #and target=='SIRAS':
          fnat = self.Get('fsigf',typ='average',col=fi,xname=crys,filetype=fsftyp)
          if fnat:
            nonanom = { 'fn':fnat, 'modn':self.Get('model',xname=crys,filetype=modtyp) }
      if anom['f+']:
        if target!='SIRAS' or (nonanom and nonanom['fn']):
          merge = dict( list(nonanom.items()) + list(anom.items()) )
          merge = dict( (k,v) for k,v in merge.items() if v is not None )
          return merge, N
        missing_obj='Fnat'
      elif missing_obj != 'model':
        missing_obj='F+ and/or F-'
    elif target=='MLHL':
      f = self.Get('fsigf',typ='average',col=fi,custom=('mlhl'),filetype=fsftyp)
      if not f:
        f = fnat  if fnat  else self.Get('fsigf',typ='average',col=fi,filetype=fsftyp)
      hl = self.Get('mapcoef',col='hla',custom='mlhl')
      phfom = self.Get('mapcoef',col=('ph','fom'),custom='mlhl')
      if not hl and not phfom:
        hl = self.Get('mapcoef',col='hla')
        phfom = self.Get('mapcoef',col=('ph','fom'))
      mdl=self.Get('model',filetype=modtyp)
      if (hl or phfom) and f and (mdl or no_model):
        return {'f':f,'hl':hl,'phfom':phfom,'mod':mdl},N
      if not mdl: missing_obj='model'
      if not hl and not phfom:  missing_obj='phase distribution'
      if not f:   missing_obj='F'
    elif target=='RICE':
      f = self.Get('fsigf',typ='average',col=fi,custom=('rice'),filetype=fsftyp)
      if not f:
        f = fnat  if fnat  else self.Get('fsigf',typ='average',col=fi,filetype=fsftyp)
      mdl=self.Get('model',typ=('partial','partial+substr'),filetype=modtyp)
      if f and (mdl or no_model):
        return {'f':f,'mod':mdl},N
      if not mdl: missing_obj='model'
      if not f:   missing_obj='F'
    else:
      common.Error("Target {0} not known.".format(target))
    if missing_obj and not no_warn:
      if missing_obj=='model' and modtyp:
        missing_obj='{0} model'.format(modtyp)
      if missing_obj!='model' or not no_model:
        common.Warning('{0} could not be found for target {1}'.format(missing_obj,target))
    return None,0
