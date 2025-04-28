import ast
import copy
import ctypes
import io
import os
import shutil
import sys
import types
import xml.etree.ElementTree as ET


class CrankError(Exception):
  def __init__(self,message=None):
    self.message=message
    Warning(message)
    Exception.__init__(self, message)

class Unsuccessful(Exception):
  def __init__(self,message=None):
    self.message=message
    Warning(message)
    Exception.__init__(self, message)

# errors/warnings calling
def Error(message,filehandle=None,debug=True,nosuccess=False):
  if filehandle is None:
    filehandle=sys.stderr
  if debug: 
    if sys.exc_info()[0]:
      filehandle.write('Error: '+message+'\n')
      raise
    elif nosuccess:
      raise Unsuccessful(message)
    else:
      raise CrankError(message)
  elif debug is not None:
    filehandle.write('Error: '+message+'\n')
    sys.exit(1)
def Warning(message,filehandle=None):
  if filehandle is None:
    filehandle=sys.stdout
  filehandle.write('Warning: '+message+'\n')
def Info(message,filehandle=None,eol=True,also_stdout=False):
  handles=[]
  if filehandle is not None:
    handles.append(filehandle)
  if (also_stdout and sys.stdout not in handles) or not handles:
    handles.append(sys.stdout)
  if eol:
    message+='\n'
  for h in handles:
    try:
      h.write(message)
    except TypeError:
      # presumably python3 only requiring bytes
      h.write(bytes(message,'utf-8'))


def SymLink(target, copy, ignore_error=True, copy_fallback=False):
  """Create symbolic link, multiplatform"""
  if copy_fallback: ignore_error=True
  try:
    if os.name=='nt':
      ctypes.windll.kernel32.CreateSymbolicLinkA( copy, target, 0 )
    elif not os.path.lexists(copy):
      os.symlink(target, copy)
  except OSError as e:
    if not ignore_error: 
      raise
    Warning('Symlink could not be created due to error: {0}'.format(e))
    if hasattr(sys,'exc_clear'): sys.exc_clear()
    if copy_fallback:
      shutil.copyfile(target, copy)
    return 1
  return 0


# Automatically convert inputted string to other types 
# (the original string returned if no conversion occured)
def AutoConvert(strng):
  try:
    return ast.literal_eval(strng.strip())
  except (ValueError,SyntaxError,AttributeError):
    if hasattr(sys,'exc_clear'): sys.exc_clear()
    return strng


# define is_string()
try:
    basestring
    def is_string(s):
        return isinstance(s, basestring)
except NameError:
    def is_string(s):
        return isinstance(s, str)

# define is_file()
try:
    file
    def is_file(f):
        return isinstance(f, file)
except NameError:
    def is_file(f):
        return isinstance(f, io.IOBase)


def ReadXML(xmlfile):
  """Reads inputted XML file using ElementTree

     Argument:
     xmlfile -- the xml file to read from

     Returns root of the elementtree
  """
  if not xmlfile:
    Error("XML file {0} not defined".format(xmlfile))
  try:
    xtree = ET.ElementTree(file=xmlfile)
    xelem=xtree.getroot()
  except:
    Error("Could not parse XML file {0}, it may be broken.".format(xmlfile))
  return xelem

def WriteXML(obj, xmlfile=None):
  """Write XML file
     Argument:
     obj - the object to be written to the XML
     xmlfile -- the file to be writen, if None then the name of the object with xml suffix is used
  """
  xroot = obj.Data2XML()
  ET_Indent(xroot)
  xtree = ET.ElementTree(xroot)
  if xmlfile is None:
    xmlfile = obj.nick+'.xml'
  try:
    xtree.write(xmlfile,xml_declaration=True,encoding='utf-8')
  except TypeError:
    xtree.write(xmlfile,encoding='utf-8')
  except:
    if os.path.isfile(xmlfile):
      os.remove(xmlfile)
    Error("Could not write XML file of object {0}".format(obj.nick))

def ET_Indent(elem, level=0):
  """In-place prettyprint formatter for ElementTree"""
  i = "\n" + level*"  "
  if len(elem):
    if not elem.text or not elem.text.strip():
      elem.text = i + "  "
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
    for elem in elem:
      ET_Indent(elem, level+1)
    if not elem.tail or not elem.tail.strip():
      elem.tail = i
  else:
    if level and (not elem.tail or not elem.tail.strip()):
      elem.tail = i



class stats(object):
  # if multiple then a list of all matching stats will be returned, otherwise the last stat only
  # if convert then automatic conversion of the stat(s) will be attempted
  def __init__(self,name=None,regexp=None,xpath=None,multiple=False,attrib=None,convert=True,accept_none=False):
    self.name=name
    self.regexp=regexp
    self.xpath=xpath
    self.multiple=multiple
    self.attrib=attrib
    self.convert=convert
    self.accept_none=accept_none
  def ConvertIfAsked(self,value,convert=None):
    if convert is None:
      convert=self.convert
    if convert:
      return AutoConvert(value)
    return value


class parameter(object):
  """general program/process parameter base class

     For program, there are keyword and argument parameters - passed directly as program keyword 
     (ie input stream to the program) or argument (in the command line)
     There are also 'virtual' parameters for process or (rarely) program  - 
     these are processed by the code (often leading to setting of keywords/arguments)

     Value can be a string, scalar, vector, bool or None.
     Internally, values are always stored in a list allowing multiple setting of a parameter.
     For virtual parameters, there is no special automatic treatment 
     (treated by the specific program/process code).
     For keyword/arguments, the following holds for value types:
       None:   keyword was not set
       bool:   if False then keyword is not used, if True it is used without anything following it
       string: keyword is followed by the string 
       scalar: keyword is followed by the scalar converted to string 
       list/tuple/vector: keyword is followed by all the elements (as strings, separated by space)
  """
  is_key=None
  is_arg=None
  accessed=False
  # True if the parameter was specified directly by the input (xml or keyworded), ie by the user
  # if the parameter was set by default (or not set) then it is False
  inputted=False

  def __init__(self, desc='', typ=None, disable_none_type=True, cap=None, share=False):
    self.value=[]
    self.description=desc
    # type(s) restriction of the parameter
    self.typ=[]
    # ignore program modifier for this parameter - only used for program arguments
    self.ignore_modif=[]
    self.capitalized=cap
    if typ is not None:
      assert not is_string(typ)
      try:
        self.typ.extend(typ)
      except (AssertionError,TypeError):
        if hasattr(sys,'exc_clear'): sys.exc_clear()
        self.typ.append(typ)
      # this has been disabled (by default) since NoneType is not pickleable in python 2.x and thus not compatible with multiprocessing
      if not disable_none_type:
        self.typ.append(types.NoneType)
    # share the value of parameter par with children of the process having a parameter with the 
    # same name supported if the par has a shared_with_children attribute (recursive)
    self.shared_with_children = share

  def Set(self, value, is_key=False, is_arg=False, append=True, ind=None, inputted=False, ignore_modif=False):
    """Set the parameter to 'value' in the parameters dictionary
       Internally, always storing the value in a list
    """
    self.inputted=inputted
    self.multuse=0
    if append:
      self.value.append(value)
      if is_arg: self.ignore_modif.append(ignore_modif)
    else:
      if ind is None:
        self.value=[value,]
        if is_arg: self.ignore_modif=[ignore_modif,]
      else:
        try:
          self.value[ind]=value
        except IndexError:
          Error('Wrong index {0} when setting parameter value "{1}"'.format(ind,value))
        else:
          if is_arg: self.ignore_modif[ind]=ignore_modif
    if is_key is not None:
      self.is_key=is_key
    if is_arg is not None:
      self.is_arg=is_arg
    if (is_key or is_arg) and value is None:
      Error('Program keywords or arguments cannot be set to None.')

  def Unset(self, is_key=False, is_arg=False):
    """Unset the parameter - clear the parameters dictionary.
    """
    # could this cause unexpected issues?  use with caution and check
    self.value=[]
    if is_arg: self.ignore_modif=[]
    if is_key is not None:
      self.is_key=is_key
    if is_arg is not None:
      self.is_arg=is_arg

  def GetAll(self, capitalized=None):
    """Returns the list of all settings (assignments) of this parameter"""
    result=[]
    for i,v in enumerate(self.value):
      result.append( self.Get(ind=i,capitalized=capitalized) )
    return result

  def Get(self, ind=-1, as_list=False, capitalized=None):
    """Get the parameter value"""
    if capitalized is None and self.capitalized:
      capitalized = True
    try:
      result=self.value[ind]
    except IndexError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      return None
    if as_list:
      result=self.ConvertToList(result)
    if capitalized:
      result=self.GetCapitalized(result)
    return result

  def GetModifIgnorance(self, ind=-1):
    """Returns whether the argument modifier was set to be ignored"""
    try:
      result=self.ignore_modif[ind]
    except IndexError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      return None
    return result


  def ConvertToList(self, val):
    # converts the value (specified by ind in case the parameter is set multiple times) to a list
    lst=val
    try:
      if not is_string(lst):
        lst=list(lst)
    except TypeError:
      if hasattr(sys,'exc_clear'): sys.exc_clear()
    if not isinstance(lst,list):
      lst=[lst,]
    return lst

  def GetCapitalized(self, val):
    # Capitalize value of the inputted parameter (if existing and if possible)
    try:
      upval=val.upper()
    except (AttributeError,TypeError):
      if hasattr(sys,'exc_clear'): sys.exc_clear()
      upval=copy.copy(val)
      try:
        for i,v in enumerate(val):
          try:
            upval[i]=v.upper()
          except (TypeError,AttributeError):
            if hasattr(sys,'exc_clear'): sys.exc_clear()
      except:
        if hasattr(sys,'exc_clear'): sys.exc_clear()
    return upval
