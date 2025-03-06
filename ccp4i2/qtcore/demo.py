from __future__ import print_function

'''
This is an incomplete example of a CInt class. This class should have the functionality of the standard Python int
and of the ccp4 CData class. The key functionality of CData is data validation and read/write to formats such as XML.
The original proposal was to subclass CInt from CData and the Python int class. But this will not work because int
is immutable - it is given a value on instantiation and then can not be changed. This is undesirable behavior for
a parameter in a gui!

So the best solution is to have the actual value of the integer as an attribute of class (called _value in example
below) and to implement all of the usual methods of int ( __add__, __cmp__ etc.) for CInt by applying the to _value.
This is done for some of the int methods below.

A required feature is that CInt can, optionally, be undefined.  I am using None as the undefined value. Many standard
methods of int will not work for an undefined value. How should we handle this?  I think the only safe approach is to
attempt the required action and let it throw an exception.
There is an isSet() method that programmers can use before using the int-like methods.

This example uses qualifiers to set parameters such as default values and validation criteria. This allows
some customisation of the class without requiring writing Python code to subclass. A data type can be defined
as the data class and the qualifiers. These data type definitions could be in an 'external' file.
'''
import types,copy
from lxml import etree

SEVERITY_OK = 0
SEVERITY_UNDEFINED = 1
SEVERITY_UNDEFINED_INVALID = 2
SEVERITY_WARNING = 3
SEVERITY_INVALID = 4

class CValidityReport():
  '''
  Holds list of errors and warnings returned by validity() methods
  '''
  def __init__(self):
    self._reports = []

  def append(self,cls=None,code=0,name='' ):
    self._reports.append ( { 'class' : cls, 'name' : name, 'code' : code } )
    
  def extend(self,other=None):
    self._reports.extend(other._reports)
    
  def prependName(self,name=''):
    #print 'prependName',name,self._reports
    for item in self._reports:
      if len(item['name'])>0:
        item['name'] = name + '_' + item['name']
      else:
        item['name'] = name
      
  def maxSeverity(self):
    maxSev = 0
    for report in self._reports:
      code = report['code']
      if code<100:
        severity = CData.VALIDITY_CODES[report['code']]['severity']
      else:
        severity = report['class'].VALIDITY_CODES[report['code']]['severity']
      maxSev = max(maxSev,severity)
    return maxSev

  def __len__(self):
    return len(self._reports)

  def __str__(self):
    text = ''
    for report in self._reports:
      if report['code']<100:
        severity = CData.VALIDITY_CODES[report['code']]['severity']
        desc = CData.VALIDITY_CODES[report['code']]['description']
      else:
        severity = report['class'].VALIDITY_CODES[report['code']]['severity']
        desc = report['class'].VALIDITY_CODES[report['code']]['description']
    
      text = text + "{0:20} {1:20} {2:3} {3:3} -{4}-\n".format(report['name'],report['class'],report['code'],severity,desc)
    return text[0:-1]

# A base class for all data classes 

class CData():
    CONTENTS = {}
    QUALIFIERS = { }
    VALIDITY_CODES = { 0 : { 'severity' : SEVERITY_OK,
                             'description' : 'OK' },
                       1 : { 'severity' : SEVERITY_UNDEFINED,
                             'description' : 'warning undefined' },
                       2 : { 'severity' : SEVERITY_UNDEFINED_INVALID,
                             'description' : 'invalid undefined' },
                       3 : { 'severity' : SEVERITY_WARNING,
                             'description' : 'warning missing data' },
                       4 : { 'severity' : SEVERITY_INVALID,
                             'description' : 'missing data' },
                       5 : { 'severity' : SEVERITY_INVALID,
                             'description' : 'wrong type' } }
    
    def __init__(self,value={},qualifiers={},**kw):
       qualifiers.update(kw)
       self.build()
       self._qualifiers = {}
       self.setQualifiers(qualifiers)
       self.setDefault()
       # Using the variable valu to overcome weird bug that a second
       # instantiation of class had values that had been passed to first
       # instance via the **kw mechanism
       if isinstance(value,dict):
         valu = {}
         valu.update(value)
         valu.update(kw)
         self.set(value=valu)
       else:
         self.set(value=value)
         
    def build(self):
        self._value = {}
        for key,defn in list(self.CONTENTS.items()):
            self._value[key] = defn['class'](qualifiers=defn.get('qualifiers',{}))
       
    def setQualifiers(self,qualifiers={},**kw):
        qualifiers.update(kw)
        #print 'CData.setQualifiers',self,qualifiers
        if 'dataType' in qualifiers:
            self._dataType = qualifiers['dataType']
            del qualifiers['dataType']
        for key,value in list(qualifiers.items()):
          objname,newkey=self.splitName(key)
          if newkey is not None:
            obj = self._value.get(objname,None)
            if obj is None:
              print('error interpreting qualifier:',key,'unknown object:',qlist[0])
              break
            else:
              obj.setQualifiers( { newkey : value } )
          else:
            if key in self.QUALIFIERS:
              self._qualifiers[key] = value
            
    def qualifiers(self,name=None):
      #return self._qualifiers.get(name,self.QUALIFIERS[name])
      #print 'CData.qualifier',self,name

      # If input name is set then return item of that name
      if name is not None:
        if name in self.QUALIFIERS:
          return self._qualifiers.get(name,self.QUALIFIERS[name])
        else:
          objname,newname = self.splitName(name)
          if newname is not None and objname in self._value:
              return self._value[objname].qualifiers(newname)
          else:
              return None

      # return 'flattened' dict of all qualifiers
      ret = {}
      # The following two lines all qualifiers (incuding class defaults)
      #for key in self.QUALIFIERS.keys():
      #  ret[key] = self._qualifiers.get(key,self.QUALIFIERS[key])
      ret.update(self._qualifiers)
      for key in list(self.CONTENTS.keys()):
        #Everything in self._value should be a CData..
        obj = self._value[key]
        if isinstance(obj,CData):
          # The following two lines give a flattened dict
          #qual = obj.qualifiers()
          #for k,v in qual.items(): ret[key+'_'+k]=v
          # Nest dict
          ret[key] = obj.qualifiers()
      return ret

    def setDefault(self):
      if 'default' in self._qualifiers:
        default =  self._qualifiers['default']
      elif 'default' in self.QUALIFIERS:
        default =  self.QUALIFIERS['default']
      else:
        if len(self.CONTENTS)>0:
          default = {}
        else:
          default = None
      self.set(default)

      
    def validity(self,arg={}):
      if isinstance(arg,CData):
        testData = arg.get()
      else:
        testData = {}
        testData.update(arg)
      #print 'validity',self,testData
      validityObj = CValidityReport()
      for key in list(self.CONTENTS.keys()):
        if key in testData:
         testValidity = self._value[key].validity(testData.get(key))
         testValidity.prependName(key)
         validityObj.extend(testValidity)
      return validityObj

    def fix(self,arg={},**kw):
      # Apply possible fixes to input data
      # This should be re-implemented in derived classes
      arg.update(kw)
      return arg

    def __getattr__(self,name):
      #print 'in CData.__getattr__',name
      if name in self.CONTENTS:
        return self._value[name]
      else:
        raise AttributeError(name)

    def splitName(self,name):
      import re
      s = re.match(r'(.*?)_(.*)',name)
      if s is None:
        return [name,None]
      else:
        return s.groups()

    def isInstance(self,other):
      # Is other an object identical to self?
      # This is inadequate test since could have different qualifiers
      return isinstance(other,self.__class__)
        
      
    def set(self,value={},**kw):
      print('CData.set',self,value)

      # Save the current value  to restore if validation fails
      saveValue = copy.deepcopy(self._value)
      
      if self.isInstance(value):
        for key in self.CONTENTS:
          other = value.get(key)
          self._value[key].set(other)
        return CValidityReport()
      else:    
        value.update(kw)
        self.copyValue(value)

        validityObj = self.validity(self._value)
        # Validity check failed - restore original
        if validityObj.maxSeverity()>1:
          self.copyValue(saveValue)
        else:
          self.dataChanged.emit()
        return validityObj

    def copyValue(self,value):
      for key,val in list(value.items()):
        objname,newkey = self.splitName(key)
        if newkey is not None:
          obj = self._value.get(objname,None)
          if obj is None:
            print('error interpreting value:',key,'unknown object:',objname)
          else:
            #print 'CData.set newkey',obj,newkey,val
            obj.set( value = { newkey : val } )
        else:
          if key in self.CONTENTS:
            self._value[key].set(val)

    def get(self,name=None):
      #print 'CData.get',self,name
      # If this is not a composite data class then just return the single
      # data item
      if len(self.CONTENTS)<=0: return self._value

      # If input name is set then return item of that name
      if name is not None:
        if name in self.CONTENTS:
          return self._value.get(name)
        else:
          objname,newname = self.splitName(name)
          if newname is not None and objname in self._value:
              return self._value[objname].get(newname)
          else:
              return None

      # return 'flattened' dict of all data
      ret = {}
      for key,value in list(self._value.items()):
        #Everything in self._value should be a CData..
        if isinstance(value,CData):
          val = value.get()
          if isinstance(val,dict):         
            for k,v in list(val.items()): ret[key+'_'+k]=v
          else:
            ret[key] = val
        else:
          print('Error - item in ',self,'is not an instance of CData:',key)
      return ret
    
    def get0(self,name=None):
      #print 'CData.get',self,name
      # If this is not a composite data class then just return the single
      # data item
      if len(self.CONTENTS)<=0: return self._value

      # If input name is set then return item of that name
      if name is not None:
        if name in self.CONTENTS:
          return self._value.get(name)
        else:
          if name in self._value:
            return self._value[name]
          
      # return nested dict of all data
      ret = {}
      for key,value in list(self._value.items()):
        #Everything in self._value should be a CData..
        if isinstance(value,CData):
          val = value.get()
          ret[key] = val
        else:
          print('Error - item in ',self,'is not an instance of CData:',key)
      return ret

    def dataType(self):
      cls = str(self.__class__)
      return cls.split('.')[-1]
      

    def getEtree(self):
      element = etree.Element(self.dataType())
      if len(self.CONTENTS) == 0:
        element.text = str(self._value)
      else:
        for key in list(self.CONTENTS.keys()):
          ele = self._value[key].getEtree()
          ele.set('id',key)
          element.append(ele)
      return element
      
    def xmlText(self):
      element = self.eTree()
      text = etree.tostring(element,pretty_print=True, xml_declaration=True)
      return text

    def setEtree(self,element):
      if len(self.CONTENTS) == 0:      
        rv = self.set(str(element.text))
      else:
        for ele in element.iterchildren():
          name = ele.tag
          ele_id = str(ele.get('id'))
          if ele_id in self.CONTENTS:
            self._value[ele_id].setEtree(ele)

    def getQualifiersEtree(self):
      element = etree.Element(self.dataType())
      #element = etree.Element('qualifiers')
      for key,value in list(self._qualifiers.items()):
        ele = etree.Element(key)
        ele.text = str(value)
        element.append(ele)
      return element
      
    def setQualifiersEtree(self,element):
      for ele in element.iterchildren():
        name = ele.tag
        if name in self.QUALIFIERS_DESCRIPTION:
          if self.QUALIFIERS_DESCRIPTION[name]['type'] in [int,float,bool,bool,str]:
            self._qualifiers[ name] = self.QUALIFIERS_DESCRIPTION[name]['type'](ele.text)
          elif self.QUALIFIERS_DESCRIPTION[name]['type'] is list:
            strList = str(ele.text).split(',')
            qList = []
            for item in strList:
              qList.append(self.QUALIFIERS_DESCRIPTION[name]['subType'](item))
            self._qualifiers[name] = qList
                                              
        elif name in self.CONTENTS:
          self._value[name].setQualifiersEtree(ele)
         

    def emitDataChanged(self):
      # Emit Qt signal
      pass
      
class CInt(CData):

    # Qualifiers determine default values and limits on allowed values.
    # They allow some customisation without requiring writing a new Python class.
    # A data type could be defined in an external file by specifying the appropriate class and qualifiers.
    #
    # enumerators are either
    #    - a list of allowed values if strict_enumerators is True
    #    - a list of recomended values if strict_enumerators is False (useful in gui)
    #
    CONTENTS = { } 
    QUALIFIERS = { 'allowUndefined' : True,
                   'max' : None,
                   'min' : None,
                   'default' : 0,
                   'enumerators' : [],
                   'strictEnumerators' : False }
    
    QUALIFIERS_DESCRIPTION = { 'allowUndefined' : { 'type' : bool},
                               'max'  : { 'type' :int},
                               'min'  : { 'type' :int},
                               'default'  : { 'type' :int},
                               'enumerators'  : { 'type' :list},
                               'strictEnumerators' : { 'type' :bool} }
    
    VALIDITY_CODES = { 101 : { 'severity' : SEVERITY_INVALID,
                              'description' : 'below minimum' },
                       102 : { 'severity' : SEVERITY_INVALID,
                              'description' : 'above maximum' },
                       103 : { 'severity' : SEVERITY_INVALID,
                              'description' : 'not one of limited allowed values' } }
                        
    def __init__(self,value=None,qualifiers={},**kw):
       qualifiers.update(kw)
       CData.__init__(self,value=value,qualifiers=qualifiers)
       
    def validity(self,arg):
        # return 0=OK 1=undefined 2+=error
        validityObj = CValidityReport()
        arg = self.fix(arg)
        if arg is None:
          if self.qualifiers('allowUndefined'):
            validityObj.append(self.__class__,1)
          else:
             validityObj.append(self.__class__,2)
          return validityObj
        if not isinstance(arg,int):
            validityObj.append(self.__class__,5)
            return validityObj
        if self.qualifiers('min') is not None and arg < self.qualifiers('min') : validityObj.append(self.__class__,101)
        if self.qualifiers('max') is not None and arg > self.qualifiers('max') : validityObj.append(self.__class__,102)
        if self.qualifiers('strictEnumerators') and self.qualifiers('enumerators').count(arg)<1: validityObj.append(self.__class__,103)
        return validityObj

    def build(self):
       self._value = None
       
    def set(self,value=None):
       #print 'CInt.set',value
       value=self.fix(value)
       v = self.validity(value)
       if v.maxSeverity()<=1:
           self._value = value
       return v

    def fix(self,arg):    
       if isinstance(arg,CInt):
         arg = arg.get()
       elif isinstance(arg,(str,float)):
         try:
           arg = int(arg)
         except:
           raise AttributeError('CInt')
       return arg

    def get(self):
      return self._value
        
    # Applying mathematical functions to self._value in place
    # These methods return a CValidityReport
    def add(self,arg):
      if self._value is not None:
        return self.set(self._value + int(arg))
    def sub(self,arg):
      if self._value is not None:
        return self.set(self._value - int(arg))
    def mul(self,arg):
      if self._value is not None:
        return self.set(self._value * int(arg))
    def div(self,arg):
      if self._value is not None:
        return self.set(self._value / int(arg))
    def abs(self):
      if self._value is not None:
        return self.set(abs(self._value))

    # Implementation of the usual methods for int
    def __abs__(self):
      return self.__class__(value=self._value.__abs__(int(arg)),qualifiers=self._qualifiers)
    def __add__(self,arg):
      return self.__class__(value=self._value.__add__(int(arg)),qualifiers=self._qualifiers)
    def __and__(self,arg):
      return self._value.__and__(int(arg))   
    def __cmp__(self,arg):
      # Will throw exception if self.value not set
      return self._value.__cmp__(int(arg))
    def __coerce__(self,arg):
      print('NOT IMPEMENTED: CInt.__coerce__')
    def __div__(self,arg):
      return self.__class__(value=self._value.__div__(int(arg)),qualifiers=self._qualifiers)
    def __divmod__(self,arg):
      return self._value.__divmod__(int(arg))
    def __float__(self):
      return self._value.__float__()
    def __floordiv__(self,arg):
      return self.__class__(value=self._value.__floordiv__(int(arg)),qualifiers=self._qualifiers)
    def  __format__(self,arg):
      return self._value.__format__(arg)
    #def __getattribute__(self,arg):
    #  print 'NOT IMPEMENTED: CInt.__getattribute__'
    def __getnewargs__(self,arg):
      return self._value.__getnewargs__()
    #def __hash__(self):
    #  print 'NOT IMPEMENTED: CInt.__hash__'
    def __hex__(self):
      return self._value.__hex__()
    def __index__(self):
      return self._value.__index__()
    def __int__(self):
      return self._value.__int__()
    def __invert__(self):
      return self.__class__(value=self._value.__invert__(),qualifiers=self._qualifiers)
    def __long__(self):
      return self._value.__long__()
    def __lshift__(self):
      print('NOT IMPEMENTED: CInt.__lshift__')
    def __mod__(self,arg):
      return self.__class__(value=self._value.__mod__(int(arg)),qualifiers=self._qualifiers)
    def __mul__(self):
      return self.__class__(value=self._value.__mul__(int(arg)),qualifiers=self._qualifiers)
    def __neg__(self):
      print('NOT IMPEMENTED: CInt.__neg__')
    def __bool__(self):
      print('NOT IMPEMENTED: CInt.__nonzero__')
    def __oct__(self):
      return self._value.__oct__()
    def __or__(self,arg):
      return self._value.__or__(int(arg))   
    def __pos__(self):
      print('NOT IMPEMENTED: CInt.__pos__')
    def __pow__(self):
      return self.__class__(value=self._value.__pow__(int(arg)),qualifiers=self._qualifiers)   
    def __radd__(self,arg):
      return self.__class__(value=self._value.__radd__(int(arg)),qualifiers=self._qualifiers)
    def __rand__(self,arg):
      return self._value.__rand__(int(arg))
    def __rdiv__(self,arg):
      return self.__class__(value=self._value.__rdiv__(int(arg)),qualifiers=self._qualifiers)
    def __rdivmod__(self,arg):
      return self._value.__rdivmod__(int(arg))
    #def __reduce__(self):
    #  print 'NOT IMPEMENTED: CInt.__reduce__'
    #def __reduce_ex__(self):
    #  print 'NOT IMPEMENTED: CInt.__reduce_ex__'
    def __repr__(self):
      return self._value.__repr__()
    def __rfloordiv__(self,arg):
      return self.__class__(value=self._value.__rfloordiv__(int(arg)),qualifiers=self._qualifiers)
    def __rlshift__(self):
      print('NOT IMPEMENTED: CInt.__rlshift__')
    def __rmod__(self):
      return self.__class__(value=self._value.__rmod__(int(arg)),qualifiers=self._qualifiers)
    def __rmul__(self):
      return self.__class__(value=self._value.__rmul__(int(arg)),qualifiers=self._qualifiers)
    def __ror__(self,arg):
      return self._value.__ror__(int(arg))   
    def __pow__(self):
      return self.__class__(value=self._value.__rpow__(int(arg)),qualifiers=self._qualifiers)   
    def __rrshift__(self):
      print('NOT IMPEMENTED: CInt.__rrshift__')
    def __rshift__(self):
      print('NOT IMPEMENTED: CInt.__rshift__')
    def __rsub__(self):
      return self.__class__(value=self._value.__rsub__(int(arg)),qualifiers=self._qualifiers)
    def __rtruediv__(self,arg):
      return self._value.__rtruediv__(int(arg))
    def __or__(self,arg):
      return self._value.__rxor__(int(arg))   
    def __str__(self):
      return self._value.__str__()
    def __sub__(self,arg):
      return self.__class__(value=self._value.__sub__(int(arg)),qualifiers=self._qualifiers)
    def __truediv__(self,arg):
      return self._value.__truediv__(int(arg))
    def __trunc__(self):
       print('NOT IMPEMENTED: CInt.__trunc__')
    def __xor__(self,arg):
      return self._value.__xor__(int(arg))   
    #'conjugate', 'denominator', 'imag', 'numerator', 'real'

class CIntRange(CData):

  CONTENTS = { 'start' :  { 'class' : CInt },
               'end' :  { 'class' : CInt } }
  QUALIFIERS = { 'compare' : None }
  QUALIFIERS_DESCRIPTION = { 'compare' : { 'type' : int } }
  VALIDITY_CODES = { 101 : { 'severity' : SEVERITY_INVALID,
                             'description' : 'End of range less than start' } }


  def validity(self,arg):

    # Demo of validity checking when class contents are not orthogonal
    # if qualifier compare is set then relation between start and end is limited
    v = CData.validity(self,arg)
    if v.maxSeverity()>0: return v
    
    if 'compare' in self._qualifiers:
      if self._value['end'].__cmp__(self._value['start']) != self._qualifiers['compare']:
        v.append(self.__class__,101)
    return v
      
  def fix(self,arg={}):
    if 'compare' in self._qualifiers:
      try:
        if arg['end'].__cmp__(arg['start']) != self._qualifiers['compare']:
          ret = {}
          ret['start'] = arg['end']
          ret['end'] = arg['start']
          return ret
      except:
        pass
    return arg
     

class CSummat(CData):

    CONTENTS =  { 'nCycles' : { 'class' : CInt ,
                                'qualifiers' : { 'allowUndefined' : False,
                                                 'default' : 3,
                                                 'max' : 20,
                                                 'min' : 1 } },
                  'cutoff' : { 'class' : CInt,
                               'qualifiers' : { 'max' : 100,
                                                'min' : 1 },
                               'default' : 5 },
                  'range' : { 'class' : CIntRange,
                              'qualifiers' :  { 'start_min' : 0,
                                                'end_max' : 100,
                                                'compare' : 1 }
                            }
                  }
    QUALIFIERS = {}

      
    pass

class foo:
    def __init__(self,valu={},**kw):
        value = {}
        value.update(valu)
        value.update(kw)
        print('foo.value',value)
     
def demo():

    t = CInt(5,max=10,min=0)
    l = ['zero','one','two','three','four','five','six']
    print('t = CInt(5,max=10,min=0)   t+6:',t + 6,'t>8:',t>8)
    print("['zero','one','two','three','four','five','six'][t]:",['zero','one','two','three','four','five','six'][t])
    t.set(12)
    print('t.set(12) t is:', t)
    x = CSummat()
    print('CSummat: ',x.nCycles,x.cutoff)
    
