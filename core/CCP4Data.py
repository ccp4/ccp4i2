from __future__ import print_function

"""
     CCP4Data.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the
     license to address the requirements of UK law.

     You should have received a copy of the modified GNU Lesser General
     Public License along with this library.  If not, copies may be
     downloaded from http://www.ccp4.ac.uk/ccp4license.php

     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
   Liz Potterton Aug 2010 - 'Generic' CCP4Data classes
"""

import sys
import re
import types
import time
import traceback

from PySide2 import QtCore

from core.CCP4ErrorHandling import *
from core.CCP4Config import QT, XMLPARSER
from core.CCP4QtObject import CObject
from lxml import etree
from xml.etree import ElementTree as ET

def isQualifier(cls, name=None):
    while issubclass(cls, CData):
        if name in cls.QUALIFIERS:
            return True
        cls = cls.__bases__[0]
    return False

def isCollectionClass(cls):
    if (CCollection,CDict,CList).count(cls):
        return True
    else:
        return issubclass(cls,(CCollection,CDict,CList))

def sortArguments(cls, args, aliases=[]):
    # Sort command line arguments into qualifiers and values and unknowns
    values = {}
    qualifiers = {}
    unknowns = {}
    for key,value in list(args.items()):
        ## BEWARE - will not work for run-time defined class contents (eg CProgramColumnGroup)
        if key in cls.CONTENTS:
            values[key] = value
        elif key in aliases:
            values[key] = value
        elif isQualifier(cls,key):
            qualifiers[key] = value
        elif key == 'build':
            pass
        else:
            objname,newkey=splitName(key)
            if newkey is None:
                unknowns[key] = value
            else:
                ####THIS WILL NOT WORK FOR RUN_TIME DEFINED CONTENTS EG CProgramColumnGroup
                subClsDefn = cls.CONTENTS.get(objname, None)
                if subClsDefn is None:
                    unknowns[key] = value
                else:
                    subCls = subClsDefn.get('class', None)
                    vs,qs,us = sortArguments(subCls, {newkey : value})
                    if len(vs) > 0:
                        values[key] = value
                    if len(qs) > 0:
                        if objname not in qualifiers:
                            qualifiers[objname] = {}
                        qualifiers[objname][newkey] = value
                    if len(us) > 0:
                        unknowns[key] = value
    return [values,qualifiers,unknowns]

def splitName(name):   # KJS : This fn. needs to be fixed or removed from the code-base.
    return [name, None]

def errorCode(cls, code):  # KJS : Revise ... this setup creates circ. dependency between data & error logging.
    while issubclass(cls, CData):
        if hasattr(cls, 'ERROR_CODES') and code in cls.ERROR_CODES:
            return cls.ERROR_CODES[code]
        cls = cls.__bases__[0]   # KJS : Revise... fail on multi-inhr.
    return {}

def errorCodeSeverity(cls, code):    # KJS - Revise
    err = errorCode(cls, code)
    if len(err) == 0:
        return -1
    elif 'severity' in err:
        return err['severity']
    else:
        return SEVERITY_ERROR # KJS : These globals need a repair job.

def errorCodeDescription(cls, code): # KJS : Needs revision.
    err = errorCode(cls, code)
    if len(err) == 0:
        return -1
    elif 'description' in err:
        return err['description']
    else:
        return ''

def baseClassList(cls):
    clsList = [cls]
    baseCls = clsList[-1].__bases__[0]
    while issubclass(baseCls, CData):
        clsList.append(baseCls)
        baseCls = clsList[-1].__bases__[0]
    clsList.reverse()
    return clsList

def classQualifier(cls, name=None):
    if isinstance(cls, dict):
        return cls.get('qualifiers',{}).get(name,None)
    else:
        return cls.QUALIFIERS.get(name,None)

def classQualifiers(cls):
    if isinstance(cls,dict):
        return cls.get('qualifiers',{})
    else:
        return cls.QUALIFIERS

class CDataQualifiers:
    PYTHONTYPE = str
    QUALIFIERS = {'allowUndefined' : True, 'default' : NotImplemented, 'toolTip' : NotImplemented,
                  'guiLabel' : NotImplemented, 'guiDefinition' : {}, 'helpFile' : NotImplemented, 'saveToDb' : False }
    QUALIFIERS_ORDER = ['allowUndefined', 'default', 'toolTip', 'guiLabel', 'guiDefinition', 'helpFile', 'saveToDb']
    QUALIFIERS_DEFINITION = {'allowUndefined' :{ 'type' : bool}, 'default' : { 'type' : dict},
                             'toolTip' : { 'type' : str}, 'guiLabel' : { 'type' : str},
                             'guiDefinition' : { 'type' : dict}, 'helpFile' : {'type' : str},
                             'saveToDb' : {'type' : bool, 'description' :'Save this data in the database'}}

    def setQualifiers(self, qualifiers={}, validateDefault=True, **kw):
        if len(qualifiers) == 0 and len(kw) == 0:
            return
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        for key,value in list(qualis.items()):
            if isQualifier(self.__class__,key) and key != 'default':
                self.setQualifier(key, value)
            elif key in self.contents() and isinstance(value, dict):
                obj = self._value.get(key, None)
                if obj is not None:
                    obj.setQualifiers(qualifiers=value)
            else:
                objname, newkey=splitName(key)
                if newkey is not None:
                    obj = self._value.get(objname, None)
                    if obj is None:
                        print('error interpreting qualifier:', key,'unknown object:', objname)  #, qlist[0] # KJS : There is no qlist !
                    else:
                        obj.setQualifier(newkey, value)
        # New set default value after validating it
        if 'default' in qualis and qualis['default'] is not NotImplemented and qualis['default'] is not None:
            default = self.coerce(qualis['default'])
            if validateDefault:
                v = self.validity(default)
                if v.maxSeverity() > SEVERITY_WARNING:
                    e = CException(CData, 6, name=self.objectPath(), details='setting default: ' + str(default))
                    e.extend(v)
                    raise e
            self.__dict__['_qualifiers']['default'] = default
            if isinstance(qualifiers['default'], dict):
                for key,value in list(qualifiers['default'].items()):
                    obj = self._value.get(key, None)
                    if obj is not None:
                        obj.setQualifiers(qualifiers={'default':value})
        else:
            pass

    def qualifierListType(self):
        return self.PYTHONTYPE

    def setQualifier(self, name, value):
        defn = self.qualifiersDefinition(name)
        if defn is None:
            raise CException(self.__class__, 17, 'Object: ' + self.objectPath() + ' qualifier: ' + name)
        dataType = defn.get('type', self.PYTHONTYPE)
        if value is None or value is NotImplemented:
            self.__dict__['_qualifiers'][name] = value
            return
        elif not isinstance(value, dataType):
            raise CException(self.__class__, 18, 'Object: ' + self.objectPath() + ' qualifier: ' + name)
        if dataType == list:
            listItemType = defn.get('listItemType', self.qualifierListType())
            for item in value:
                if not isinstance(item, (listItemType, type(NotImplemented), type(None))):
                    raise CException(self.__class__, 19, 'Object: ' + self.objectPath() + ' qualifier: ' + name)
        if value == self.qualifiers(name, custom=False):
            if name in self.__dict__['_qualifiers']:
                del self.__dict__['_qualifiers'][name]
        else:
            self.__dict__['_qualifiers'][name] = value

    def qualifiers(self, name=None, default=True, custom=True, contentQualifiers=True):
        # Return class default or custom qualifier value dependent on flags
        # If input name is set then return item of that name
        if name is not None:
            if custom and name in self.__dict__['_qualifiers']:
                return self.__dict__['_qualifiers'][name]
            if default:
                cls = self.__class__
                while issubclass(cls, CData):
                    if name in cls.QUALIFIERS:
                        return cls.QUALIFIERS[name]
                    else:
                        cls = cls.__bases__[0]
                objname,newname = splitName(name)
                if newname is not None and objname in self._value:
                    return self._value[objname].qualifiers(newname, custom=custom, default=default)
                else:
                    return NotImplemented
            else:
                return NotImplemented
        ret = {}
        if default:
            clsList = baseClassList(self.__class__)
            for cls in clsList:
                ret.update(cls.QUALIFIERS)
        if custom:
            ret.update(self.__dict__['_qualifiers'])
        if contentQualifiers:
            for key in self.dataOrder():
                if self._value[key] is not None:
                    objQualifiers = self._value[key].qualifiers(default=default, custom=custom)
                    if len(objQualifiers) > 0:
                        ret[key] = objQualifiers
        return ret

    def qualifiersDefinition(self, name=None):
        if name is None:
            defn = {}
            for cls in baseClassList(self.__class__):
                defn.update(cls.QUALIFIERS_DEFINITION)
            return defn
        else:
            cls = self.__class__
            while issubclass(cls,CData):
                if name in cls.QUALIFIERS_DEFINITION:
                    return cls.QUALIFIERS_DEFINITION[name]
                else:
                    cls = cls.__bases__[0]
            objname, newname = splitName(name)
            if newname is not None and objname in self._value:
                return self._value[objname].qualifiersDefinition(newname)
            else:
                if self.__dict__.get('_subItemObject', None) is not None:
                    return self.__dict__['_subItemObject'].qualifiersDefinition(name)
                else:
                    return {}

    def qualifiersOrder(self):
        order = []
        for cls in baseClassList(self.__class__):
            for item in cls.QUALIFIERS_ORDER:
                if not order.count(item):
                    order.append(item)
        return order

    def qualifiersEtree(self, customOnly=True, tag=None, recurse=True):
        error = CErrorReport()
        if tag is None:
            tag = self.objectName()
        if tag is None or len(tag) == 0:
            tag = self.__class__.__name__
            error.append(self.__class__, 14, name=self.objectPath())
        element = ET.Element(str(tag))
        if customOnly:
            qualis = self.qualifiers(default=False, contentQualifiers=False)
        else:
            qualis = self.qualifiers(contentQualifiers=False)
        for key, value in list(qualis.items()):
            try:
                ele = ET.Element(key)
                if self.qualifiersDefinition(key).get('type', str) == list:
                    txt = ''
                    for item in value:
                        txt = txt + str(item) + ','
                    if len(txt) > 0:
                        ele.text = txt[0:-1]
                elif self.qualifiersDefinition(key).get('type', str) == dict:
                    ele, err = self.dictQualifersEtree(ele, value)
                    error.extend(err)
                else:
                    ele.text = str(value)
                element.append(ele)
            except:
                error.append(self.__class__, 15, 'Object: ' + str(self.objectName()) + ' qualifier: ' + key)
        # Get qualifiers for content items
        if recurse:
            for key in self.dataOrder():
                if self._value[key] is not None:
                    ele,err = self._value[key].qualifiersEtree(customOnly=customOnly)
                    if len(ele) > 0:
                        element.append(ele)
                    if len(err) > 0:
                        error.extend(err)
        return element, error

    def dictQualifersEtree(self, root, inputDict):
        error = CErrorReport()
        for key,value in list(inputDict.items()):
            root.append(ET.Element(key))
            if isinstance(value,list):
                txt = ''
                for item in value:
                    txt = txt + str(item) + ','
                if len(txt) > 0:
                    root[-1].text = txt[0:-1]
            elif isinstance(value, dict):
                err,root = self.dictQualifersEtree(root, value)
                error.extend(err)
            else:
                root[-1].text = str(value)
        return root, error

    def qualifiersXmlText(self, pretty_print=True, xml_declaration=False):
        tree, error = self.qualifiersEtree()
        if pretty_print:
            ET.indent(tree)
        text = ET.tostring(tree, xml_declaration=xml_declaration)
        return text

    def setQualifiersEtree(self, element):
        rv = CErrorReport()
        qualifiers = {}
        for ele in element:
            name = str(ele.tag)
            value = str(ele.text)
            if value == 'None':
                value = None
            qualifierType = self.qualifiersDefinition(name).get('type', self.PYTHONTYPE)
            if qualifierType in [int, float, bool, str]:
                if value is None:
                    qualifiers[name] = None
                elif qualifierType == bool:
                    if ['0', 'False', 'false'].count(value):
                        qualifiers[name] = False
                    else:
                        qualifiers[name] = True
                else:
                    try:
                        qualifiers[name] = qualifierType(value)
                    except:
                        rv.append(self.__class__, 8, 'Error interpreting qualifier: ' + name, name=self.objectPath())
            elif qualifierType is list:
                strList = str(ele.text).split(',')
                qList = []
                listItemType = self.qualifiersDefinition(name).get('listItemType', self.qualifierListType())
                for item in strList:
                    try:
                        if item == 'None':
                            qList.append(None)
                        else:
                            qList.append(listItemType(item))
                    except:
                        rv.append(self.__class__, 8, 'Error interpreting qualifier: ' + name, name=self.objectPath())
                qualifiers[name] = qList
            elif qualifierType is dict:
                qualifiers[name] = self.eTreeToDict(ele)
            elif qualifierType is type:
                from core import CCP4DataManager
                cls = CCP4DataManager.DATAMANAGER().getClass(value)
                if cls is None:
                    rv.append(self.__class__, 13, 'Qualifier: ' + name, name=self.objectPath())
                else:
                    qualifiers[name] = cls
            elif name in self.contents():
                rv.extend(self._value[name].setQualifiersEtree(ele))
            elif self.__dict__.get('_subItemObject',None) is not None:
                self.__dict__['_subItemObject'].setQualifiersEtree(ele)
            else:
                rv.append(self.__class__, 7, 'Unidentified qualifier: ' + name, name=self.objectPath())
        if len(qualifiers)>0:
            try:
                self.setQualifiers(qualifiers=qualifiers, validateDefault=False)
            except CException as e:
                rv.extend(e)
            except Exception as e:
                rv.append(self.__class__, 21, name=self.objectPath())
        if self.__dict__.get('_subItemObject', None) is not None:
            self.__dict__['_subItemObject'].setQualifiersEtree(element)
        return rv

    def eTreeToDict(self, root):
        ret = {}
        for ele in root:
            if len(ele) > 0:
                ret[ele.tag] = self.eTreeToDict(ele)
            else:
                ret[ele.tag] = ele.text
        return ret

    def pythonType(self):
        return self.PYTHONTYPE

    def getDataByKey(self, keyName):
        # keyName is name of a qualifier which contains the name of another data object in the same container
        dataObjName = self.qualifiers(keyName)
        if dataObjName is None:
            return None
        dobj = self
        while dobj is not None and isinstance(dobj, CData):
            if isinstance(dobj.__dict__['_value'], dict) and \
                                dataObjName in dobj.__dict__['_value']:
                return dobj.__dict__['_value'][dataObjName]
            else:
                dobj = dobj.parent()
        return None


# A base class for all data classes
class CData(CObject, CDataQualifiers):
    VERSION = '0.0'
    CONTENTS = {}
    CONTENTS_ORDER = []
    PROPERTIES = {}
    ERROR_CODES = {0 : {'severity' : SEVERITY_OK, 'description' : 'OK'},
                   1 : {'severity' : SEVERITY_UNDEFINED, 'description' : 'Data has undefined value'},
                   2 : {'severity' : SEVERITY_UNDEFINED_ERROR, 'description' : 'Data has undefined value'},
                   3 : {'severity' : SEVERITY_WARNING, 'description' : 'Missing data'},
                   4 : {'description' : 'Missing data'},
                   5 : {'description' : 'Attempting to set data of wrong type'},
                   6 : {'description' : 'Default value does not satisfy validity check'},
                   7 : {'severity' : SEVERITY_WARNING, 'description' : 'Unrecognised qualifier in data input'},
                   8 : {'severity' : SEVERITY_WARNING, 'description' : 'Attempting to get inaccessible attribute:'},
                   9 : {'description' : 'Failed to get property'},
                   10 : {'severity' : SEVERITY_WARNING, 'description' : 'Attempting to set inaccessible attribute:'},
                   11 : {'description' : 'Failed to set property:'},
                   12 : {'description' : 'Undetermined error setting value from XML'},
                   13 : {'description' : 'Unrecognised class name in qualifier'},
                   14 : {'severity' : SEVERITY_WARNING, 'description' : 'No object name when saving qualifiers to XML'},
                   15 : {'description' : 'Error saving qualifier to XML'},
                   16 : {'severity' : SEVERITY_WARNING, 'description' : 'Unrecognised item in XML data file'},
                   17 : {'description' : 'Attempting to set unrecognised qualifier'},
                   18 : {'description' : 'Attempting to set qualifier with wrong type'},
                   19 : {'description' : 'Attempting to set qualifier with wrong list item type'},
                   20 : {'description' : 'Error creating a list/dict item object'},
                   21 : {'description' : 'Unknown error setting qualifiers from Xml file'},
                   22 : {'description' : 'Unknown error testing validity'},
                   23 : {'description' : 'Error saving data object to XML'},
                   24 : {'description' : 'Unable to test validity of default','severity' : SEVERITY_WARNING},
                   300 : {'description' : 'Compared objects are the same','severity' : SEVERITY_OK},
                   315 : {'description' : 'Both compared objects are null','severity' : SEVERITY_OK},
                   301 : {'description' : 'Unable to compare this class of data','severity' : SEVERITY_WARNING},
                   302 : {'description' : 'Other data has null value'},
                   303 : {'description' : 'My data has null value'},
                   304 : {'description' : 'Data has different values'}}

    def __init__(self, value=[], qualifiers={}, parent=None, name=None, build=True, **kw):
        CObject.__init__(self, parent=parent, name=name)
        #Sort the 'loose' arguments to values or qualifiers
        qualis = {}
        qualis.update(qualifiers)
        if len(kw) > 0:
            vs, qs, us = sortArguments(self.__class__, kw)
            if len(us) > 0:
                for key in list(us.keys()):
                    print(self.className(), 'input argument unrecognised:', key)
            qualis.update(qs)
        else:
            vs = {}
        # Beware sub-classes that need some qualifiers for build may have
        # already initialised _qualifiers
        if '_qualifiers' not in self.__dict__:
            self.__dict__['_qualifiers'] = {}
        self.__dict__['_value'] = {}
        if len(qualis) > 0:
            self.setQualifiers(qualis)
        if not build:
            return
        self.build(qualifiers=qualis)
        self.setDefault()
        # Using the variable valu to overcome weird bug that a second
        # instantiation of class had values that had been passed to first
        # instance via the **kw mechanism
        if isinstance(value, CData):
            self.set(value=value)
        else:
            valu = {}
            if isinstance(value, dict):
                valu.update(value)
            if len(vs) > 0:
                valu.update(vs)
            if len(valu) > 0:
                self.set(value=valu)

    def build(self, qualifiers={}):
        for key, defn in list(self.contents().items()):
            qualis = defn.get('qualifiers',{})
            try:
                if key in qualifiers:
                    qualis.update(qualifiers[key])
            except:
                print('Error in CData.build', qualifiers, key, qualis)
            default = qualifiers.get('default', {}).get(key, None)
            if default is not None:
                qualis['default'] = default
            # Nov13 - the following line was previously uncommented - but it does not make sense
            if 'subItem' in defn:
                self.buildItem(key=key, cls=defn['class'], qualifiers=qualis, subItem=defn['subItem'])
            else:
                self.buildItem(key=key, cls=defn['class'], qualifiers=qualis)

    def buildItem(self, key=None, cls=None, qualifiers={}, subItem={}, emitDataChanged=True):
        if len(subItem) > 0:
            self.__dict__['_value'][key] = cls(qualifiers=qualifiers, parent=self, name=key, subItem=subItem)
        else:
            self.__dict__['_value'][key] = cls(qualifiers=qualifiers, parent=self, name=key)
        if emitDataChanged:
            try:
                self.__dict__['_value'][key].dataChanged.connect(self.dataChanged.emit)
            except:
                raise
                print("Fail")

    def objectPath(self, ifContainer=True):
        if not ifContainer:
            from core import CCP4Container
        obj = self
        path = ''
        sep = ''
        while obj is not None and isinstance(obj, CData) and (ifContainer or not isinstance(obj, CCP4Container.CContainer)):
            name = obj.objectName()
            if isinstance(obj.parent(),CList) and name != 'subItemObject':
                try:
                    indx = str(obj.parent().index(obj))
                except:
                    indx = '?'
                path = '[' + indx + ']' + sep + path
                sep = ''
            else:
                path = obj.objectName() + sep + path
                sep = '.'
            obj = obj.parent()
        return path

    def getVersion(self):
        return self.VERSION

    PROPERTIES = {'version' : {'fget' : getVersion}}

    def Property(self, name=None, key='fget'):
        cls = self.__class__
        while issubclass(cls, CData):
            if hasattr(cls, 'PROPERTIES') and name in cls.PROPERTIES:
                if key in cls.PROPERTIES[name]:
                    return cls.PROPERTIES[name][key]
                else:
                    return None
            else:
                cls = cls.__bases__[0]
        return None

    def guiLabel(self):
        label = self.qualifiers('guiLabel')
        if label is not None and label is not NotImplemented:
            return label
        else:
            return str(self.objectName())

    def contents(self, name=None):
        if name is None:
            return self.__class__.CONTENTS
        elif name in self.__class__.CONTENTS:
            return self.__class__.CONTENTS[name]
        else:
            return None

    def dataOrder(self):
        if len(self.__class__.CONTENTS_ORDER) == len(self.__dict__['_value']):
            return self.__class__.CONTENTS_ORDER
        else:
            keyList = list(self.__dict__['_value'].keys())
            return keyList

    def setDefault(self, force=False):
        for name, obj in list(self.__dict__['_value'].items()):
            if obj is not None:
                if hasattr(obj,"setDefault"):
                    obj.setDefault(force=force)

    def validity(self, arg={}):
        # Reimplement in subclasses that have some dependency between
        # the different items of data
        return self.itemValidity(arg)

    def itemValidity(self, arg={}, keys=None):
        validityObj = CErrorReport()
        # test validity of all items in a complex data object
        if isinstance(arg, CData):
            testData = arg.get()
        else:
            testData = {}
            testData.update(arg)
        if keys is None:
            keys = list(self.contents().keys())
        for key in keys:
            if key in testData:
                if key in self._value:
                    try:
                        testValidity = self._value[key].validity(testData[key])
                    except:
                        print('ERROR testing', key, 'has value:', testData[key], 'of type:', type(testData[key]))
                        testValidity = CErrorReport(self.__class__, 22, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                else:
                    print('ERROR unable to test because item not built',key)
                    testValidity = CException(self.__class__, 24, name=self.objectPath(),label=self.qualifiers('guiLabel'), stack=False)
                validityObj.extend(testValidity, label=key)
        return validityObj

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        if allSet:
            for key in list(self.contents().keys()):
                if not self._value[key].isSet(allowUndefined=allowUndefined, allowDefault=allowDefault):
                    return False
            return True
        else:
            for key in list(self.contents().keys()):
                try:
                    if self._value[key].isSet(allowUndefined=allowUndefined, allowDefault=allowDefault, allSet=allSet):
                        return True
                except:
                    print('ERROR in CData.isSet', self.objectPath(), list(self.contents().keys()), list(self._value.keys()))
            return False

    def fix(self, arg=None):
        # Apply possible fixes to data that has failed validity test
        # This should be re-implemented in derived classes
        return arg

    def coerce(self, arg):
        # Apply possible coersion to input data
        # This should be re-implemented in derived classes
        return arg

    def __getattr__(self, name):
        if name in self.__dict__['_value']:
            return self.__dict__['_value'][name]
        else:
            propFget = self.Property(name, 'fget')
            if propFget is not None:
                try:
                    return propFget(self)
                except Exception as e:
                    raise AttributeError("%r object has no attribute %r" % (self.__class__.__name__, name))
                    """
                    raise CException(CData, 9, details=name + '\n' + str(e), exc_info=sys.exc_info())
                    print('Inaccessible attribute',name,str(e))
                    traceback.print_stack(limit=5)
                    """
            else:
                raise AttributeError("%r object has no attribute %r" % (self.__class__.__name__, name))
                """
                raise CException(CData, 8, name)
                print('Inaccessible attribute', name, str(e))
                traceback.print_stack(limit=5)
                """

    def __setattr__(self, name, value):
        if '_dataOrder' in self.__dict__ and name in self.dataOrder():
            if hasattr(self.__dict__['_value'][name],"set"):
                return self.__dict__['_value'][name].set(value)
        else:
            propFset = self.Property(name, 'fset')
            if propFset is not None:
                return propFset(self, value)
            else:
                pass
                #raise CException(CData, 10, name)

    def __str__(self):
        d = {}
        for item in self.dataOrder():
            d[item] = self.__dict__['_value'][item].__str__()
        return d.__str__()

    def getTextItem(self):
        '''Return a str value to display on the gui'''
        if self.isSet():
            return self.__str__()
        else:
            return '--'

    def getTableTextItems(self):
        '''Return list of str values to appear in a list widget'''
        items = []
        for key in self.dataOrder():
            if self.__dict__['_value'][key].isSet():
                items.append( self.__dict__['_value'][key].__str__())
            else:
                items.append('')
        return items

    def orderedDictText(self):
        # return  dict of all data
        ret = '{'
        for key in self.dataOrder():
            #Everything in self._value should be a CData..
            value = self.__dict__['_value'][key].__str__()
            ret = ret + "'" + key + "':" + value + ","
        if ret[-1] == ',':
            ret = ret[0:-1] + '}'
        else:
            ret = ret + '}'
        return ret

    def xmlText(self, pretty_print=True, xml_declaration=False):
        tree = self.getEtree()
        if pretty_print:
            ET.indent(tree)
        text = ET.tostring(tree, xml_declaration=xml_declaration)
        return text

    def __len__(self):
        return self.__dict__['_value'].__len__()

    def set(self, value={}, **kw):
        myException = None
        saveValue = self.get0()
        if isinstance(value, self.__class__):
            for key in self.contents():
                try:
                    other = value.get(key)
                    self._value[key].set(other)
                except CException as e:
                    if myException is None:
                        myException = e
                    else:
                        myException.extend(e)
        else:
            valueIn = {}
            if len(value) > 0:
                valueIn.update(value)
            if len(kw) > 0:
                valueIn.update(kw)
            validityObj = self.copyValue(self.coerce(valueIn))
            if validityObj.maxSeverity() < SEVERITY_ERROR:
                validityObj = self.validity(self._value)
            # Validity check failed - restore original
        self.updateData()

    def copyValue(self,value):
        ''' Copy data without validity check '''
        rv = CErrorReport()
        for key, val in list(value.items()):
            objname, newkey = splitName(key)
            if newkey is not None:
                obj = self._value.get(objname, None)
                if obj is None:
                    print('error interpreting value:', key, 'unknown object:', objname)
                else:
                    print('CData.copyValue newkey', obj, newkey,val)
                    try:
                        obj.set(value={newkey : val})
                    except CException as e:
                        rv.extend(e)
            else:
                if key in self.contents():
                    try:
                        self.__dict__['_value'][key].set(val)
                    except CException as e:
                        rv.extend(e)
        return rv

    # Reimplement in sub-class if require something updating after
    # calling set() or setEtee (NB for complex classes setEtree just
    # calls the component class setEtree and does not call set)
    def updateData(self):
        self.dataChanged.emit()

    def get(self, name=None, default=None):
        # If input name is set then return item of that name
        if name is not None:
            if name in self.__dict__['_value']:
                return self.__dict__['_value'].get(name)
            else:
                objname, newname = splitName(name)
                if newname is None:
                    return default
                elif objname in self._value:
                    return self._value[objname].get(newname)
                #Try possibility that first item in name is actually reference to self
                elif objname == str(self.objectName()):
                    objname,newname = splitName(newname)
                    if newname is None:
                        return default
                    elif objname in self._value:
                        return self._value[objname].get(newname)
        # return 'flattened' dict of all data
        ret = {}
        for key in self.dataOrder():
            #Everything in self._value should be a CData.
            value = self.__dict__['_value'].get(key,None)
            if isinstance(value,CData):
                val = value.get()
                ret[key] = val
            else:
                pass
        return ret

    def get0(self, name=None):
        # If input name is set then return item of that name
        if name is not None:
            if name in self.contents():
                return self._value.get(name)
            else:
                if name in self._value:
                    return self._value[name]
        # return nested dict of all data
        ret = {}
        for key, value in list(self.__dict__['_value'].items()):
            #Everything in self._value should be a CData..
            if isinstance(value, CData):
                val = value.get()
                ret[key] = val
            else:
                pass
        return ret

    def getDataObjects(self, name=None):
        if name is not None:
            if name in self.contents():
                return self._value[name]
            else:
                return None
        else:
            ret = {}
            for key in list(self.contents().keys()):
                ret[key] = self._value[key]
            return ret

    def unSet(self):
        self.blockSignals(True)
        for key in list(self.contents().keys()):
            self._value[key].unSet()
        self.blockSignals(False)
        self.updateData()

    def getEtree(self, excludeUnset=True, name=None, useLXML=True):
        if name is None:
            name = self.objectName()
        if name is None or len(name) == 0:
            name = self.className()
        if useLXML:
            element = ET.Element(name)
        else:
            element = ET.Element(name)
        for key in self.dataOrder():
            if not excludeUnset or self._value[key].isSet(allSet=False):
                ele = self._value[key].getEtree(useLXML=useLXML)
                element.append(ele)
        return element

    def setEtree(self, element, checkValidity=True):
        rv = CErrorReport()
        for ele in element:
            name = ele.tag
            if name in self.dataOrder():
                try:
                    self._value[name].setEtree(ele, checkValidity=checkValidity)
                except CException as e:
                    rv.extend(e)
                except:
                    rv.append(self.__class__, 12, name)
            else:
                eleId = str(ele.get('id'))
                rv.append(self.__class__, 16, 'Item: ' + name + ' id: ' + eleId, name=self.objectPath())
        self.updateData()
        return rv

    def saveDataToXml(self, fileName=None):
        errorReport = CErrorReport()
        from core import CCP4File
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        if getattr(self, 'header', None) is not None:
            f.header.set(self.header)
        f.header.function.set(self.__class__.__name__)
        f.header.setCurrent()
        bodyEtree = self.getEtree()
        try:
            f.saveFile(bodyEtree=bodyEtree)
        except CException as e:
            errorReport.append(e)
        except:
            errorReport.append(self.__class__, 23, fileName)
        return errorReport


    def blockSignals(self,mode):
        # Propagate blockSignals to component objects
        CObject.blockSignals(self,mode)
        for obj in list(self._value.values()):
            obj.blockSignals(mode)

    def getTableTextItems(self):  # KJS : Flagged as duplicate. Look into this.
        return [self.__str__()]

    def saveToDb(self):
        return [], None, {}

    def testComparisonData(self):
        return self.saveToDb()

    def parentContainer(self):
        from core import CCP4Container
        try:
            obj = self
            while isinstance(obj, CData):
                obj = obj.parent()
                if isinstance(obj, CCP4Container.CContainer):
                    return obj
            return None
        except:
            return None

    def assertSame(self, arg=None):
        return CErrorReport(self.__class__, 301, name=self.objectPath(), details=str(self) + ' : ' + str(arg))

    def removeUnsetListItems(self):
        # Go through child items looking for lists and remove any
        # list items that are unset
        for key in self.dataOrder():
            try:
                self.__dict__['_value'][key].removeUnsetListItems()
            except:
                print('Error in CData.removeUnsetListItems', self.objectPath())

    def locateElement(self, objectPath):
        baseElement = self
        pathElements = objectPath.split('.')
        arrayFinder = re.compile(r'(?P<baseElement>.*)\[(?P<index>\d+)\]$')
        for pathElement in pathElements:
            print('pathElement', pathElement)
            matches = arrayFinder.match(pathElement)
            if matches is not None or isinstance(baseElement, CList):
                if matches is not None:
                    elementList = getattr(
                        baseElement, matches.groupdict()['baseElement'])
                    elementIndex = int(matches.groupdict()['index'])
                elif isinstance(baseElement, CList):
                    elementList = getattr(baseElement, pathElement)
                    elementIndex = 0
                while len(elementList) <= elementIndex:
                    elementList.append(elementList.makeItem())
                try:
                    baseElement = elementList[elementIndex]
                except TypeError as err:
                    print(f'Failed deindexing {objectPath} at {pathElement}')
                    raise err
            else:
                baseElement = getattr(baseElement, pathElement)
        return baseElement
'''
A base class for classes that are equivalent to Python types.
This implies that the CONTENT definition is empty.
'''

class CBaseData(CData):
    '''Base class for simple classes'''

    PYTHONTYPE = str
    QUALIFIERS = {'default' : NotImplemented, 'charWidth' : 10}
    QUALIFIERS_DEFINITION = {'charWidth' : {'type' : int}}
    QUALIFIERS_ORDER = ['charWidth']

    def __init__(self, value=NotImplemented, qualifiers={}, parent=None, name=None,**kw):
        #Sort the 'loose' arguments to values or qualifiers
        # Do limited functions from CData.__init__
        CObject.__init__(self, parent=parent, name=name)
        if '_qualifiers' not in self.__dict__:
            self.__dict__['_qualifiers'] = {}
        self.build()
        qualis = {}
        qualis.update(qualifiers)
        if len(kw) > 0:
            vs, qs, us = sortArguments(self.__class__, kw)
            if len(us) > 0:
                for key in list(us.keys()):
                    print(self.className(), 'input argument unrecognised:', key)
            qualis.update(qs)
        self.__dict__['_qualifiers'] = {}
        if len(qualis) > 0:
            self.setQualifiers(qualis)
        self.setDefault()
        if value is not NotImplemented:
            self.set(value=value)

    def validity(self, arg=None):
        e = CErrorReport()
        if arg is None:
            e.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return e

    def coerce(self, arg):
        if arg is None:
            return arg
        elif isinstance(arg, self.__class__):
            arg = arg.get()
        if isinstance(arg, self.PYTHONTYPE) or arg is None:
            pass
        elif arg == '':
            arg = None
        else:
            try:
                arg = self.PYTHONTYPE(arg)
            except:
                pass
        return arg

    def getValue(self, arg):
        try:
            # MN a NoneType value cannot be coerced to the corresponding type, hence return None
            if arg.get() is None:
                return None
            return self.pythonType()(arg.get())
        except:
            return arg

    def build(self, qualifiers={}):
        default = self.qualifiers('default')
        if default is not NotImplemented:
            self.__dict__['_value'] = default
        else:
            self.__dict__['_value'] = None

    def set(self, value=None, checkValidity=True):
        value=self.coerce(value)
        if checkValidity:
            v = self.validity(value)
            if v.maxSeverity() < SEVERITY_ERROR:
                if self.__dict__['_value'] != value:
                    self.__dict__['_value'] = value
                    self.updateData()
            else:
                e = CException()
                e.extend(v)
                raise e
        else:
            if self.__dict__['_value'] != value:
                self.__dict__['_value'] = value
                self.updateData()

    def get(self):
        return self.__dict__['_value']

    def get0(self):   # KJS : WTH ?!? What does this 0 even mean ?
        return self.__dict__['_value']

    def setDefault(self, force=False):
        default = self.qualifiers('default')
        if default is NotImplemented:
            if force:
                self.__dict__['_value'] = None
            else:
                return
        else:
            self.__dict__['_value'] = default

    def __str__(self):
        return self.__dict__['_value'].__str__()

    def __getattr__(self, name=None):
        if name is None:
            return self.__dict__['_value']
        else:
            propFget = self.Property(name, 'fget')
            if propFget is not None:
                try:
                    return propFget(self)
                except:
                    raise AttributeError("%r object has no attribute %r" % (self.__class__.__name__, name))
                    """
                    raise CException(CData, 9, name)
                    print('Inaccessible attribute')
                    traceback.print_stack(limit=5)
                    """
            else:
                raise AttributeError("%r object has no attribute %r" % (self.__class__.__name__, name))
                """
                raise CException(CData, 8, name)
                print('Inaccessible attribute')
                traceback.print_stack(limit=5)
                """

    def __setattr__(self, name=None, value=NotImplemented):
        if value is NotImplemented:
            # The 'name' is actually the value that we want to set
            return self.set(name)
        else:
            propFset = self.Property(name, 'fset')
            if propFset is not None:
                return propFset(self, value)
            else:
                pass
                #raise CException(CData, 10, name)


    def isDefault(self):
        default = self.qualifiers('default')
        if default is None and self._value is None:
            return True
        elif  default is None or self._value is None:
            return False
        else:
            return self._value == default

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        # test whether the data value is set
        # If allowUndefined true and  qualifier allowUndefined true then return 'set'
        # If the data has the default value and  allowDefault==False then return as 'unset'
        # allSet is unused - for consistency with CData.isSet()
        if self._value is None or self._value is NotImplemented:
            if not allowUndefined or not self.qualifiers('allowUndefined'):
                return False
            else:
                return True
        else:
            if not allowDefault and self._value == self.qualifiers('default'):
                return False
        return True

    def unSet(self):
        value = self.qualifiers('default')
        if value is NotImplemented:
            self.set(value=None)
        else:
            self.set(value=self.qualifiers('default'))

    def getEtree(self, excludeUnset=True, name=None, useLXML=True):
        if name is None:
            name = self.objectName()
        if name is None or len(name) == 0:
            name = self.className()
        if useLXML:
            element = ET.Element(str(name))
        else:
            element = ET.Element(str(name))
        if hasattr(self.__str__(),"decode"):
            text = self.__str__().decode('ISO-8859-1')
        else:
            text = self.__str__()
        if text not in ['None', 'NotImplemented']:
            element.text = text
        return element

    def setEtree(self, element, checkValidity=True):
        rv = CErrorReport()
        if element.text is None:
            value = None
        else:
            value = str(element.text)
        try:
            self.set(value, checkValidity=checkValidity)
        except CException as e:
            rv.extend(e)
        except:
            rv.append(self.__class__, 12, element.tag, name=self.objectPath())
        return rv

    def getPersistent(self):
        from CCP4Persistent import CPersistentBaseData   # KJS : Check this is in there ok !
        p = CPersistentBaseData(name=str(self.objectName()), className = self.__class__.__name__,
                                version = self.version, value = self.get())
        return p

    def setPersistent(self, p):
        # Should be checking version !!!
        self.setObjectName(p.name)
        self.set(p.value)

    def dataOrder(self):
        return []

    def blockSignals(self, mode):
        # Override CData.blockSignals that propagates blockSignals to component objects
        CObject.blockSignals(self, mode)

    def assertSame(self, arg):
        if isinstance(arg, self.__class__):
            other = arg.__dict__['_value']
        else:
            other = arg

        if self.__dict__['_value'] is None:
            if other is None:
                return CErrorReport()
            else:
                return CErrorReport(self.__class__, 303, name=self.objectPath(), details=str(self) + ' : ' + str(other))
        elif other is None:
            return CErrorReport(self.__class__, 302, name=self.objectPath(), details=str(self) + ' : ' + str(other))
        if other == self.__dict__['_value']:
            return CErrorReport()
        else:
            return CErrorReport(self.__class__, 304, name=self.objectPath(), details=str(self) + ' : ' + str(other))

    def removeUnsetListItems(self):
        return


class CInt(CBaseData):
    '''An integer'''
    # Qualifiers determine default values and limits on allowed values.
    # They allow some customisation without requiring writing a new Python class.
    # A data type could be defined in an external file by specifying the appropriate class and qualifiers.
    CONTENTS = {}
    PYTHONTYPE = int
    QUALIFIERS = {'default' : NotImplemented, 'max' : None, 'min' : None,
                  'enumerators' : [], 'menuText' : [], 'onlyEnumerators' : False,}
    QUALIFIERS_ORDER = ['min', 'max', 'onlyEnumerators', 'enumerators', 'menuText']
    QUALIFIERS_DEFINITION = {'default' : {'type' :int},
                             'max' : {'type' :int, 'description' : 'The inclusive minimum allowed value'},
                             'min' : {'type' :int,'description' : 'The inclusive maximum allowed value'},
                             'enumerators' : {'type' :list, 'listItemType' :int,
                                              'description' : 'A Python list of allowed or recommended values - see onlyEnumerators'},
                             'menuText' : {'type' : list, 'listItemType' : str,
                                           'description' : 'A Python list of strings, matching items in enumerators list, to appear on GUI menu'},
                             'onlyEnumerators': {'type' : bool,
                                                 'description' : 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'}}
    ERROR_CODES = {101 : {'description' : 'below minimum'},
                   102 : {'description' : 'above maximum'},
                   103 : {'description' : 'not one of limited allowed values'}}

    def validity(self, arg):
        validityObj = CErrorReport()
        if arg is None or (isinstance(arg,CInt) and arg.__dict__['_value'] is None ):
            if self.qualifiers('allowUndefined'):
                validityObj.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                validityObj.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return validityObj
        if not isinstance(arg, (int, CInt)):
            validityObj.append(self.__class__, 5, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return validityObj
        if self.qualifiers('min') is not None and self.qualifiers('min') is not NotImplemented and arg < self.qualifiers('min'):
            validityObj.append(self.__class__, 101, name=self.objectPath(), details='of ' + str(self.qualifiers('min')), label=self.qualifiers('guiLabel'), stack=False)
        if self.qualifiers('max') is not None and self.qualifiers('max') is not NotImplemented and arg > self.qualifiers('max'):
            validityObj.append(self.__class__, 102, name=self.objectPath(), details='of ' + str(self.qualifiers('max')), label=self.qualifiers('guiLabel'), stack=False)
        if self.qualifiers('onlyEnumerators') and len(self.qualifiers('enumerators')) > 0 and self.qualifiers('enumerators').count(arg) < 1:
            validityObj.append(self.__class__, 103, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return validityObj

    def getMenuValue(self):
        if 'menuText' in self._qualifiers and len(self._qualifiers['menuText']) > 0 and self._value in self._qualifiers['enumerators']:
            return self._qualifiers['menuText'][self._qualifiers['enumerators'].index(self._value)]
        else:
            return self._value

    # Applying mathematical functions to self._value in place
    # These methods return a CErrorReport
    def add(self, arg):
        if self._value is not None:
            self.set(self._value + self.getValue(arg))

    def sub(self, arg):
        if self._value is not None:
            self.set(self._value - self.getValue(arg))

    def mul(self, arg):
        if self._value is not None:
            self.set(self._value * self.getValue(arg))

    def div(self, arg):
        if self._value is not None:
            self.set(self._value / self.getValue(arg))

    def abs(self):
        if self._value is not None:
            self.set(abs(self._value))

    def pow(self, arg):
        if self._value is not None:
            self.set(self._value.__pow__(self.getValue(arg)))

    def neg(self, arg):
        if self._value is not None:
            self.set((self._value.__neg__()))

    # Implementation of the usual methods for int
    def __abs__(self):
        return self.__class__(value=self._value.__abs__(),qualifiers=self._qualifiers)

    def __add__(self, arg):
        return self.__class__(value=self._value.__add__(self.getValue(arg)),qualifiers=self._qualifiers)

    def __and__(self, arg):
        return self._value.__and__(self.getValue(arg))

    def __round__(self):
        return round(self._value)

    def __cmp__(self, arg):
        # Will throw exception if self.value not set
        #MN Allow for equality of both self._value and other._value are None
        if self._value is None and self.getValue(arg) is None:
            return 0
        return self._value.__cmp__(self.getValue(arg))

#These methods are required by Python3 which does not honour __cmp__ anymore : https://portingguide.readthedocs.io/en/latest/comparisons.html
    def __eq__(self, arg):
        if self._value is None and self.getValue(arg) is None:
            return 0
        return self._value == self.getValue(arg)

    def __ne__(self, arg):
        if self._value is None and self.getValue(arg) is not None:
            return 1
        if self._value is not None and self.getValue(arg) is None:
            return 1
        return self._value != self.getValue(arg)

    def __cmp__(self, arg):
        # Will throw exception if self.value not set
        #MN Allow for equality of both self._value and other._value are None
        if self._value is None and self.getValue(arg) is None:
            return 0
        return self._value.__cmp__(self.getValue(arg))

#These methods are required by Python3 which does not honour __cmp__ anymore : https://portingguide.readthedocs.io/en/latest/comparisons.html
    def __eq__(self, arg):
        if self._value is None and self.getValue(arg) is None:
            return 0
        return self._value == self.getValue(arg)

    def __ne__(self, arg):
        if self._value is None and self.getValue(arg) is not None:
            return 1
        if self._value is not None and self.getValue(arg) is None:
            return 1
        return self._value != self.getValue(arg)

    def __lt__(self, arg):
        if self._value is None or self.getValue(arg) is None:
            raise TypeError("'<' not supported between instances of '"+str(type(self._value))+"' and '"+str(type(self.getValue(arg)))+"'")
        return self._value < self.getValue(arg)

    def __le__(self, arg):
        if self._value is None or self.getValue(arg) is None:
            raise TypeError("'<=' not supported between instances of '"+str(type(self._value))+"' and '"+str(type(self.getValue(arg)))+"'")
        return self._value <= self.getValue(arg)

    def __gt__(self, arg):
        if self._value is None or self.getValue(arg) is None:
            raise TypeError("'>' not supported between instances of '"+str(type(self._value))+"' and '"+str(type(self.getValue(arg)))+"'")
        return self._value > self.getValue(arg)

    def __ge__(self, arg):
        if self._value is None or self.getValue(arg) is None:
            raise TypeError("'>=' not supported between instances of '"+str(type(self._value))+"' and '"+str(type(self.getValue(arg)))+"'")
        return self._value >= self.getValue(arg)

    def __coerce__(self, arg):
        return NotImplemented

    def __div__(self, arg):
        return self.__class__(value=self._value.__div__(self.pythonType()(arg)), qualifiers=self._qualifiers)

    def __divmod__(self, arg):
        return self._value.__divmod__(self.getValue(arg))

    def __float__(self):
        return self._value.__float__()

    def __floordiv__(self, arg):
        return self.__class__(value=self._value.__floordiv__(self.getValue(arg)), qualifiers=self._qualifiers)

    def  __format__(self, arg):
        return self._value.__format__(arg)

    def __hex__(self):
        return self._value.__hex__()

    def __index__(self):
        return self._value.__index__()

    def __int__(self):
        return self._value.__int__()

    def __invert__(self):
        return self.__class__(value=self._value.__invert__(), qualifiers=self._qualifiers)

    def __long__(self):
        return self._value.__long__()

    def __lshift__(self):
        print('NOT IMPEMENTED: CInt.__lshift__')

    def __mod__(self,arg):
        return self.__class__(value=self._value.__mod__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __mul__(self,arg):
        return self.__class__(value=self._value.__mul__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __neg__(self):
        return self.__class__(value=self._value.__neg__(), qualifiers=self._qualifiers)

    def __bool__(self):
        if self._value is None:
            return False
        elif self._value:
            return True
        else:
            return False

    def __nonzero__(self):
        return self.__bool__()

    def __oct__(self):
        return self._value.__oct__()

    def __or__(self,arg):
        return self._value.__or__(self.getValue(arg))

    def __pos__(self):
        return self.__class__(value=self._value.__pos__(), qualifiers=self._qualifiers)

    def __pow__(self, arg):
        return self.__class__(value=self._value.__pow__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __radd__(self, arg):
        return self.__class__(value=self._value.__radd__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rand__(self, arg):
        return self._value.__rand__(arg)

    def __rdiv__(self, arg):
        return self.__class__(value=self._value.__rdiv__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rdivmod__(self, arg):
        return self._value.__rdivmod__(self.getValue(arg))

    def __rfloordiv__(self, arg):
        return self.__class__(value=self._value.__rfloordiv__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rlshift__(self):
        print('NOT IMPEMENTED: CInt.__rlshift__')

    def __rmod__(self,arg):
        return self.__class__(value=self._value.__rmod__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rmul__(self, arg):
        return self.__class__(value=self._value.__rmul__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __ror__(self, arg):
        return self._value.__ror__(self.getValue(arg))

    def __rpow__(self, arg):
        return self.__class__(value=self._value.__rpow__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rrshift__(self):
        print('NOT IMPEMENTED: CInt.__rrshift__')

    def __rshift__(self):
        print('NOT IMPEMENTED: CInt.__rshift__')

    def __rsub__(self, arg):
        return self.__class__(value=self._value.__rsub__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rtruediv__(self, arg):
        return self._value.__rtruediv__(self.getValue(arg))

    def __rxor__(self, arg):
        return self._value.__rxor__(self.getValue(arg))

    def __str__(self):
        return self._value.__str__()

    def __sub__(self, arg):
        return self.__class__(value=self._value.__sub__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __truediv__(self, arg):
        return self._value.__truediv__(self.getValue(arg))

    def __trunc__(self):
        print('NOT IMPEMENTED: CInt.__trunc__')

    def __xor__(self, arg):
        return self._value.__xor__(self.getValue(arg))


class CString(CBaseData):
    '''A string'''
    PYTHONTYPE = str
    CONTENTS = {}
    QUALIFIERS = {'minLength' : None, 'maxLength' : None, 'enumerators' : [],
                  'menuText' : [], 'onlyEnumerators' : False,
                  'charWidth' : -1, 'allowedCharsCode' : 0}
    QUALIFIERS_ORDER = ['minLength','maxLength','onlyEnumerators','enumerators','menuText','allowedCharsCode']
    QUALIFIERS_DEFINITION = {'default' : {'type' :str},
                             'maxLength' : {'type' :int,
                                            'description' : 'Maximum length of string'},
                             'minLength' : {'type' :int,
                                            'description' : 'Minimum length of string'},
                             'enumerators' : {'type' :list,
                                              'description' : 'A list of allowed or recommended values for string'},
                             'menuText' : {'type' :list,
                                           'description' : 'A list of strings equivalent to the enumerators that will appear in the GUI'},
                             'onlyEnumerators' : {'type' : bool,
                                                  'description' : 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'},
                             'allowedCharsCode' : {'type' :int,
                                                   'description' : 'Flag if the text is limited to set of allowed characters'}}
    ERROR_CODES = {101 : {'description' : 'String too short'},
                   102 : {'description' : 'String too long'},
                   103 : {'description' : 'not one of limited allowed values'},
                   104 : {'description' : 'Contains disallowed characters'}}

    RE_PATTERN_WHITESPACE = None

    def validity(self, arg):
        validityObj = CErrorReport()
        if arg is None:
            if self.qualifiers('allowUndefined'):
                validityObj.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                validityObj.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return validityObj
        if not isinstance(arg,(str, CString)):
            validityObj.append(self.__class__, 5, name=self.objectPath(),details='Type is:' + str(type(arg)), label=self.qualifiers('guiLabel'), stack=False)
            return validityObj
        if self.qualifiers('minLength') is not None and len(arg) < self.qualifiers('minLength'):
            validityObj.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        if self.qualifiers('maxLength') is not None and len(arg) > self.qualifiers('maxLength'):
            validityObj.append(self.__class__, 102, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        if self.qualifiers('onlyEnumerators') and len(self.qualifiers('enumerators')) > 0 and self.qualifiers('enumerators').count(arg) < 1:
            validityObj.append(self.__class__, 103, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        if self.qualifiers('allowedCharsCode') > 0:
            from core import CCP4Utils
            if str(arg) != CCP4Utils.safeOneWord(str(arg)):
                validityObj.append(self.__class__, 104, details=str(arg), name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return validityObj

    def fix(self,arg):
        if self.qualifiers('allowedCharsCode') > 0:
            from core import CCP4Utils
            return CCP4Utils.safeOneWord(str(arg))

    def reWhiteSpacePattern(self):
        if CString.RE_PATTERN_WHITESPACE is None:
            import string
            pat = ''
            for item in string.whitespace:
                pat = pat + repr(item)[1:-1] + '|'
            CString.RE_PATTERN_WHITESPACE = re.compile(pat[0:-1])
        return  CString.RE_PATTERN_WHITESPACE

    def removeWhiteSpace(self,arg):
        p = self.reWhiteSpacePattern()
        arg = p.sub('',arg)
        return arg

    def getTextItem(self):
        if self.isSet():
            return self.__dict__['_value']
        else:
            return '--'

    def getMenuValue(self):
        if 'menuText' in self._qualifiers and len(self._qualifiers['menuText']) > 0 and self._value in self._qualifiers['enumerators']:
            return self._qualifiers['menuText'][self._qualifiers['enumerators'].index(self._value)]
        else:
            return self._value

    def add(self, arg):
        other = self.getValue(arg)
        return self.set(self._value.__add__(other))

    def __delitem__(self, indx):
        newString = self._value[0:indx] + self._value[indx+1:]
        return self.set(newString)

    def __delslice__(self, indx1, indx2):
        newString = self._value[0:indx1] + self._value[indx2:]
        return self.set(newString)

    def insert(self,indx,subString):
        subString = self.getValue(subString)
        newString = self._value[0:indx] + subString + self._value[indx:]
        return self.set(newString)

    def Capitalize(self):
        return self.set(self._value.capitalize())

    def Center(self, arg1, arg2=' '):
        return self.set(self._value.center(arg1, arg2))

    def Expandtabs(self, size=8):
        return self.set(self._value.expandtabs(size))

    def Ljust(self, width, fillChar=' '):
        return self.set(self._value.ljust(width, fillChar))

    def Lower(self):
        return self.set(self._value.lower())

    def Lstrip(self, subString=None):
        return self.set(self._value.lstrip(subString))

    def Join(self, arg):
        other = self.getValue(arg)
        return self.set(self._value.join(other))

    def Replace(self, arg1, arg2):
        return self.set(self._value.replace(arg1, arg2))

    def Rjust(self, width, fillChar=' '):
        return self.set(value=self._value.rjust(width, fillChar))

    def Rstrip(self, arg=None):
        return self.set(self._value.rstrip(arg))

    def Strip(self, arg=None):
        return self.set(self._value.strip(arg))

    def Swapcase(self):
        return self.set(self._value.swapcase())

    def Title(self):
        return self.set(self._value.title())

    # translate
    def Upper(self):
        return self.set(self._value.upper())

    def Zfill(self, width):
        return self.set(self._value.zfill(width))

    def __add__(self, arg):
        other = self.getValue(arg)
        return self.__class__(value=self._value.__add__(other), qualifiers=self._qualifiers)

    def __contains__(self, arg):
        other = self.getValue(arg)
        return self._value.__contains__(other)

    def __eq__(self, arg):
        other = self.getValue(arg)
        if self._value is None:
            if other is None:
                return True
            else:
                return False
        else:
            return self._value.__eq__(other)

    def __format__(self, arg):
        return self._value.__format__(arg)

    def __ge__(self, arg):
        other = self.getValue(arg)
        return self._value.__ge__(other)

    def __getitem__(self, arg):
        return self._value.__getitem__(arg)

    def __getslice__(self, arg1, arg2):
        return self._value.__getslice__(arg1, arg2)

    def __gt__(self, arg):
        other = self.getValue(arg)
        return self._value.__gt__(other)

    def __hash__(self):
        return self._value.__hash__()

    def __le__(self, arg):
        other = self.getValue(arg)
        return self._value.__le__(other)

    def __len__(self):
        if self._value is None:
            return 0
        else:
            return self._value.__len__()

    def __lt__(self, arg):
        other = self.getValue(arg)
        return self._value.__lt__(other)

    #__mod__ __mul__
    def __ne__(self,arg):
        other = self.getValue(arg)
        if self._value is None:
            if other is None:
                return False
            else:
                return True
        else:
            return self._value.__ne__(other)

    # __reduce__ __reduce_ex__
    # __rmod__ __rmul__ __setattr__
    def __sizeof__(self):
        return self._value.__sizeof__()

    def __str__(self):
        if self._value is None:
            return ''
        else:
            return self._value.__str__()

    # __subclasshook__ __formatter_field_name_split__ __formatter_parser__
    def capitalize(self):
        return self.__class__(value=self._value.capitalize(), qualifiers=self._qualifiers)

    def center(self,arg1,arg2=' '):
        return self.__class__(value=self._value.center(arg1, arg2), qualifiers=self._qualifiers)

    def count(self, sub, start=0, end=None):
        if end is None:
            end=len(self._value)
        return self._value.count(sub, start, end)

    # decode encode
    def endswith(self, arg):
        return self._value.endswith(arg)

    def find(self, sub, start=0, end=None):
        if end is None:
            end=len(self._value)
        return self._value.find(sub, start, end)

    def expandtabs(self, size=8):
        return self.__class__(value=self._value.expandtabs(size), qualifiers=self._qualifiers)

    # format index
    def isalnum(self):
        return self._value.isalnum()

    def isalpha(self):
        return self._value.isalpha()

    def isdigit(self):
        return self._value.isdigit()

    def islower(self):
        return self._value.islower()

    def isspace(self):
        return self._value.ispace()

    def istitle(self):
        return self._value.istitle()

    def isupper(self):
        return self._value.isupper()

    def join(self, arg):
        return self.__class__(value=self._value.join(arg), qualifiers=self._qualifiers)

    def ljust(self,width,fillChar=' '):
        return self.__class__(value=self._value.ljust(width, fillChar), qualifiers=self._qualifiers)

    def lower(self):
        return self.__class__(value=self._value.lower(), qualifiers=self._qualifiers)

    def lstrip(self):
        return self.__class__(value=self._value.lstrip(), qualifiers=self._qualifiers)

    def partition(self, arg):
        return self._value.partition(arg)

    def replace(self, arg1, arg2):
        return self.__class__(value=self._value.replace(arg1, arg2), qualifiers=self._qualifiers)

    def rfind(self, arg):
        return self._value.rfind(arg)

    # rindex
    def rpartition(self, arg):
        return self._value.rpartition(arg)

    def rjust(self, width, fillChar=' '):
        return self.__class__(value=self._value.rjust(width, fillChar), qualifiers=self._qualifiers)

    def rsplit(self, arg1='', arg2=-1):
        return self._value.rsplit(arg1, arg2)

    def rstrip(self, arg=None):
        return self.__class__(value=self._value.rstrip(arg),qualifiers=self._qualifiers)

    def split(self, arg1=' ', arg2=-1):
        return self._value.split(arg1, arg2)

    def splitlines(self, arg=0):
        return self._value.splitlines(arg)

    def startswith(self, arg):
        return self._value.startswith(arg)

    def strip(self, arg=None):
        if arg is None:
            return ''
        return self.__class__(value=self._value.strip(arg), qualifiers=self._qualifiers)

    def swapcase(self):
        return self.__class__(value=self._value.swapcase(), qualifiers=self._qualifiers)

    def title(self):
        return self.__class__(value=self._value.title(), qualifiers=self._qualifiers)

    # translate
    def upper(self):
        return self.__class__(value=self._value.upper(), qualifiers=self._qualifiers)

    def zfill(self, width):
        return self.__class__(value=self._value.zfill(width), qualifiers=self._qualifiers)


class COneWord(CString):
    '''A single word string - no white space'''
    ERROR_CODES = {201 : {'description' : 'Word contains white space item'}}

    def __init__(self, value=NotImplemented, qualifiers={}, parent=None, name=None,**kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        CString.__init__(self, value=value, qualifiers=qualis, parent=parent, name=name)

    def validity(self, arg):
        # This is space tab and newline
        err= CString.validity(self, arg)
        if arg is not None and isinstance(arg,str):
            s = re.search(self.reWhiteSpacePattern(), arg)
            if s is not None:
                err.append(self.__class__, 201, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return err

    def fix(self, arg):
        if arg is not None and isinstance(arg, str):
            return re.sub(self.reWhiteSpacePattern(), '', arg)
        else:
            return arg


class CRangeSelection(CString):
    ERROR_CODES = {201 : {'description' : 'Range selection contains invalid character'},
                   202 : {'description' : 'Range selection contains bad syntax'}}

    def validity(self, arg):
        err = CErrorReport()
        if arg is None:
            if not self.qualifiers('allowUndefined'):
                err.append(CData, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return err
        arg = self.removeWhiteSpace(arg)
        s = re.search('[^0-9,\,,\-]',arg)
        if s is not None:
            err.append(self.__class__, 201, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        else:
            rList = arg.split(',')
            for r in rList:
                if len(r) < 1:
                    err.append(self.__class__, 202, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                elif r.count('-') > 1:
                    err.append(self.__class__, 202, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                elif r.count('-') == 1 :
                    rr = r.split('-')
                    try:
                        if int(rr[0])>int(rr[1]):
                            err.append(self.__class__, 202, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                    except:
                        err.append(self.__class__, 202, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return err


class CFloat(CBaseData):
    '''A float'''

    CONTENTS = {}
    PYTHONTYPE = float
    QUALIFIERS = {'max' : None, 'min' : None, 'enumerators' : [], 'menuText' : [], 'onlyEnumerators' : False }
    QUALIFIERS_ORDER = ['min','max','onlyEnumerators','enumerators','menuText']
    QUALIFIERS_DEFINITION = {'default' : { 'type' :float },
                             'max'  : { 'description' : 'The inclusive maximum value' },
                             'min'  : {'description' : 'The inclusive minimum value'  },
                             'enumerators'  : {'type' :list,
                                               'description' : 'A Python list of allowed or recommended values - see onlyEnumerators'},
                             'menuText' :  {'type' :list, 'listItemType' : str,
                                            'description' : 'A Python list of strings, matching items in enumerators list, to appear on GUI menu'},
                             'onlyEnumerators' : {'type' : bool,
                                                  'description' : 'If this is true then the enumerators are obligatory - otherwise they are treated as recommended values'}}
    ERROR_CODES = {101 : {'description' : 'below minimum'},
                   102 : {'description' : 'above maximum'},
                   103 : {'description' : 'not one of limited allowed values'}}

    def validity(self,arg):
        validityObj = CErrorReport()
        if arg is None or (isinstance(arg,CFloat) and arg.__dict__['_value'] is None):
            if self.qualifiers('allowUndefined'):
                validityObj.append(self.__class__, 1, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            else:
                validityObj.append(self.__class__, 2, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return validityObj
        if not isinstance(arg, (self.PYTHONTYPE, CFloat)):
            validityObj.append(self.__class__, 5, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
            return validityObj
        if self.qualifiers('min') is not None and arg < self.qualifiers('min'):
            validityObj.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), details='of '+str(self.qualifiers('min')), stack=False)
        if self.qualifiers('max') is not None and arg > self.qualifiers('max'):
            validityObj.append(self.__class__, 102, name=self.objectPath(), details='of '+str(self.qualifiers('max')), label=self.qualifiers('guiLabel'), stack=False)
        if self.qualifiers('onlyEnumerators') and len(self.qualifiers('enumerators')) > 0 and self.qualifiers('enumerators').count(arg) < 1:
            validityObj.append(self.__class__,103,name=self.objectPath(),label=self.qualifiers('guiLabel'),stack=False)
        return validityObj

    def getMenuValue(self):
        if 'menuText' in self._qualifiers and len(self._qualifiers['menuText']) > 0 and self._value in self._qualifiers['enumerators']:
            return self._qualifiers['menuText'][self._qualifiers['enumerators'].index(self._value)]
        else:
            return self._value

    # Applying mathematical functions to self._value in place
    # These methods return a CErrorReport
    def add(self,arg):
        if self._value is not None:
            self.set(self._value + self.getValue(arg))

    def sub(self,arg):
        if self._value is not None:
            self.set(self._value - self.getValue(arg))

    def mul(self,arg):
        if self._value is not None:
            self.set(self._value * self.getValue(arg))

    def div(self,arg):
        if self._value is not None:
            self.set(self._value / self.getValue(arg))

    def abs(self):
        if self._value is not None:
            self.set(abs(self._value))

    def pow(self,arg):
        if self._value is not None:
            self.set(self._value.__pow__(self.getValue(arg)))

    def neg(self,arg):
        if self._value is not None:
            self.set((self._value.__neg__()))

    # Implementation of the usual methods for int
    def __abs__(self):
        return self.__class__(value=self._value.__abs__(),qualifiers=self._qualifiers)

    def __add__(self,arg):
        return self.__class__(value=self._value.__add__(self.getValue(arg)),qualifiers=self._qualifiers)

    # Ho hum. Python float objects do not have a __cmp__ method
    # But we really need CFloat to support this so CInt & CFloat can be handled identically
    def __cmp__(self,arg):
        if self._value == arg:
            return 0
        elif self._value > arg:
            return 1
        else:
            return -1

#These methods are required by Python3 which does not honour __cmp__ anymore : https://portingguide.readthedocs.io/en/latest/comparisons.html
    def __eq__(self, arg):
        if self._value is None and self.getValue(arg) is None:
            return 0
        return self._value == self.getValue(arg)

    def __ne__(self, arg):
        if self._value is None and self.getValue(arg) is not None:
            return 1
        if self._value is not None and self.getValue(arg) is None:
            return 1
        return self._value != self.getValue(arg)

    def __lt__(self, arg):
        if self._value is None or self.getValue(arg) is None:
            raise TypeError("'<' not supported between instances of '"+str(type(self._value))+"' and '"+str(type(self.getValue(arg)))+"'")
        return self._value < self.getValue(arg)

    def __le__(self, arg):
        if self._value is None or self.getValue(arg) is None:
            raise TypeError("'<=' not supported between instances of '"+str(type(self._value))+"' and '"+str(type(self.getValue(arg)))+"'")
        return self._value <= self.getValue(arg)

    def __gt__(self, arg):
        if self._value is None or self.getValue(arg) is None:
            raise TypeError("'>' not supported between instances of '"+str(type(self._value))+"' and '"+str(type(self.getValue(arg)))+"'")
        return self._value > self.getValue(arg)

    def __ge__(self, arg):
        if self._value is None or self.getValue(arg) is None:
            raise TypeError("'>=' not supported between instances of '"+str(type(self._value))+"' and '"+str(type(self.getValue(arg)))+"'")
        return self._value >= self.getValue(arg)

    def __coerce__(self,arg):
        return NotImplemented

    def __div__(self,arg):
        return self.__class__(value=self._value.__div__(self.pythonType()(arg)),qualifiers=self._qualifiers)

    def __divmod__(self,arg):
        return self._value.__divmod__(self.getValue(arg))

    def __float__(self):
        return float(self._value)

    def __round__(self,nsig=0):
        return float(round(self._value,nsig))

    def __floordiv__(self,arg):
        return self.__class__(value=self._value.__floordiv__(self.getValue(arg)),qualifiers=self._qualifiers)

    def  __format__(self,arg):
        return self._value.__format__(arg)

    def __int__(self):
        return self._value.__int__()

    def __invert__(self):
        return self.__class__(value=self._value.__invert__(), qualifiers=self._qualifiers)

    def __long__(self):
        return self._value.__long__()

    def __mod__(self,arg):
        return self.__class__(value=self._value.__mod__(self.getValue(arg)),qualifiers=self._qualifiers)

    def __mul__(self,arg):
        return self.__class__(value=self._value.__mul__(self.getValue(arg)),qualifiers=self._qualifiers)

    def __neg__(self):
        return self.__class__(value=self._value.__neg__(),qualifiers=self._qualifiers)

    def __bool__(self):
        if self._value is None:
            return False
        elif self._value:
            return True
        else:
            return False

    def __nonzero__(self):
        return self.__bool__()

    def __oct__(self):
        return self._value.__oct__()

    def __or__(self,arg):
        return self._value.__or__(self.getValue(arg))

    def __pos__(self):
        return self.__class__(value=self._value.__pos__(), qualifiers=self._qualifiers)

    def __pow__(self, arg):
        return self.__class__(value=self._value.__pow__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __radd__(self,arg):
        return self.__class__(value=self._value.__radd__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rand__(self,arg):
        return self._value.__rand__(arg)

    def __rdiv__(self,arg):
        return self.__class__(value=self._value.__rdiv__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rdivmod__(self,arg):
        return self._value.__rdivmod__(self.getValue(arg))

    def __rfloordiv__(self,arg):
        return self.__class__(value=self._value.__rfloordiv__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rmod__(self,arg):
        return self.__class__(value=self._value.__rmod__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rmul__(self, arg):
        return self.__class__(value=self._value.__rmul__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rpow__(self, arg):
        return self.__class__(value=self._value.__rpow__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rsub__(self, arg):
        return self.__class__(value=self._value.__rsub__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __rtruediv__(self,arg):
        return self._value.__rtruediv__(self.getValue(arg))

    def __str__(self):
        return self._value.__str__()

    def __sub__(self,arg):
        return self.__class__(value=self._value.__sub__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __truediv__(self,arg):
        return self._value.__truediv__(self.getValue(arg))

    def __trunc__(self):
        return NotImplemented


class CUUID(CString):
    def pyType(self):
        #Convert to python type - presently string
        return self.__str__()


def varToUUID(var):
    ret = var.__str__()
    if isinstance(ret, str):
        ret = ret.encode('ascii', 'ignore')
    if not isinstance(ret, str):
        if sys.version_info > (3,0):
            if type(ret) == bytes:
                try:
                    return ret.decode()
                except:
                    pass
        print('CCP4Data.varToUUID', ret, type(ret))
    return ret


class CBoolean(CBaseData):
    '''A Boolean'''
    CONTENTS = {}
    PYTHONTYPE = bool
    QUALIFIERS = {'menuText' : [NotImplemented,NotImplemented]}
    QUALIFIERS_DEFINITION = {'default' : {'type' : bool},
                             'menuText' : {'type' : list,'listItemType' : str,
                                           'description' : 'A list of two string descriptions for true and false'}}
    ERROR_CODES = {101 : {'description' : 'not allowed value'}}

    def set(self, value, checkValidity=True):
        if isinstance(value,str):
            value = value.strip()
            if ['True','true','1'].count(value):
                self.__dict__['_value'] = True
            else:
                self.__dict__['_value'] = False
        else:
            self.__dict__['_value'] = value
        self.updateData()

    def __abs__(self):
        return self._value.__abs__()

    def __add__(self, arg):
        return self._value.__add__(self.getValue(arg))

    def __and__(self, arg):
        return self.__class__(value=self._value.__add__(self.getValue(arg)), qualifiers=self._qualifiers)

    def __cmp__(self, arg):
        return self._value.__cmp__(self.getValue(arg))

    def __corece__(self):
        return NotImplemented

    def __div__(self, arg):
        return self._value.__div__(self.getValue(arg))

    def __divmod__(self, arg):
        return self._value.__divmod__(self.getValue(arg))

    def __float__(self):
        return self._value.__float__()

    def __floordiv__(self, arg):
        return self._value.__floordiv__(self.getValue(arg))

    def __format__(self):
        return NotImplemented

    def __hash__(self):
        return self._value.__hash__()

    def __hex__(self):
        return self._value.__hex__()

    def __index__(self):
        return self._value.__index__()

    def __int__(self):
        return self._value.__int__()

    def __invert__(self):
        return self._value.__invert__()

    def __long__(self):
        return self._value.__long__()

    def __lshift__(self,arg):
        return self._value.__lshift__(arg)

    def __mod__(self,arg):
        return self._value.__mod__(arg)

    def __mul__(self,arg):
        return self._value.__mul__(arg)

    def __neg__(self):
        return self._value.__neg__()

    def __bool__(self):
        if self._value is None:
            return False
        elif self._value:
            return True
        else:
            return False

    def __nonzero__(self):
        return self.__bool__()

    def __oct__(self):
        return self._value.__oct__()

    def __or__(self, arg):
        return self._value.__or__(self.getValue(arg))

    def __pos__(self):
        return self._value.__pos__()

    def __pow__(self, arg):
        return self._value.__pow__(self.getValue(arg))

    def __radd__(self, arg):
        return self._value.__radd__(self.getValue(arg))

    def __rand__(self, arg):
        return self._value.__rand__(self.getValue(arg))

    def __rdiv__(self, arg):
        return self._value.__rdiv__(self.getValue(arg))

    def __rdivmod__(self, arg):
        return self._value.__rdivmod__(self.getValue(arg))

    # reduce reduce_ex
    def __rfloordiv__(self,arg):
        return self._value.__rfloordiv__(self.getValue(arg))

    def __rlshift__(self,arg):
        return self._value.__rlshift__(self.getValue(arg))

    def __rmod__(self,arg):
        return self._value.__rmod__(self.getValue(arg))

    def __rmul__(self,arg):
        return self._value.__rmul__(self.getValue(arg))

    def __ror__(self,arg):
        return self._value.__ror__(self.getValue(arg))

    def __rpow__(self,arg):
        return self._value.__rpow__(self.getValue(arg))

    def __rrshift__(self,arg):
        return self._value.__rrshift__(self.getValue(arg))

    def __rshift__(self,arg):
        return self._value.__rshift__(self.getValue(arg))

    def __rsub__(self,arg):
        return self._value.__rsub__(self.getValue(arg))

    def __rtruediv__(self,arg):
        return self._value.__rtruediv__(self.getValue(arg))

    def __rxor__(self,arg):
        return self._value.__rxor__(self.getValue(arg))

    def __str__(self):
        return self._value.__str__()

    def __sub__(self,arg):
        return self._value.__sub__(self.getValue(arg))

    # subclasshook
    def __truediv__(self,arg):
        return self._value.__truediv__(self.getValue(arg))

    def __trunc__(self):
        return self._value.__trunc__()

    def __xor__(self,arg):
        return self._value.__xor__(self.getValue(arg))
    # conjugate denominator imag numerator real

class CCollection(CData):
    SUBITEM = {'class' : CString}

    def __init__(self, value={}, qualifiers={}, parent=None, name=None, subItem={}, **kw):
        qualis = {}
        qualis.update(qualifiers)
        for key,val in list(kw.items()):
            if key[0:7] != 'subItem':
                qualis[key] = val
        self.__dict__['_subItemObject'] = None
        itemDef = {}
        itemDef.update(self.SUBITEM)
        itemDef.update(subItem)
        for key,val in list(kw.items()):
            if key == 'subItemClass':
                itemDef['class'] = kw['subItemClass']
            elif key == 'subItemClassName':
                itemDef['className'] = kw['subItemClassName']
            elif key[0:8] == 'subItem_':
                if 'qualifiers' not in itemDef:
                    itemDef['qualifiers'] = {}
                itemDef['qualifiers'][key[8:]] = val
        if itemDef.get('className',None) is not None:
            from core import CCP4DataManager
            itemDef['class'] = CCP4DataManager.DATAMANAGER().getClass(itemDef['className'])
        CData.__init__(self, qualifiers=qualis, parent=parent, name=name, build=False)
        self.setSubItem(itemDef)
        # Reparent the _subItemObject - could not make it child of self until after
        self.build()
        if len(value) > 0:
            self.set(value)

    def dataOrder(self):
        return self.__dict__['_dataOrder']

    def setSubItem(self,subItem={}):
        if subItem.get('class',None) is None:
            self.__dict__['_subItemObject'] = None
        else:
            self.__dict__['_subItemObject']=subItem['class'](parent=self,name='subItemObject',qualifiers=subItem.get('qualifiers',{}))

    def subItemClassName(self):
        if self.__dict__.get('_subItemObject',None) is None:
            return None
        else:
            return self.__dict__['_subItemObject'].__class__.__name__

    def subItemClass(self):
        if self.__dict__.get('_subItemObject', None) is None:
            return self.SUBITEM.get('class', None)
        else:
            return self.__dict__['_subItemObject'].__class__

    def subItemQualifiers(self, name=None, default=True, custom=True):
        if self.__dict__.get('_subItemObject', None) is None:
            return {}
        else:
            return self.__dict__['_subItemObject'].qualifiers(name=name, default=default, custom=custom)

    def subItemObject(self):
        return self.__dict__['_subItemObject']

    def qualifierListType(self):
        return self.subItemClass().PYTHONTYPE

    def setQualifiers(self, qualifiers={}, **kw):
        qualis = {}
        qualis.update(qualifiers)
        qualis.update(kw)
        for key, value in list(qualis.items()):
            # Do NOT use isQualifier method here - we dont want to use the CBaseData.QUALIFIERS
            if key in self.QUALIFIERS:
                self.setQualifier(key, value)
            elif key == 'collectionDefault':
                CData.setQualifier('default', value)
            elif self.__dict__['_subItemObject'] is not None and isQualifier(self.subItemClass(), key):
                self.__dict__['_subItemObject'].setQualifier(key,value)

    def setQualifier(self, name, value):
        if name in self.QUALIFIERS:
            CData.setQualifier(self, name, value)
        elif name == 'collectionDefault':
            CData.setQualifier(self, 'default', value)
        elif self.__dict__['_subItemObject'] is not None and isQualifier(self.subItemClass(), name):
            self.__dict__['_subItemObject'].setQualifier(name, value)

    def subItemEtree(self, customOnly=True):
        errors = CErrorReport()
        element = ET.Element('subItem')
        if self.__dict__['_subItemObject'] is None:
            errors.append(self.__class__, 112, name=self.objectPath())
            return element, errors
        classEle = ET.Element('className')
        classEle.text = self.subItemClassName()
        element.append(classEle)
        qualiEle,errs = self.__dict__['_subItemObject'].qualifiersEtree(customOnly=True, tag='qualifiers')
        errors.extend(errs)
        element.append(qualiEle)
        return element,errors

    def itemValidity(self, arg):
        if isinstance(arg, self.subItemClass()):
            if arg.qualifiers() != self.subItemQualifiers():
                return 1
            else:
                return 0
        else:
            if isinstance(arg, dict):
                return 3
            elif isinstance(arg, self.subItemClass().PYTHONTYPE):
                return 2
            else:
                return 4

    def makeItem(self, value= NotImplemented ):
        qualifiers = self.subItemQualifiers(default=False)
        if value is NotImplemented:
            obj = self.subItemClass()(qualifiers=qualifiers, parent=self)
        else:
            obj = self.subItemClass()(value=value, qualifiers=qualifiers, parent=self)
        return obj

    def validItem(self, arg):
        if arg is None:
            return self.makeItem()
        iV = self.itemValidity(arg)
        if iV == 0:
            return arg
        elif iV == 4:
            raise CException(self.__class__, 104, name=self.objectPath())
        elif iV > 0:
            if iV == 1:
                # arg is right CData class but wrong qualifiers
                obj = self.makeItem(arg.get())
            else:
                # arg is right Python type or a dict
                obj = self.makeItem(arg)
            if obj is None:
                raise CException(self.__class__, 106, name=self.objectPath())
            else:
                return obj

    def __getattr__(self, name=None):
        if name is None:
            return self.__dict__['_value']
        elif name == 'SUBITEM':
            return self.__dict__['_subItemObject']
        else:
            propFget = self.Property(name, 'fget')
            if propFget is not None:
                try:
                    return propFget(self)
                except:
                    raise AttributeError("%r object has no attribute %r" % (self.__class__.__name__, name))
            else:
                raise AttributeError("%r object has no attribute %r" % (self.__class__.__name__, name))


class CDict(CCollection):

    PYTHONTYPE = dict
    QUALIFIERS = {'default' : {}}
    CONTENTS = {}
    ERROR_CODES = {101 : {'description' : 'Attempting to access unknown item'},
                   102 : {'description' : 'Unknown error trying to create new item'},
                   103 : {'description' : 'Attempting to add item which is not appropriate class'}}

    def build(self,qualifiers={}):
        self.__dict__['_value'] = {}
        self.__dict__['_dataOrder'] = []

    def __contains__(self,arg):
        return self._value.__contains__(arg)

    def keys(self):
        return list(self._value.keys())

    def has_key(self,arg):
        return arg in self._value

    def items(self):
        return list(self._value.items())

    def __getitem__(self, arg1):
        return self._value.__getitem__(arg1)

    def __delitem__(self, arg):
        if self._dataOrder.count(arg):
            self._dataOrder.remove(arg)
        return self._value.__delitem__(arg)

    def __setitem__(self,arg1,arg2):
        # Should be similar to __setattr__ but without CData.__setattr__ option
        # NB must use buildItem rather than just assign or could end up with any
        # crap as Python dict is not fussy
        if arg1 in self._value:
            return self._value[arg1].set(arg2)
        else:
            # Try making an item called arg1 and setting
            self.buildItem(arg1)
            try:
                return self._value[arg1].set(arg2)
            except CException as e:
                self._value.__delitem__(arg1)
                raise e
            except Exception as e:
                self._value.__delitem__(arg1)
                raise CException(self.__class__, 102, arg1, name=self.objectPath())

    def __getattr__(self, arg):
        if arg in self._value:
            return self._value.__getitem__(arg)
        else:
            return CData.__getattr__(self, arg)

    def __setattr__(self, arg1, arg2):
        # Should be same as __setitem__ but with option to try CData.__setattr__
        if arg1 in self._value:
            return self._value[arg1].set(arg2)
        else:
            try:
                return CData.__setattr__(self, arg1, arg2)
            except:
                # Try making an item called arg1 and setting
                self.buildItem(arg1)
                return self._value[arg1].set(arg2)

    def dataOrder(self):
        return list(self.__dict__['_value'].keys())

    def buildItem(self,key=None):
        self._value[key] = self.subItemClass()(parent=self, qualifiers=self.subItemQualifiers(default=False), name=key)
        self._dataOrder.append(key)

    def clear(self):
        self.__dict__['_dataOrder'] = []
        self.__dict__['_value'] = {}

    def set(self, value={}, **kw):
        myException = CException()
        valueIn = {}
        if len(value) > 0:
            valueIn.update(value)
        if len(kw) > 0:
            valueIn.update(kw)
        for key,val in list(valueIn.items()):
            if key in self._value:
                try:
                    self._value[key].set(val)
                except CException as e:
                    myException.extend(e)
            else:
                self.buildItem(key)
                try:
                    self._value[key].set(val)
                except CException as e:
                    self._value.__delitem__(key)
                    myException.extend(e)
        self.updateData()
        return myException

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        if self.__dict__['_value'] is None or len(self.__dict__['_value']) == 0:
            return False
        else:
            return True

    def setEtree(self, element, checkValidity=True):
        rv = CErrorReport()
        for ele in element:
            key = ele.find('key').text
            value =  ele.find('value').text
            if key not in self._value:
                self.buildItem(key=key)
            try:
                self._value[key].set(value)
            except CException as e:
                rv.extend(e)
            except:
                rv.append(self.__class__, 12, key, name=self.objectPath())
        self.updateData()
        return rv

    def getEtree(self, excludeUnset=True, name=None, useLXML=True):
        if name is None:
            name = self.objectName()
        if name is None or len(name) == 0:
            name = self.className()
        if useLXML:
            element = ET.Element(name)
        else:
            element = ET.Element(name)
        for key in self.dataOrder():
            if useLXML:
                element.append(ET.Element('item'))
            else:
                element.append(ET.Element('item'))
            if useLXML:
                ele = ET.Element('key')
            else:
                ele = ET.Element('key')
            ele.text = key
            element[-1].append(ele)
            element[-1].append(self.__dict__['_value'][key].getEtree(name='value'))
        return element

    def setDefault(self, force=False):
        self.__dict__['_value'] = {}
        default = self.qualifiers('default')
        if default is None or default is NotImplemented or not isinstance(default, dict):
            return
        for key, value in list(default.items()):
            self.__setitem__(key, value)


class CList(CCollection):
    '''A list with all items of one CData sub-class'''

    itemAdded = QtCore.Signal(int)
    itemDeleted = QtCore.Signal(int)

    PYTHONTYPE = list
    CONTENTS = {}
    QUALIFIERS = {'default' : NotImplemented, 'listMinLength' : 0, 'listMaxLength' : NotImplemented,
                  'listCompare' : NotImplemented}
    QUALIFIERS_ORDER = ['listMinLength', 'listMaxLength', 'listCompare']
    QUALIFIERS_DEFINITION = {'default' : {'type' :list},
                             'listMaxLength' : {'type' :int,
                                                'description' : 'Inclusive maximum length of list'},
                             'listMinLength' : {'type' :int,
                                                'description' : 'Inclusive minimum length of list'},
                             'listCompare' : {'type' :int,
                                              'description': 'If has value 1/-1 consecutive items in list must be greater/less than preceeding item. The list item class must have a __cmp__() method.'}}

    ERROR_CODES = {101 : {'description' : 'List shorter than required minimum length'},
                   102 : {'description' : 'List longer than required maximum length'},
                   103 : {'description' : 'Consecutive values in list fail comparison test'},
                   104 : {'description' : 'Attempting to add object of wrong type'},
                   105 : {'description' : 'Attempting to add object of correct type but wrong qualifiers'},
                   106 : {'description' : 'Attempting to add data which does not satisfy the qualifiers for a list item'},
                   107 : {'description' : 'Deleting item will reduce list below minimum length'},
                   108 : {'description' : 'Adding item will extend list beyond maximum length'},
                   109 : {'description' : 'Invalid item class'},
                   110 : {'description' : 'etree (XML) list item of wrong type'},
                   112 : {'description' : 'No list item object set for list'}}

    def coerce(self, arg):
        newArg = []
        # Input is not a list - if it is an item of appropriate type
        # then convert to list length 1
        if not isinstance(arg,(list, CList)):
            testArg = [arg]
        else:
            testArg = arg
        if 1:
            for item in testArg:
                iV = self.itemValidity(item)
                #print 'CList.coerce itemValidity',self.objectName(),iV
                if iV == 0:
                    newArg.append(item)
                else:
                    if iV == 1:
                        itemObj = self.makeItem(item.get())
                    elif iV < 4:
                        itemObj = self.makeItem(item)
                    else:
                        itemObj = None
                    if itemObj is None:
                        return arg
                    newArg.append(itemObj)
        try:
            if self.qualifiers('listMinLength') is not NotImplemented:
                while len(newArg) < self.qualifiers('listMinLength'):
                    newArg.append(self.makeItem())
        except:
            return arg
        return newArg


    def validity(self, arg):
        myException = CErrorReport()
        if not isinstance(arg,(list, CList)):
            myException.append(self.__class__, 5, 'Data not a list', name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        else:
            pythonType = self.subItemObject().pythonType()
            for indx in range(min(len(arg), len(self.__dict__['_value']))):
                if not isinstance(arg[indx],(dict, self.subItemClass(), pythonType)):
                    errorMsg = 'List item is not appropriate type: ' + str(type(arg[indx]))
                    myException.append(self.__class__, 5, errorMsg, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                    break
                err = self.__dict__['_value'][indx].validity(arg[indx])
                myException.extend(err)
        if myException.maxSeverity() > SEVERITY_WARNING:
            return myException
        nUnset = self.removeUnsetItems(ifApply=False)
        listMinLength = self.qualifiers('listMinLength')
        listMaxLength = self.qualifiers('listMaxLength')
        if len(arg) - nUnset < listMinLength:
            myException.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        if listMaxLength is not NotImplemented and listMaxLength is not None and len(arg)>listMaxLength:
            myException.append(self.__class__, 102, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        if self.qualifiers('listCompare') is not NotImplemented:
            cmpValue = self.qualifiers('listCompare')
            for ii in range(1,len(arg)):
                if arg[ii].__cmp__(arg[ii-1]) != cmpValue:
                    myException.append(self.__class__, 103, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
                    break
        return myException

    def removeUnsetListItems(self):
        self.removeUnsetItems(ifApply=True)
        CData.removeUnsetListItems(self)

    def removeUnsetItems(self, ifApply=True):
        # Find all list items that are totally unset (count items with default value as unset)
        # If apply is True then delete the unset items
        count = 0
        for indx in range(self.__dict__['_value'].__len__() - 1, -1, -1):
            if not self.__dict__['_value'][indx].isSet(allowUndefined=True, allowDefault=False, allSet=False):
                count += 1
                if ifApply:
                    del self.__dict__['_value'][indx]
        return count

    def set(self, value=[], validate=False):
        value = self.coerce(value)
        if validate:
            v = self.validity(value)
            if v.maxSeverity() > SEVERITY_WARNING:
                raise v
        self.__dict__['_value'] = []
        for item in value:
            self.__dict__['_value'].append(item)
        self.updateData()


    def unSet(self):
        self.__dict__['_value'] = []
        self.set(self.fix())

    def isSet(self, allowUndefined=False, allowDefault=True, allSet=True):
        if allSet:
            #testing if every item to be set
            for item in self.__dict__['_value']:
                if not item.isSet(allowUndefined=allowUndefined, allowDefault=allowDefault):
                    return False
            return True
        else:
            # testing if any item is set
            for item in self.__dict__['_value']:
                if item.isSet(allowUndefined=allowUndefined, allowDefault=allowDefault, allSet=False):
                    return True
            return False

    def build(self,qualifiers={}):
        self.__dict__['_value'] = []
        fx = self.fix([])
        for item in fx:
            self.__dict__['_value'].append(item)
            s = self.__dict__['_value'][-1].dataChanged
            s.emit()
            @QtCore.Slot()
            def dataChangedHandler():
                self.dataChanged.emit()
            s.connect(dataChangedHandler)

    def fix(self,arg=[]):
        minLength = self.qualifiers('listMinLength')
        rv = []
        if minLength is not NotImplemented and minLength is not None :
            if minLength > 0 and len(arg) < minLength:
                for i in range(len(arg), minLength):
                    rv.append(self.makeItem())
        return rv

    def getEtree(self, excludeUnset=True, name=None, useLXML=True):
        if name is None:
            name = self.objectName()
        if name is None or len(name) == 0:
            name = self.className()
        if useLXML:
            element = ET.Element(name)
        else:
            element = ET.Element(name)
        if len(self.__dict__['_value']) == 0:
            #This seems to be necessary to get empty list written out properly
            element.text = ''
        else:
            for item in self.__dict__['_value']:
                ele = item.getEtree(useLXML=useLXML)
                element.append(ele)
        return element

    def setEtree(self, element, checkValidity=True):
        rv = CErrorReport()
        itemClassName = self.subItemClassName()
        self.unSet()
        ii = 0
        self.blockSignals(True)
        for ele in element:
            ii = ii + 1
            name = str(ele.tag).split('.')[-1]
            if name != itemClassName and name != 'item':
                rv.append(self.__class__, 110, name=self.objectPath())
            else:
                while ii > self.__len__():
                    self.addItem()
                self.__dict__['_value'][ii-1].setEtree(ele)
        self.blockSignals(False)
        self.updateData()
        return rv

    def addItem(self,value= NotImplemented, index = -1):
        obj = self.makeItem(value=value)
        if index < 0:
            self.__dict__['_value'].append(obj)
        else:
            self.__dict__['_value'].insert(index, obj)
        print(obj)
        print(type(obj))
        print(obj.dataChanged)
        obj.dataChanged.connect(self.dataChanged.emit)
        self.itemAdded.emit(index)
        self.updateData()
        return obj

    def setDefault(self, force=False):
        default = self.qualifiers('default')
        if default is NotImplemented:
            return
        self.__dict__['_value'] = []
        if isinstance(default, list):
            for item in default:
                self._value.append(item)

    def getValue(self, arg):
        if isinstance(arg, self.__class__):
            return arg.get()
        else:
            return arg

    def get(self):
        ret = []
        for obj in self.__dict__['_value']:
            ret.append(obj.get())
        return ret

    def get0(self):
        return self.__dict__['_value']

    def dataOrder(self):
        return []

    def __add__(self, arg):
        value = self._value + self.getValue(arg)
        return self.__class__(value=value, qualifiers=self.qualifiers(),
                              subItem={ 'class' : self.subItemClass(), 'qualifiers' : self.subItemQualifiers()})

    # This will only return true if arg is one of the actual objects in the list
    # It will not return true if there is an object with the same value
    def __contains__(self, arg):
        return self._value.__contains__(arg)

    def containsValue(self, arg):
        arg = self.getValue(arg)
        if len(self._value) == 0:
            return False
        if hasattr(self._value[0], '__eq__'):
            for item in self._value:
                if item.__eq__(arg):
                    return True
        elif hasattr(self._value[0], '__cmp__'):
            for item in self._value:
                if item.__cmp__(arg) == 0:
                    return True
        return False

    # __delattr__
    def __delitem__(self, arg, validate=False):
        if self.qualifiers('listMinLength') is not NotImplemented and len(self._value) - 1 < self.qualifiers('listMinLength'):
            raise CException(self.__class__, 107, name=self.objectPath())
        save = self._value[arg]
        rv = self.__dict__['_value'].__delitem__(arg)
        if validate:
            # the delete should only be prevented if it is going below listMinLength
            # the following causes too many issues from other incomplete list items which is
            # a likely condition in the gui
            v = self.validity(self._value)
            if v.maxSeverity() > SEVERITY_WARNING:
                self.__dict__['_value'].insert(arg, save)
                myException = CException()
                myException.extend(v)
                raise myException
        self.itemDeleted.emit(arg)
        self.updateData()
        return rv

    def __delslice__(self, arg1, arg2, validate=False):
        if self.qualifiers('listMinLength') is not NotImplemented and len(self._value) - (arg2 - arg1) < self.qualifiers('listMinLength'):
            raise CException(self.__class__, 107, name=self.objectPath())
        save = self._value[arg1, arg2]
        rv = self.__dict__['_value'].__delslice__(arg1, arg2)
        if validate:
            v = self.validity(self._value)
            if v.maxSeverity() > 1:
                self.__dict__['_value'].insert(arg1, arg2, save)
                raise v
        for indx in range(arg1, arg2):
            self.itemDeleted.emit(indx)
        self.updateData()
        return rv

    def __eq__(self, arg):
        return self._value.__eq__(self.getValue(arg))

    # __format__ __ge__ __getattribute__
    def __getitem__(self, arg):
        return self.__dict__['_value'].__getitem__(arg)

    def __getslice__(self, arg1, arg2):
        return self._value.__getslice__(arg1, arg2)

    # __gt__ __hash__ __iadd__ __imul__
    def __iter__(self):
        return self._value.__iter__()

    # __le__
    def __len__(self):
        return self._value.__len__()

    # __lt__
    def __mul__(self, arg):
        return self.__class__(value=self._value.__mul__(arg), qualifiers=self.qualifiers(),
                              subItem={'class' : self.subItemClass(), 'qualifiers' : self.subItemQualifiers()})

    # __ne__ __reduce__ __reduce_ex__
    def __reversed__(self, arg):
        return self.__class__(value=self._value.__reversed__(arg), qualifiers=self.qualifiers(),
                              subItem={'class' : self.subItemClass(), 'qualifiers' : self.subItemQualifiers()})

    # __rmul__ __setattr__
    def __setitem__(self, indx, arg, validate=False):
        obj = self.validItem(arg)
        save = self._value[indx].get()
        rv = self.__dict__['_value'].__setitem__(indx, obj)
        if validate:
            v = self.validity(self._value)
            if v.maxSeverity() > SEVERITY_WARNING:
                self._value[indx].set(save)
                raise v
        self.updateData()
        return rv

    def __setslice__(self, indx1, indx2, val, validate=False):
        valin = []
        save = []
        for item in val:
            valin.append(self.validItem(item))
        for item in self._value[indx1:indx2]:
            save.append(item.get())
        rv = self.__dict__['_value'].__setslice__(indx1, indx2, valin)
        if validate:
            v = self.validity(self._value)
            if v.maxSeverity() > SEVERITY_WARNING:
                for ii in range(len(save)):
                    idx = indx1 + ii
                    self.__dict__['_value'][idx].set(save[ii])
                raise v
        self.updateData()
        return rv

    # __sizeof__
    def __str__(self):
        text = '['
        if len(self.__dict__['_value']) > 0:
            for item in self.__dict__['_value']:
                text = text + item.__str__() + ','
            text = text[0:-1] + ']'
        else:
            text = text + ']'
        return text
    # __subclasshook__

    def append(self, arg, validate=False):
        if self.qualifiers('listMaxLength') is not NotImplemented and len(self._value)+1>self.qualifiers('listMaxLength'):
            raise CException(self.__class__, 108, name=self.objectPath())
        iV = self.itemValidity(arg)
        if iV > 0:
            obj = self.validItem(arg)
            rv = self.__dict__['_value'].append(obj)
        else:
            rv = self.__dict__['_value'].append(arg)
        if validate:
            v = self.validity(self.__dict__['_value'])
            if v.maxSeverity() > SEVERITY_WARNING:
                del self.__dict__['_value'][-1]
                raise v
        self.__dict__['_value'][-1].dataChanged.connect(self.dataChanged.emit)
        self.updateData()
        return rv

    # This will count number of times the exact same object is in list..
    def count(self, arg):
        return self.__dict__['_value'].count(arg)

    # Count the number of items in list with the same value
    def countValue(self, arg):
        arg = self.getValue(arg)
        if len(self._value) == 0:
            return 0
        rv = 0
        if hasattr(self._value[0], '__eq__'):
            for item in self._value:
                if item.__eq__(arg):
                    rv = rv + 1
        elif hasattr(self._value[0], '__cmp__'):
            for item in self._value:
                if item.__cmp__(arg) == 0:
                    rv = rv + 1
        return rv

    def extend(self, arg, validate=False):
        if not isinstance(arg, (CList, list)):
            raise CException(self.__class__, 106, name=self.objectPath())
        if self.qualifiers('listMaxLength') is not NotImplemented and len(self._value) + len(arg) > self.qualifiers('listMaxLength'):
            raise CException(self.__class__, 108, name=self.objectPath())
        oldLength = len(self._value)
        newList = []
        for item in arg:
            obj = self.validItem(item)
            newList.append(obj)
            obj.dataChanged.connect(self.dataChanged.emit)
        rv = self.__dict__['_value'].extend(newList)
        if validate:
            v = self.validity(self._value)
            if v.maxSeverity() > SEVERITY_WARNING:
                del self.__dict__['_value'][oldLength:-1]
                raise v
        self.updateData()
        self.itemAdded.emit(-1)
        return rv

    def index(self,arg):
        return self.__dict__['_value'].index(arg)

    def indexValue(self,arg):
        ii = -1
        for item in self.__dict__['_value']:
            ii = ii + 1
            if item.__cmp__(arg) == 0:
                return ii
        return -1

    def insert(self, indx, arg, validate=False):
        if self.qualifiers('listMaxLength') is not NotImplemented and len(self._value) + 1 > self.qualifiers('listMaxLength'):
            raise CException(self.__class__, 108, name=self.objectPath())
        obj = self.validItem(arg)
        rv = self.__dict__['_value'].insert(indx, obj)
        if validate:
            v = self.validity(self._value)
            if v.maxSeverity() > SEVERITY_WARNING:
                del self.__dict__['_value'][indx]
                raise v
        obj.dataChanged.connect(self.dataChanged.emit)
        self.itemAdded.emit(indx)
        self.updateData()
        return rv

    def pop(self, arg, validate=False):
        if self.qualifiers('listMinLength') is not NotImplemented and len(self._value) - 1 < self.qualifiers('listMinLength'):
            raise CException(self.__class__, 107, name=self.objectPath())
        #self.disConnectSignal(self._value[arg], 'dataChanged', self.dataChanged.emit)
        save = self._value[arg]
        rv = self._value.pop(arg)
        if validate:
            v = self.validity(self._value)
            if v.maxSeverity() > SEVERITY_WARNING:
                self.__dict__['_value'].insert(arg,save)
                raise v
        self.itemDeleted.emit(arg)
        self.updateData()
        return rv

    def remove(self, arg, validate=False):
        if self.qualifiers('listMinLength') is not NotImplemented and len(self._value) - 1 < self.qualifiers('listMinLength'):
            raise CException(self.__class__, 107, name=self.objectPath())
        indx =  self._value.index(arg)
        rv = self.__dict__['_value'].remove(arg)
        if validate:
            v = self.validity(self._value)
            if v.maxSeverity() > SEVERITY_WARNING:
                self.__dict__['_value'].insert(indx, arg)
                raise v
        self.itemDeleted.emit(indx)
        self.updateData()
        return rv

    def reverse(self, validate=False):
        rv = self.__dict__['_value'].reverse()
        if validate:
            v = self.validity(self._value)
            if v.maxSeverity() > SEVERITY_WARNING:
                self.__dict__['_value'].reverse()
                raise v
        self.updateData()
        return rv

    def sort(self, ccmp=None, key=None, reverse=False):
        save = []
        for item in self._value:
            save.append(item)
        rv = self.__dict__['_value'].sort(ccmp, key, reverse)
        v = self.validity(self._value)
        if v.maxSeverity() > SEVERITY_WARNING:
            self.__dict__['_value'] = []
            for item in save:
                self.__dict__['_value'].append(item)
            raise v
        self.updateData()
        return rv

    def blockSignals(self, mode):
        CObject.blockSignals(self, mode)
        for obj in self._value:
            obj.blockSignals(mode)

    def saveToDb(self):
        ifSave = self.subItemQualifiers('saveToDb', default=False)
        if ifSave:
            objList = []
            for item in self.__dict__['_value']:
                objList.append(item)
            return objList, None, {}
        else:
            return [], None, {}

    # Support Qt QAbstractItemModel
    def data(self, column, role):
#FIXME PYQT - or maybe None? This used to return QVariant.
        return None

    def child(self,row):
        return self.__dict__['_value'].__getitem__(row)

    def childCount(self):
        return self.__dict__['_value'].__len__()

    def columnCount(self):
        return 1

    def abstractModelParent(self):
        return self.__dict__['_abstractModelParent']

    def setAbstractModelParent(self,value):
        self.__dict__['_abstractModelParent'] = value

class COutputFileList(CList):
    '''A list with all items of one CData sub-class'''

    PYTHONTYPE = list
    CONTENTS = {}
    QUALIFIERS = {'default' : NotImplemented, 'listMinLength' : 0, 'listMaxLength' : 10,
        'listCompare' : NotImplemented, 'nameRoot' : NotImplemented}
    QUALIFIERS_ORDER = ['listMinLength', 'listMaxLength', 'listCompare', 'nameRoot']
    QUALIFIERS_DEFINITION = {'default' : {'type' :list},
                             'listMaxLength' : {'type' :int,
                                                'description' : 'Inclusive maximum length of list'},
                             'listMinLength' : {'type' :int,
                                                'description' : 'Inclusive minimum length of list'},
                             'listCompare' : {'type' :int,
                                              'description': 'If has value 1/-1 consecutive items in list must be greater/less than preceeding item. The list item class must have a __cmp__() method.'},
                              'nameRoot' : {'type':str,
                                            'description':'Name hint for the base name of output files'}
                              }

    def getEtree(self, excludeUnset=True, name=None, useLXML=True):
        if name is None:
            name = self.objectName()
        if name is None or len(name) == 0:
            name = self.className()
        if useLXML:
            element = ET.Element(name)
        else:
            element = ET.Element(name)
        if len(self.__dict__['_value']) == 0:
            #This seems to be necessary to get empty list written out properly
            element.text = ''
        else:
            for item in self.__dict__['_value']:
                ele = item.getEtree(useLXML=useLXML)
                element.append(ele)
        print(ET.tostring(element))
        return element


class CRange(CData):
    '''Base class for CIntRange and CFloatRange'''
    CONTENTS_ORDER = ['start', 'end']
    QUALIFIERS = {'compare' : NotImplemented }
    QUALIFIERS_ORDER = ['compare']
    QUALIFIERS_DEFINITION = {'compare' : {'type' : int, 'description' : 'If value is  1/-1 the end value must be greater/less than start.'}}
    ERROR_CODES = {101 : {'description' : 'End of range less than start'},
                   102 : {'description' : 'End of range greater than start'}}

    def fix(self, arg={}):
        compare = self.qualifiers('compare')
        if compare is not NotImplemented:
            try:
                if arg['end'].__cmp__(arg['start']) != compare:
                    ret = {}
                    ret['start'] = arg['end']
                    ret['end'] = arg['start']
                    return ret
            except:
                pass
        return arg

    def getValue(self, arg):
        if isinstance(arg, self.__class__):
            return arg.get()
        else:
            return arg

    def __cmp__(self, arg):
        arg = self.getValue(arg)
        startCmp = self._value['start'].__cmp__(arg['start'])
        if startCmp == 0:
            return self._value['end'].__cmp__(arg['end'])
        else:
            return startCmp


class CIntRange(CRange):
    '''Two integers defining start and end of range'''
    CONTENTS = {'start' : {'class' : CInt},
                'end' : {'class' : CInt}}

    def validity(self, arg):
        # Demo of validity checking when class contents are not orthogonal
        # if qualifier compare is set then relation between start and end is limited
        v = self.itemValidity(arg)
        if v.maxSeverity() > 0:
            return v
        compare = self.qualifiers('compare')
        if compare is NotImplemented:
            return v
        if  not arg.get('end').isSet() or not arg.get('start').isSet():
            return v
        if compare > 0:
            if arg.get('end').__cmp__(arg.get('start')) < 0:
                v.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        else:
            if arg.get('end').__cmp__(arg.get('start')) > 0:
                v.append(self.__class__, 102, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return v

class CFloatRange(CRange):
    '''Two floats defining start and end of range'''
    CONTENTS = {'start' : {'class' : CFloat},
                'end' : {'class' : CFloat}}

    def validity(self, arg):
        # Demo of validity checking when class contents are not orthogonal
        # if qualifier compare is set then relation between start and end is limited
        v = self.itemValidity(arg)
        if v.maxSeverity() > 0:
            return v
        compare = self.qualifiers('compare')
        if compare is NotImplemented:
            return v
        if compare > 0:
            if arg.get('end').__lt__(arg.get('start')):
                v.append(self.__class__, 101, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        else:
            if arg.get('end').__gt__(arg.get('start')):
                v.append(self.__class__, 102, name=self.objectPath(), label=self.qualifiers('guiLabel'), stack=False)
        return v

    def __cmp__(self,other):
        for item in ['start','end']:
            c = self.__dict__['_value'][item].__cmp__(other.__dict__['_value'][item])
            if c!= 0:
                return c
        return 0


class CFollowFromJob(CUUID):
    pass


class CPatchSelection(CData):
    CONTENTS = {'taskName' : { 'class' : CString},
                'patch' : { 'class' : CString}}
    CONTENTS_ORDER = ['taskName', 'patch']

    def __init__(self, value=[], qualifiers={}, parent=None, name=None, build=True, **kw):
        CData.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name, build=build, **kw)
        self.__dict__['patchsForTask'] = []

    def fix(self,arg={}):
        if arg.get('taskName',None) is not None:
            from core import CCP4Modules
            self.__dict__['patchsForTask'] = CCP4Modules.COMFILEPATCHMANAGER().patchForTaskName(taskName=arg['taskName'])
            if arg.get('patch',None) not in self.getPatchList():
                arg['patch'] = None
        return arg

    def getPatchList(self):
        nameList = []
        for name, title in self.__dict__['patchsForTask']:
            nameList.append(name)
        return nameList

    def getPatches(self):
        return self.__dict__['patchsForTask']


class CJobTitle(CString):
    pass


class CJobStatus(CInt):
    pass


class CI2DataType(CString):
    QUALIFIERS = {'enumerators' : ['CPdbDataFile', 'CSeqDataFile', 'CObsDataFile', 'CPhsDataFile', 'CMapCoeffsDataFile',
                                    'CFreeRDataFile', 'CMtzDataFile', 'CDictDataFile', 'CDataFile', 'CInt', 'CFloat', 'CString', 'CRefmacKeywordFile'],
                  'menuText' : []}

    def __init__(self, value=[], qualifiers={}, parent=None, name=None, **kw):
        CString.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name,**kw)
        if len(CI2DataType.QUALIFIERS['menuText']) == 0:
            self.makeMenuText()

    def makeMenuText(self):
        from core import CCP4DataManager
        menu = []
        for name in CI2DataType.QUALIFIERS['enumerators']:
            cls = CCP4DataManager.DATAMANAGER().getClass(className=name)
            if cls is None or cls.QUALIFIERS.get('guiLabel', NotImplemented) is NotImplemented:
                menu.append(name)
            else:
                menu.append(cls.QUALIFIERS['guiLabel'] + ' - ' + name)
        CI2DataType.QUALIFIERS['menuText'] = menu

    def validate(self,arg):
        from core import CCP4DataManager
        cls = CCP4DataManager.DATAMANAGER().getClass(className=arg)
        if cls is None:
            return None
        else:
            return arg

    def isCDataFile(self):
        from core import CCP4DataManager
        from core import CCP4File
        cls = CCP4DataManager.DATAMANAGER().getClass(className=self.__dict__['_value'])
        if cls is not None and issubclass(cls, CCP4File.CDataFile):
            return True
        else:
            return False

    def getMenu(self):
        from core import CCP4DataManager
        menu = []
        for className in self.qualifiers('enumerators'):
            cls = CCP4DataManager.DATAMANAGER().getClass(className=className)
            if cls is not None:
                from core import CCP4File
                if issubclass(cls,CCP4File.CDataFile):
                    label = cls.QUALIFIERS.get('mimeTypeDescription')
                else:
                    label = cls.QUALIFIERS.get('guiLabel')
                menu.append(className,label)
        return menu


#===========================================================================================================
import unittest
def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testCListAppend)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testCListAssorted))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testQObject))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testDict))
    # suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testTable))  # KJS : This is broken
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testCListAppend(unittest.TestCase):
    def setUp(self):
        self.l = CList([2,3,4],listMinLength=3,listCompare=1, subItem= {'class' : CInt, 'qualifiers' : {'default':0,'min':0}})

    def testAppend0(self):
        self.l.append(5)
        self.assertEqual(self.l.get(),[2,3,4,5])

    def testAppend1(self):
        self.l.append(CInt(5))
        self.assertEqual(self.l.get(),[2,3,4,5])

    def testAppend2(self):
        # Switch off the validity checks to allow an 'uninitialised' CInt to append
        self.l.setQualifiers(min=NotImplemented,listCompare= NotImplemented)
        print('testAppend2 qualifiers',self.l.qualifiers(),self.l.subItemQualifiers())
        self.l.append(CInt(default=0))
        self.assertEqual(self.l.get(),[2,3,4,0])

    def testAppend3(self):
        # Should not append value < min
        # Expect CInt 101
        #self.failUnlessRaises(CException,self.l.append,-1)
        try:
            self.l.append(-1)
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in CList.append that should fail item quailfier test')
            self.assertEqual(e[0]['code'],101,'Unexpected exception in CList.append that should fail item qualifier test')
        else:
            self.fail('No exception in CList.append that should fail item quailfier test')

    def testAppend4(self):
        # Should not append value < last item in list
        #self.failUnlessRaises(CException,self.l.append,3)
        try:
            self.l.append(3)
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in CList.append that should fail comparison test')
            self.assertEqual(e[0]['code'],103,'Unexpected exception in CList.append that should fail comparison test')
        else:
            self.fail('No exception in CList.append that should fail comparison test')

class testCListAssorted(unittest.TestCase):
    def setUp(self):
        self.l = CList(subItemClass=CIntRange,listMaxLength=4)

    def testList1(self):
        self.l.set({'start' : 2, 'end':6})
        self.l.append({'start' : 12, 'end':16})
        self.assertEqual(self.l.get(),[{'start' : 2, 'end':6},{'start' : 12, 'end':16}])
        self.assertEqual(self.l[1],{'start' : 12, 'end':16})

    def testList2(self):
        self.l.set({'start' : 2, 'end':6})
        n =  CList({'start' : 12, 'end':16},subItemClass=CIntRange,listMaxLength=4)
        m = self.l + n
        self.assertEqual(m.get(),[{'start' : 2, 'end':6},{'start' : 12, 'end':16}])

    def testList3(self):
        self.l.set({'start' : 2, 'end':6})
        self.l.insert(0,{'start' : 12, 'end':16})
        self.l.reverse()
        self.assertEqual(self.l.get(),[{'start' : 2, 'end':6},{'start' : 12, 'end':16}])

    def testList4(self):
        self.l.set({'start' : 2, 'end':6})
        m =  [  CIntRange(start=4,end=14),CIntRange(start=5,end=15)]
        n = [CIntRange(start=6,end=16),CIntRange(start=7,end=17)]
        self.l.extend(m)
        self.assertEqual(self.l.get(),[{'start' : 2, 'end':6},{'start' : 4, 'end':14},{'start' : 5, 'end':15}])
        try:
            self.l.extend(m)
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in CList.extend')
            self.assertEqual(e[0]['code'],108,'Unexpected exception in CList.extend')
        else:
            self.fail('No exception in CList.extend that should fail')
        self.l.remove({'start' : 4, 'end':14})
        self.assertEqual(self.l[1].end,15,'Fail after CList.remove')

    def testList5(self):
        testXML = '''<CList>
  <CIntRange>
    <start>2</start>
    <end>6</end>
  </CIntRange>
  <CIntRange>
    <start>4</start>
    <end>14</end>
  </CIntRange>
  <CIntRange>
    <start>5</start>
    <end>15</end>
  </CIntRange>
</CList>
'''
        self.l.set([{'start' : 2, 'end':6},{'start' : 4, 'end':14},{'start' : 5, 'end':15}])
        element = self.l.getEtree()
        text = ET.tostring(element,pretty_print=True)
        #print text
        self.assertEqual(text,testXML,'Failed writing XML comparison')
        m = CList(subItemClass=CIntRange,listMaxLength=4)
        m.setEtree(element)
        self.assertEqual(self.l.get(),m.get(),'Failed write/read XML etree')

    def testList6(self):
        self.l.set({'start' : 2, 'end':6})
        j = CList(subItemClass=CIntRange,listMaxLength=4)
        j.append({'start' : 4, 'end':8})
        self.assertEqual(len(j),1,'With two lists - second list wrong length')

class testQObject(unittest.TestCase):
#class testQObject():
    def setUp(self):
        self.bleeped = False
        if QT():
            from PySide2 import QtCore
            self.app = QtCore.QCoreApplication(sys.argv)
            self.master = QtCore.QObject(self.app)

    @QtCore.Slot()
    def bleep(self):
        print('BLEEP!!!')
        self.bleeped = True

    def test1(self):
        if QT():
            from PySide2 import QtCore
            f = CFloat(parent=self.master)
            f.dataChanged.connect(self.bleep)
            f.set(12.0)
            self.assertTrue(self.bleeped,'dataChanged signal not connected')
        else:
            return

class testDict(unittest.TestCase):
    def test1(self):
        d = CDict(subItem={'class' : CIntRange, 'qualifiers' : {'compare' : -1, 'start' : {'min':0}}})
        e = d.set({'foo' : {'start' : 10, 'end' : 5}})
        #print 'd.foo',d.foo
        #print 'error',e
        self.assertEqual(d.foo.start, 10, 'Failed to set CDict')

    def test2(self):
        d = CDict(subItem = {'class' : CIntRange, 'qualifiers' : {'compare' : -1, 'start' : {'min':0}}})
        e = d.set({'foo' : {'start' : -10, 'end' : -20}})
        print('testDict.test2',d,e)
        if len(e) > 0:
            self.assertEqual(e[0]['code'],101,'Setting incorrect dict item does not give correct error code 101')
        else:
            self.fail('Setting incorrect dict item does not give error')

    def test3(self):
        d = CDict(subItem = { 'class' : CIntRange, 'qualifiers' : { 'compare' : -1, 'start' : {'min':0}}})
        e = d.set( { 'foo' : { 'start' : 10, 'end' : 20}} )
        #print 'testDict.test2',e.report()
        self.assertEqual(e[0]['code'],102,'Setting incorrect dict item does not give correct error code 102')

    def test4(self):
        from core.CCP4File import CDataFile
        d = CDict(subItemClass=CDataFile)
        d.PDBIN = { 'project' : 'FOO' , 'baseName' : 'bar.pdb' }
        d.PDBINX = '/foo/bar_x.pdb'
        self.assertEqual(d.dataOrder(),['PDBIN','PDBINX'],'Failed creating Dict using __setattr__')
        self.assertEqual(str(d.PDBINX),'/foo/bar_x.pdb','Failed creating Dict using __setattr__ - 2')
