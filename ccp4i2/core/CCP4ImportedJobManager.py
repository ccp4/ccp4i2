"""
     CCP4ImportedJobManager.py: CCP4 GUI Project




     Liz Potterton Sept 2013 - report job performed outside ccp4i2
"""

import os
import re

from . import CCP4Container
from . import CCP4CustomManager
from . import CCP4Data
from . import CCP4File
from .CCP4ErrorHandling import CErrorReport


def IMPORTEDJOBMANAGER():
    if CImportedJobManager.insts is None:
        CImportedJobManager.insts = CImportedJobManager()
    return CImportedJobManager.insts


class CImportedJobManager(CCP4CustomManager.CCustomManager):

    ERROR_CODES = {}
    ERROR_CODES.update(CCP4CustomManager.CCustomManager.ERROR_CODES)
    ERROR_CODES.update({})
    
    insts = None
    
    def __init__(self, parent=None):
        CCP4CustomManager.CCustomManager.__init__(self, parent, 'importedjobs')

    def getMergedMtzs(self, paramList, function=None):
        mergedMtzs = {}
        for parDef in paramList:
            if parDef.mergeTo.isSet():
                mergeToObj, indx = paramList.getItemByName(parDef.mergeTo.__str__())
                if mergeToObj is not None and (function is None or mergeToObj.function == function):
                    if parDef.mergeTo.__str__() not in mergedMtzs:
                        mergedMtzs[parDef.mergeTo.__str__()] = [ ]
                    if function == 'input':
                        mergedMtzs[parDef.mergeTo.__str__()].append(parDef.name.__str__())
                    else:
                        mergedMtzs[parDef.mergeTo.__str__()].append([parDef.name.__str__(), parDef.splitColumns.get()])
        return mergedMtzs

    def getDefFile(self, name=None, mustExist=True):
        defFile = os.path.join(self.getDirectory(name=name), name + '.def.xml')
        if not mustExist or os.path.exists(defFile):
            return defFile
        else:
            return None

    def extractParams(self, tag, text):
        tagLen = len(tag)
        hits = re.findall(tag + r'[0-9,a-z,A-Z,\-,_]*', text)
        ret = []
        for item in hits:
            ret.append(item[tagLen:])
        return ret

    def substituteParams(self, text='', tag='#', container=None, mergedMtzFileNames={}):
        err = CErrorReport()
        paramNameList = self.extractParams(tag, text)
        for paramName in paramNameList:
            obj = container.find(paramName)
            if obj is None:
                if paramName in mergedMtzFileNames:
                    text = re.sub(tag + paramName, mergedMtzFileNames[paramName], text)
                else:
                    err.append(self.__class__, 203, paramName)
            elif not obj.isSet():
                text = re.sub(tag + paramName, '', text)
                if not obj.qualifier('allowUndefined'):
                    err.append(self.__class__, 204, paramName)
            else:
                text = re.sub(tag + paramName,obj.__str__(), text)
        return text, err


class CImportedJobDefinition(CCP4Container.CContainer):

    def __init__(self, parent=None, name=None, title=None):
        CCP4Container.CContainer.__init__(self, parent=parent, name=None, title=None)
        header = self.addHeader()
        header.setCurrent()
        header.function.set('IMPORTEDJOB')
        header.pluginName.set(name)
        header.pluginTitle.set(title)
        # Duplicate name/title to simplify handling by gui
        self.addObject(CCP4Data.COneWord(parent=self, name='name'))
        self.addObject(CCP4Data.CString(parent=self, name='title'))
        self.addObject(CCP4Data.CBoolean(parent=self, name='ifImportDir'))
        self.addObject(CCP4File.CDataFile(parent=self, name='commandFile', qualifiers={'mustExist' : True}))
        self.addObject(CCP4File.CDataFile(parent=self, name='logFile', qualifiers={'mustExist' : True}))
        self.addObject(CImportedJobDataList(parent=self, name='inputFileDefinitionList'))
        self.addObject(CImportedJobDataList(parent=self, name='outputFileDefinitionList'))


class CImportedJobData(CCP4Data.CData):
    CONTENTS = {'name' : {'class' : CCP4Data.COneWord}, 'dataType' : {'class' : CCP4Data.CI2DataType,
                                                                      'qualifiers' : {'default': 'CPdbDataFile'}},
               'label' : {'class' : CCP4Data.CString},
               'fileName' : {'class' : CCP4File.CDataFile, 'qualifiers' : {'mustExist' : True, 'saveToDb' : True, 'allowUndefined' : False}}}
                #'splitTo' : { 'class' : CCP4Data.CList, 'subItem' : CCP4Data.CInt, 'qualifiers' : { 'listMinLength' : 4, 'listMaxLength': 4, 'default' : [0,0,0,0] } }

    def resetFileNameClass(self, cls):
        if cls is None:
            return
        self.__dict__['_value']['fileName'].deleteLater()
        self.__dict__['_value']['fileName'] = cls(parent=self, name='fileName')


class CImportedJobDataList(CCP4Data.CList):
    SUBITEM = {'class' : CImportedJobData}
    QUALIFIERS = {'listMinLength' : 1}


