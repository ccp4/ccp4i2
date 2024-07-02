from __future__ import print_function

from PySide6 import QtCore

"""
     CCP4CustomTaskManager.py: CCP4 GUI Project
     Copyright (C) 2013 STFC

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
     Liz Potterton July 2013 - create and manage custom tasks
"""

import os
import re
from core import CCP4Data
from core import CCP4Container
from core import CCP4File
from core import CCP4CustomManager
from core import CCP4DataManager
from core.CCP4ErrorHandling import *


class CCustomTaskManager(CCP4CustomManager.CCustomManager):

    ERROR_CODES = {}
    ERROR_CODES.update(CCP4CustomManager.CCustomManager.ERROR_CODES )
    ERROR_CODES.update({201 : {'description' : 'Data class not recognised'},
                        202 : {'description' : 'Error creating data object'},
                        203 : {'description' : 'Can not find paramter to substitute into command'},
                        204 : {'description' : 'Paramter to substitute into command is unset'},
                        205 : {'description' : ''}})

    insts = None

    def __init__(self, parent=None):
        CCP4CustomManager.CCustomManager.__init__(self, parent, 'task')

    def createCustomTask(self, name=None, title=None, container=None, overwrite=False):
        from core import CCP4XtalData
        err = CErrorReport()
        self.createDirectory(name=name, overwrite=overwrite)
        container.header.pluginName = name
        container.header.pluginTitle = title
        taskFile = self.getCustomFile(name=name, mustExist=False)
        #print 'CCustomTaskManager.createCustomTask', name, title, taskFile, container.paramList
        container.saveDataToXml(fileName=taskFile)
        if title is None:
            import copy
            title = copy.deepcopy(name)
        paramsContainer = CCP4Container.CContainer()
        header = paramsContainer.addHeader()
        header.setCurrent()
        header.function.set('DEF')
        header.pluginName.set(name)
        header.pluginTitle.set(title)
        paramsContainer.addParamsSubContainers()
        mergedMtzs = self.getMergedMtzs(container.paramList)
        for parDef in container.paramList:
            if parDef.function.__str__() in ['input', 'output', 'control parameter'] and \
                    not (parDef.function.__str__() == 'input' and parDef.name.__str__() in mergedMtzs):
                cls = CCP4DataManager.DATAMANAGER().getClass(parDef.dataType.__str__())
                if cls is None:
                    err.append(self.__class__, 201, parDef.dataType.__str__() + ' for ' + parDef.name.__str__())
                else:
                    #try:
                    if 1:
                        qualifiers = {'allowUndefined' : not(parDef.obligatory.__bool__()), 'saveToDb' : parDef.saveDataToDb.__bool__()}
                        if issubclass(cls,CCP4File.CDataFile):
                            qualifiers['fromPreviousJob'] = True
                        if parDef.dataType.__str__() in ['CObsDataFile','CPhsDataFile'] and False in parDef.requiredContentType.get():
                            requiredContentType = parDef.requiredContentType.get()
                            if parDef.dataType.__str__() == 'CObsDataFile':
                                flagList = [CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR, CCP4XtalData.CObsDataFile.CONTENT_FLAG_FPAIR,
                                            CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN, CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN]
                            else:
                                flagList = [CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL, CCP4XtalData.CPhsDataFile.CONTENT_FLAG_PHIFOM]
                            qualifiers['requiredContentFlag'] = []
                            for indx in range(len(flagList)):
                                if requiredContentType[indx]:
                                    qualifiers['requiredContentFlag'].append(flagList[indx])
                        #print 'CCustomTaskManager.createCustomTask',parDef.name.__str__(),qualifiers
                        obj = cls(name=parDef.name.__str__(), qualifiers=qualifiers)
                    #except Exception as e:
                    #  err.append(self.__class__,202,parDef.name.__str__()+' '+str(e))
                    #else:
                        if parDef.function == 'input' :
                            paramsContainer.inputData.addObject(obj, name=parDef.name.__str__())
                        elif parDef.function == 'output' :
                            paramsContainer.outputData.addObject(obj, name=parDef.name.__str__())
                        elif parDef.function == 'control parameter' :
                            paramsContainer.controlParameters.addObject(obj, name=parDef.name.__str__())
        defFile = os.path.join(self.getDirectory(name=name), name + '.def.xml')
        #print 'CCustomTaskManager.createCustomTask defFile', defFile
        err.extend(paramsContainer.saveContentsToXml(fileName=self.getDefFile(name, False), header=header))
        self.listChanged.emit()
        return err

    def getMergedMtzs(self, paramList, function=None):
        mergedMtzs = {}
        for parDef in paramList:
            if parDef.mergeTo.isSet():
                mergeToObj, indx = paramList.getItemByName(parDef.mergeTo.__str__())
                if mergeToObj is not None and (function is None or mergeToObj.function == function):
                    if parDef.mergeTo.__str__() not in mergedMtzs:
                        mergedMtzs[parDef.mergeTo.__str__()] = []
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
        hits = re.findall(tag + '[0-9,a-z,A-Z,\-,_]*', text)
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
                    text = re.sub(tag+paramName, mergedMtzFileNames[paramName], text)
                else:
                    err.append(self.__class__, 203, paramName)
            elif not obj.isSet():
                text = re.sub(tag + paramName, '', text)
                if not obj.qualifier('allowUndefined'):
                    err.append(self.__class__, 204, paramName)
            else:
                text = re.sub(tag + paramName,obj.__str__(), text)
        return text, err


class CCustomComFile(CCP4Data.CData):

    CONTENTS = {'text' : { 'class' : CCP4Data.CString},
                'name' : { 'class' : CCP4Data.CString, 'qualifiers' : {'default' : './com.txt'}}}

    def getTableTextItems(self):
        try:
            text = self.__dict__['_value']['text'].__str__().split('\n')[0]
        except:
            text = ''
        return [ str(self.__dict__['_value']['name']) , text ]


class CCustomComFileList(CCP4Data.CList):

    SUBITEM = {'class' : CCustomComFile}


class CCustomTaskDefinition(CCP4Container.CContainer):

    def __init__(self, parent=None, name=None, title=None):
        CCP4Container.CContainer.__init__(self, parent=parent, name=None, title=None)
        header = self.addHeader()
        header.setCurrent()
        header.function.set('CUSTOMTASK')
        header.pluginName.set(name)
        header.pluginTitle.set(title)
        # Duplicate name/title to simplify handling by gui
        self.addObject(CCP4Data.COneWord(parent=self, name='name'))
        self.addObject(CCP4Data.CString(parent=self, name='title'))
        self.addObject(CCP4Data.CString(parent=self, name='tagCharacter', qualifiers= { 'default' : '#' }))
        self.addObject(CCP4Data.CString(parent=self, name='comLine'))
        self.addObject(CCustomComFileList(parent=self, name='comFileList'))
        #print 'CCustomTaskDefinition.__init__',len(self.comFileList)
        self.comFileList.append({'comFileName' : './com.txt', 'comFileText' : ''})
        self.addObject(CCustomTaskParamList(parent=self, name='paramList'))
        self.paramList.addItem()

    '''
    def setEtree(self,element,checkValidity=True):
        CCP4Container.CContainer.setEtree(self,element,checkValidity=checkValidity)
        print 'CCustomTaskDefinition.setEtree',
        self.name.set(self.header.pluginName.__str__())
        self.title.set(self.header.pluginTitle.__str__())
        
    def getEtree(self,element):
        self.header.pluginName.set(self.name.__str__())
        self.header.pluginTitle.set(self.title.__str__())
        return CCP4Container.CContainer.getEtree(self)
    '''


class CCustomTaskFileFunction(CCP4Data.CString):
    QUALIFIERS = {'enumerators' : ['unknown', 'input' , 'output', 'control parameter', 'log']}


class CCustomTaskParam(CCP4Data.CData):

    CONTENTS = {'name' : {'class' : CCP4Data.COneWord},
                'dataType' : {'class' : CCP4Data.CI2DataType, 'qualifiers' : {'default': 'CPdbDataFile'}},
                'label' : {'class' : CCP4Data.CString},
                'obligatory' : {'class' : CCP4Data.CBoolean, 'qualifiers' : {'default': True}},
                'saveDataToDb' : {'class' : CCP4Data.CBoolean,'qualifiers' : {'default': False}},
                'function' : {'class' : CCustomTaskFileFunction, 'qualifiers' : {'default': 'input'}},
                'mergeTo' : {'class' : CCP4Data.CString},
                'splitColumns' : {'class' : CCP4Data.CString},
                'requiredContentType' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CCP4Data.CBoolean}},  # list length CObsDataFile=4, CPhsDataFile=2, others=0
                'outputFilePath' : {'class' : CCP4Data.CString}}

    def __init__(self, value=None, qualifiers={}, parent=None, name=None, **kw):
        CCP4Data.CData.__init__(self, value=value, qualifiers=qualifiers, parent=parent, name=name, **kw)
        self.__dict__['_value']['dataType'].dataChanged.connect(self.handleDataTypeChange)
        self.handleDataTypeChange()

    @QtCore.Slot()
    def handleDataTypeChange(self):
        # Set requiredContentType list to appropriate length
        if self.__dict__['_value']['dataType'] == 'CObsDataFile':
            self.__dict__['_value']['requiredContentType'].set([True, True, True, True])
        elif self.__dict__['_value']['dataType'] == 'CPhsDataFile':
            self.__dict__['_value']['requiredContentType'].set([True, True])
        else:
            self.__dict__['_value']['requiredContentType'].unSet()
        #print 'CCustomTaskParam.handleDataTypeChange',self.__dict__['_value']['dataType'] ,self.__dict__['_value']['requiredContentType']

    def getTableTextItems(self):
        if self.__dict__['_value']['dataType'].isCDataFile():
            fileFunction = self.__dict__['_value']['function'].__str__()
        else:
            fileFunction = '_'
        if self.__dict__['_value']['obligatory']:
            oblig = 'yes'
        else:
            oblig = 'no'
        if self.__dict__['_value']['saveDataToDb']:
            save = 'yes'
        else:
            save = 'no'
        contentType = ''
        n = 0
        if self.__dict__['_value']['dataType'] in ['CObsDataFile', 'CPhsDataFile']:
            if  self.__dict__['_value']['dataType'] == 'CObsDataFile':
                dataTypeList = ['I+/-', 'F+/-', 'Imean', 'Fmean']
            else:
                dataTypeList = ['HL coeffs', 'Phi/FOM']
            for idx in range(len(dataTypeList)):
                if self.__dict__['_value']['requiredContentType'][idx]:
                    contentType = contentType + dataTypeList[idx] + ' '
                    n = n + 1
            if n == len(dataTypeList):
                contentType = ''
        return [str(self.__dict__['_value']['name']), str(self.__dict__['_value']['label']),
                str(self.__dict__['_value']['dataType']), fileFunction, oblig, save,
                self.__dict__['_value']['mergeTo'].__str__(),
                self.__dict__['_value']['outputFilePath'].__str__(), contentType,
                self.__dict__['_value']['splitColumns'].__str__()]


class CCustomTaskParamList(CCP4Data.CList):

    SUBITEM = {'class' : CCustomTaskParam}

    def getItemByName(self, name=None):
        if name is None:
            return None, None
        for ii in range(len(self.__dict__['_value'])):
            if self.__dict__['_value'][ii].name == name:
                return self.__dict__['_value'][ii],ii
        return None, None

