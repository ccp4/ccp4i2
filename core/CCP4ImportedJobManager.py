from __future__ import print_function

"""
     CCP4ImportedJobManager.py: CCP4 GUI Project
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
     Liz Potterton Sept 2013 - report job performed outside ccp4i2
"""

import os
import re
from core import CCP4Data
from core import CCP4Container
from core import CCP4File
from core import CCP4CustomManager
from core import CCP4DataManager
from core.CCP4ErrorHandling import *


class CImportedJobManager(CCP4CustomManager.CCustomManager):

    ERROR_CODES = {}
    ERROR_CODES.update(CCP4CustomManager.CCustomManager.ERROR_CODES)
    ERROR_CODES.update({})
    
    insts = None
    
    def __init__(self, parent=None):
        CCP4CustomManager.CCustomManager.__init__(self, parent, 'importedjobs')

    def createImportedJob(self, name=None, title=None, container=None, overwrite=False):
        from core import CCP4XtalData
        err = CErrorReport()
        self.createDirectory(name = name, overwrite=overwrite)
        container.header.pluginName = name
        container.header.pluginTitle = title
        taskFile = self.getCustomFile(name=name, mustExist=False)
        #print 'CCustomTaskManager.createCustomTask',name,title,taskFile,container.paramList
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
                        qualifiers = {'allowUndefined' : not(parDef.obligatory.__bool__()),
                                      'saveToDb' : parDef.saveDataToDb.__bool__()}
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
                        obj = cls(name=parDef.name.__str__(), qualifiers = qualifiers )
                    #except Exception as e:
                    #  err.append(self.__class__, 202, parDef.name.__str__() + ' ' + str(e))
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


