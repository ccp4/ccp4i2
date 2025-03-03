from __future__ import print_function

from PySide2 import QtCore

"""
     CCP4WorkflowManager.py: CCP4 GUI Project
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
     Liz Potterton July 2013 - create and manage workflows
"""

import os
import re
import glob
from core import CCP4Modules
from core import CCP4Data
from core import CCP4Container
from core import CCP4File
from core import CCP4CustomManager
from .CCP4ErrorHandling import *
from .CCP4TaskManager import TASKMANAGER

class CWorkflowManager(CCP4CustomManager.CCustomManager):

    listChanged = QtCore.Signal()

    ERROR_CODES = {101 : {'description' : 'Error attempting to create directory to save workflow'},
                   102 : {'description' : 'Selected job for workflow did not finish successfully'},
                   103 : {'description' : 'The job parameters file for workflow is not found'},
                   104 : {'description' : 'Error attempting to overwrite directory to save workflow'},
                   105 : {'description' : 'Error attempting to open compressed file for write'},
                   106 : {'description' : 'Error attempting to save workflow to compressed file'},
                   107 : {'description' : 'Error opening compressed file'},
                   108 : {'description' : 'Error extracting from compressed file'},
                   109 : {'description' : 'Compressed file does not have expected content' },
                   110 : {'description' : 'Workflow in compressed file has same name as existing workflow' },
                   111 : {'description' : 'Error attemting to remove a workflow directory' },
                   112 : {'description' : ''},
                   113 : {'description' : ''},}

    insts = None

    def __init__(self,parent=None):
        CCP4CustomManager.CCustomManager.__init__(self, parent, 'workflow')

    def getDefFile(self, name, mustExist=True):
        fileName = os.path.join(self.getDirectory(name=name), 'job_0.def.xml')
        if not mustExist or os.path.exists(fileName):
            return fileName
        else:
            return None

    def createWorkflow(self, projectId=None, jobList=[], name=None, title=None, overwrite=False):
        #print 'CWorkflowMananger.createWorkflow',projectId,jobList,name,overwrite,title
        from dbapi import CCP4DbApi
        workflowDir = self.createDirectory(name, overwrite=overwrite)
        db = CCP4Modules.PROJECTSMANAGER().db()
        projectDir = CCP4Modules.PROJECTSMANAGER().getProjectDirectory(projectId=projectId)
        jobZeroInputParams = CCP4Container.CContainer(parent=self, name='inputData')
        jobZeroOutputParams = CCP4Container.CContainer(parent=self, name='outputData')
        workflowParams = CWorkflowDefinition(name=name)
        workflowParams.jobDef.buildItem('job_0')
        workflowParams.jobDef['job_0'].taskName = 'workflow'
        inputFiles = {}
        nInputFiles = 0
        nOutputFiles = 0
        for jobIndex in range(len(jobList)):
            jobId = jobList[jobIndex]
            jobInfo = db.getJobInfo(jobId=jobId, mode=['status', 'jobnumber', 'taskname'])
            if jobInfo['status'] != 'Finished':
                raise CException(self.__class__, 102, 'Job number ' + str(jobInfo['jobnumber']))
            # Copy the controlParemters of the task def file to the workflow
            self.copyJobDefFile(name, jobIndex, jobInfo['taskname'])
            defFile = os.path.join(projectDir, 'CCP4_JOBS', 'job_' + str(jobInfo['jobnumber']), 'params.xml')
            jobLabel = 'job_' + str(jobIndex + 1)
            #print 'CWorkflowMananger.createWorkflow defFile',defFile
            if not os.path.exists(defFile):
                raise CException(self.__class__, 103, defFile)
            container = CCP4Container.CContainer()
            container.loadDataFromXml(fileName=defFile)
            workflowParams.jobDef.buildItem(jobLabel)
            workflowParams.jobDef[jobLabel].taskName = jobInfo['taskname']
            for key in container.inputData.dataOrder():
                obj0 = container.inputData.__getattr__(key)
                objList, xmlText, keyValues = obj0.saveToDb()
                #print 'createWorkflow input obj0,objList',jobLabel,key,obj0,objList
                for obj in objList:
                    #print 'obj,isSet',key,obj,obj.isSet()
                    if obj.isSet():
                        fileId = CCP4DbApi.UUIDTYPE(obj.dbFileId)
                        if fileId in inputFiles:
                            workflowParams.jobDef[jobLabel].input.append({'toKey': obj.objectName(), 'fromJob' : 'job_0', 'fromKey' : inputFiles[fileId] })
                        else:
                            try:
                                fileinfo = db.getFileInfo(fileId=fileId,mode=['jobid', 'importid', 'jobparamname'])
                            except:
                                fileinfo = {'jobid': None, 'importid' : None, 'jobparamname' : None }
                            #print 'createWorkflow jobIndex key,obj',jobIndex,key,obj,fileinfo
                            if fileinfo['importid'] is not None or fileinfo['jobid'] not in jobList:
                                fromKey = self.uniqueKey(key, jobZeroInputParams)
                                qualifiers = obj.qualifiers(default=False, custom=True)
                                #print 'CWorkflowMananger.createWorkflow qualifiers',obj.objectName(),qualifiers
                                for item in ['contentFlag', 'subType']:
                                    if item in qualifiers:
                                        del qualifiers[item]
                                #print 'WorkFlowManager newObj',key,qualifiers.get('default','NO DEFAULT')
                                newObj = obj.__class__(parent=jobZeroInputParams, name=fromKey)
                                newObj.setQualifiers(qualifiers=qualifiers, validateDefault=False)
                                jobZeroInputParams.addObject(newObj)
                                #print 'CWorkflowMananger.createWorkflow append',str({ 'toKey': obj.objectName(),'fromJob' : 'job_0', 'fromKey' : fromKey })
                                workflowParams.jobDef[jobLabel].input.append({'toKey': obj.objectName(), 'fromJob' : 'job_0', 'fromKey' : fromKey})
                                inputFiles[fileId] = fromKey
                            else:
                                if fileinfo['jobid'] in jobList:
                                    workflowParams.jobDef[jobLabel].input.append({'toKey': obj.objectName(), 'fromJob' : 'job_' + str(jobList.index(fileinfo['jobid']) + 1), 'fromKey' : fileinfo['jobparamname']})
                                else:
                                    print('CWorkflowMananger.createWorkflow Unsure of provenance', obj.objectName())
                # Unset inputData and outputData params and save to the workflow dir
                if len(objList) > 0:
                    obj0.unSet()
            for key in container.outputData.dataOrder():
                obj0 = container.outputData.__getattr__(key)
                objList, xmlText, keyValues = obj0.saveToDb()
                #if len(objList)>0 and obj0.isSet() and obj0.exists():
                if len(objList)>0 and obj0.isSet():
                    workflowParams.jobDef[jobLabel].allOutputFiles.append({'key' : key, 'className' : obj0.__class__.__name__})
                    obj0.unSet()
                    if jobIndex == len(jobList) - 1:
                        toKey = self.uniqueKey(key, jobZeroOutputParams)
                        workflowParams.jobDef['job_0'].output.append({'toKey': toKey, 'fromJob' : 'job_' + str(jobIndex + 1), 'fromKey' : key})
                        qualifiers = obj0.qualifiers(default=False, custom=True)
                        for item in ['contentFlag','subType','sameCrystalAs']:
                            if item in qualifiers:
                                del qualifiers[item]
                        newObj = obj0.__class__(parent=container.outputData, name=toKey, qualifiers=qualifiers)
                        jobZeroOutputParams.addObject(newObj)
                obj0.unSet()
            self.saveJobParams(container=container, name=name, jobLabel=jobLabel)
        #print 'CWorkflowMananger.createWorkflow inputFiles',inputFiles
        workflowParams.header.pluginTitle = title
        #print 'CWorkflowMananger.createWorkflow',self.getCustomFile(name=name,mustExist=False),workflowParams.jobDef.get()
        workflowParams.saveDataToXml(self.getCustomFile(name=name ,mustExist=False), function='WORKFLOW')
        err = self.saveJobZeroDef(name=name, inputDataContainer=jobZeroInputParams, outputDataContainer=jobZeroOutputParams)
        self.listChanged.emit()
        return CErrorReport()

    def uniqueKey(self, keyIn, container):
        import copy
        nn = 0
        key = copy.deepcopy(keyIn)
        while 1:  # KJS : Potential infinite loop if things wrong ?
            try:
                obj = container.__getattr__(key)
            except:
                # No obj with that name so use it
                return key
            else:
                nn += 1
                key = keyIn + str(nn)

    def saveJobParams(self, container=None, name=None, jobLabel=None):
        header = container.header
        #print 'saveJobParams',container.header
        header.setCurrent()
        header.function.set('PARAMS')
        xmlFilePath = os.path.join(self.getDirectory(name), jobLabel + '.params.xml')
        container.saveDataToXml(fileName=xmlFilePath)

    def saveJobZeroDef(self, name=None, title=None, inputDataContainer=None, outputDataContainer=None):
        header = CCP4File.CI2XmlHeader(parent=self)
        header.setCurrent()
        header.function.set('DEF')
        header.pluginName = name
        header.pluginTitle = title
        defContainer = CCP4Container.CContainer(parent=self)
        defContainer.setObjectName(name)
        defContainer.addObject(inputDataContainer)
        defContainer.addObject(CCP4Container.CContainer(name='controlParameters'))
        defContainer.addObject(outputDataContainer)
        xmlFilePath = os.path.join(self.getDirectory(name=name), 'job_0.def.xml')
        defContainer.saveContentsToXml(fileName=xmlFilePath, header=header)

    def copyJobDefFile(self, workflowName, index, taskName, version=None):
        defFile = TASKMANAGER().lookupDefFile(taskName, version=version)
        if defFile is None:
            return
        defContainer = CCP4Container.CContainer(parent=self)
        defContainer.loadContentsFromXml(fileName=defFile)
        defContainer.inputData.clear()
        defContainer.outputData.clear()
        xmlFilePath = os.path.join(self.getDirectory(name=workflowName), 'job_' + str(index+1) + '.def.xml')
        defContainer.saveContentsToXml(fileName=xmlFilePath, header=defContainer.header)


class CWorkflowDataFlow(CCP4Data.CData):
    CONTENTS = {'fromJob' : {'class' : CCP4Data.CString}, 'fromKey' : {'class' : CCP4Data.CString},
                'toKey' : {'class' : CCP4Data.CString}, 'annotation' : {'class' : CCP4Data.CString}}


class CWorkflowDataFlowList(CCP4Data.CList):
    SUBITEM = {'class' : CWorkflowDataFlow}


class CWorkflowFileOut(CCP4Data.CData):
    CONTENTS = {'key' : {'class' : CCP4Data.CString}, 'className' : {'class' : CCP4Data.CString}}
#               'ifOverallOutput' : { 'class' : CCP4Data.CBoolean , 'qualifiers' : { 'default' : False } } }


class CWorkflowJobDefinition(CCP4Data.CData):
    CONTENTS = {'taskName' : { 'class' : CCP4Data.CString}, 'input' : {'class' : CWorkflowDataFlowList},
                'allOutputFiles' : {'class' : CCP4Data.CList, 'subItem' : {'class' : CWorkflowFileOut}},
                'output' : {'class' : CWorkflowDataFlowList}}

    def isSet(self, **args):
        return self.__dict__['_value']['taskName'].isSet()


class CWorkflowJobDefinitionDict(CCP4Data.CDict):
    SUBITEM = { 'class' : CWorkflowJobDefinition }

    def isSet(self,**args):
        #print 'CWorkflowJobDefinitionDict.isSet',self.dataOrder()
        for key in self.dataOrder():
            #print 'CWorkflowJobDefinitionDict.isSet',key,self.__dict__['_value'][key].get('taskName'),self.__dict__['_value'][key].isSet(),
            if self.__dict__['_value'][key].isSet():
                return True
        return False

  
class CWorkflowDefinition(CCP4Container.CContainer):

    def __init__(self, parent=None, name=None, title=None):
        CCP4Container.CContainer.__init__(self, parent=parent, name=name)
        header = self.addHeader()
        header.setCurrent()
        header.function.set('WORKFLOW')
        header.pluginName.set(name)
        header.pluginTitle.set(title)
        self.addObject(CWorkflowJobDefinitionDict(parent=self, name='jobDef'))
        self.addObject(CCP4Data.CString(parent=self, name='title'))

class CWorkflowContainerList(CCP4Data.CDict):
    SUBITEM = {'class' : CCP4Container.CContainer}

    def loadContentsForWorkflow(self, workflowDir):
        # get list of the job.params.xml files - sort and exclude the first (job_0)
        paramFileList = glob.glob(os.path.join(workflowDir, '*.params.xml'))
        numList = []
        for item in paramFileList:
            numList.append(int(os.path.split(item)[1][4:-11]))
        #print 'CContainerList.loadContentsForWorkflow',numList
        numList.sort()
        for i in range(len(numList)):
            numList[i] = 'job_' + str(numList[i])
        # Load container from def.xml file for 'job_0'
        self.buildItem(key = 'job_0')
        self._value['job_0'].loadContentsFromXml(fileName=os.path.join(workflowDir, 'job_0.def.xml'))
        #print 'CWorkflowContainerList.loadContentsForWorkflow numList',numList
        for num in numList:
            self.buildItem(key=num)
            fileName = os.path.join(workflowDir, num + '.params.xml')
            # Kludge unset container name so it is reread from the file header
            # so it is then set to correct taskname for loading contents definition
            self._value[num].setObjectName('')
            rv = self._value[num].loadDataFromXml(fileName, check=False)
        self.__dict__['numList'] = numList
  
    def inputFileList(self,workflowDef=None):
        fileList = []
        for key in self._value['job_0'].inputData.dataOrder():
            dobj = self._value['job_0'].inputData.__getattr__(key)
            info = {'name' : key, 'className' : dobj.className(), 'label' : dobj.qualifiers('guiLabel')}
            for jobKey in workflowDef.jobDef.dataOrder():
                for item in workflowDef.jobDef.__getattr__(jobKey).input:
                    if item.fromJob == 'job_0' and item.fromKey == key:
                        info['taskName'] = str(workflowDef.jobDef.__getattr__(jobKey).taskName)
                        info['paramName'] = str(item.toKey)
                        break
                if 'taskName' in info:
                    break
            #print 'inputFileList key',key,info
            fileList.append(info)
        return fileList
