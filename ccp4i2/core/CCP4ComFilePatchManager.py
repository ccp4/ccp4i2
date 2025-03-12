"""
Copyright (C) 2013 STFC
Liz Potterton July 2013 - create and manage com file patches
"""

import os

from PySide2 import QtCore

from . import CCP4Container
from . import CCP4CustomManager
from . import CCP4Data
from . import CCP4File
from ..googlecode import diff_match_patch_py3
from .CCP4ErrorHandling import CErrorReport, CException


def COMFILEPATCHMANAGER():
    if CComFilePatchManager.insts is None:
        CComFilePatchManager.insts = CComFilePatchManager()
    return CComFilePatchManager.insts


class CComFilePatchManager(CCP4CustomManager.CCustomManager):

    ERROR_CODES = {201 : {'description' : 'The job does not have a com file'},
                   202 : {'description' : 'ERROR while analysing differences in command file'},
                   203 : {'description' : 'ERROR extracting patch from file'},
                   204 : {'description' : 'ERROR applying patch to script'},
                   205 : {'description' : 'ERROR failed to find def file containing default parameters'},
                   206 : {'description' : 'ERROR no job has being selected'}}

    insts = None

    def __init__(self, parent=None):
        CCP4CustomManager.CCustomManager.__init__(self, parent, 'comfilepatch')
        self.lookup = {}
        self.listChanged.connect(self.makeLookup)
        self.makeLookup()

    @QtCore.Slot()
    def makeLookup(self):
        nameList = self.getList()
        self.lookup = {}
        for name in nameList:
            try:
                container = CPatchDefinition(self, name)
                container.loadDataFromXml(fileName=self.getCustomFile(name=name))
                taskName = container.taskNameList[0].__str__()
                if taskName not in self.lookup:
                    self.lookup[taskName] = []
                title = self.getTitle(name, withTaskName=False)
                self.lookup[taskName].append([name, title])
            except:
                print('CComFilePatchManager.makeLookup error loading:', name)
        #print 'CComFilePatchManager.makeLookup', self.lookup

    def patchForTaskName(self,taskName):
        return self.lookup.get(taskName,[])

    def createPatch(self, name, title, taskNameList, projectId, jobId, text1, text2, useControlParams=False, overwrite=False):
        patchDir = self.createDirectory(name, overwrite=overwrite)
        dmp =  diff_match_patch_py3.diff_match_patch()
        try:
            diffs = dmp.diff_main(text1, text2)
            patches = dmp.patch_make(text1, diffs)
            #print 'CComFilePatchManager.createPatch diffs', diffs
            #print 'CComFilePatchManager.createPatch patches', patches
        except:
            raise CException(self.__class__, 202)
        container = CPatchDefinition(parent=self, name=name, title=title)
        container.taskNameList.set(taskNameList)
        container.projectId.set(projectId)
        container.jobId.set(jobId)
        container.text1.set(text1)
        container.text2.set(text2)
        if useControlParams:
            if jobId is None:
                raise CException(self.__class__, 206)
            from core import CCP4TaskManager
            jobInfo = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode=['taskname', 'taskversion'])
            defFile = CCP4TaskManager.TASKMANAGER().lookupDefFile(name=jobInfo['taskname'], version=jobInfo['taskversion'])
            taskContainer = CCP4Container.CContainer()
            taskContainer.loadContentsFromXml(defFile)
            paramsFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=jobId, mode='JOB_INPUT')
            #print 'CComFilePatchManager.createPatch paramsFile',paramsFile
            if not os.path.exists(paramsFile):
                raise CException(self.__class__, 205, paramsFile)
            taskContainer.loadDataFromXml(paramsFile)
            #print 'CComFilePatchManager.createPatch taskContainer.controlParameters',taskContainer.controlParameters
            container.__dict__['_value']['controlParameters'] = taskContainer.controlParameters
        else:
            if 'controlParameters' in container.__dict__['_value']:
                container.__dict__['_value']['controlParameters'].unSet()
        container.saveDataToXml(fileName=self.getCustomFile(name, mustExist=False))
        self.listChanged.emit()
        return CErrorReport()

    def applyPatches(self,name,text):
        container = CPatchDefinition(parent=self, name=name)
        container.loadDataFromXml(fileName=self.getCustomFile(name, mustExist=True))
        dmp =  diff_match_patch_py3.diff_match_patch()
        try:
            diffs = dmp.diff_main(str(container.text1),str(container.text2))
            patches = dmp.patch_make(str(container.text1),diffs)
        except:
            raise CException(self.__class__, 203, name)
        try:
            newText, results = dmp.patch_apply(patches, text)
            #print 'CComFilePatchManager.applyPatch',newText,results
        except:
            raise CException(self.__class__, 204, name)
        return newText, results

    def getComFileText(self, jobId):
        if jobId is None:
            return ''
        try:
            jobDir = CCP4Modules.PROJECTSMANAGER().jobDirectory(jobId=jobId)
            comFilePath = os.path.join(jobDir, 'com.txt')
        except:
            return ''
        if os.path.exists(comFilePath):
            from . import CCP4Utils
            text = CCP4Utils.readFile(comFilePath)
            return text
        else:
            return ''

    def getComFileControlParameters(self, name):
        fileName=self.getCustomFile(name, mustExist=True)
        container = CPatchDefinition(parent=self, name=name)
        container.loadDataFromXml(fileName=fileName, check=False, loadHeader=True)
        # Now need to load the contents of the 'controlParameters' container
        from . import CCP4TaskManager
        defFile = CCP4TaskManager.TASKMANAGER().lookupDefFile(name=container.taskNameList[-1].__str__())
        f = CCP4File.CI2XmlDataFile(fullPath=defFile)
        f.loadFile()
        eleList = f.getBodyEtree().xpath("./container[@id='controlParameters']")
        if len(eleList) > 0:
            container.controlParameters.loadContentsFromEtree(eleList[0])
            container.loadDataFromXml(fileName=fileName, check=False)
        #print 'CComFilePatchManager.getComFileControlParameters name,container', name, container
        return container.controlParameters

    def getTitle(self, name, withTaskName=True):
        from . import CCP4TaskManager
        fileName = self.getCustomFile(name)
        if fileName is None:
            return name
        if not os.path.exists(fileName):
            return name
        container = CPatchDefinition(parent=self, name=name)
        container.loadDataFromXml(fileName=fileName, check=False, loadHeader=True)
        if withTaskName:
            taskName = container.taskNameList[0].__str__()
            #print 'CComFilePatchManager.getTitle',container.header,container.header.pluginTitle
            title = CCP4TaskManager.TASKMANAGER().getTitle(taskName) + ': '
        else:
            title = ''
        if container.header.pluginTitle.isSet():
            title += container.header.pluginTitle.__str__()
        else:
            title += name
        return title

class CPatchDefinition(CCP4Container.CContainer):

    def __init__(self,parent=None, name=None, title=None):
        CCP4Container.CContainer.__init__(self, parent=parent, name=name, title=title)
        taskNameList = CCP4Data.CList(parent=self, name='taskNameList', subItem={'class' : CCP4Data.CString})
        self.addObject(taskNameList)
        self.addObject(CCP4Data.CUUID(parent=self, name='projectId'))
        self.addObject(CCP4Data.CUUID(parent=self, name='jobId'))
        self.addObject(CCP4Data.CString(parent=self, name='text1'))
        self.addObject(CCP4Data.CString(parent=self, name='text2'))
        self.addObject(CCP4Data.CString(parent=self, name='diffs'))
        self.addObject(CCP4Container.CContainer(parent=self, name='controlParameters'))
        header = self.addHeader()
        header.setCurrent()
        header.function.set('COMFILEPATCH')
        header.pluginName.set(name)
        header.pluginTitle.set(title)
