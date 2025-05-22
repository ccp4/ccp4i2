"""
Copyright (C) 2011 University of York
Liz Potterton Mar 2011 - utilities for working with CCP4DbApi
"""

import copy
import functools
import glob
import os
import shutil
import time
import xml.etree.ElementTree as ET

from PySide2 import QtCore

from . import CCP4DbApi
from ..core import CCP4Data
from ..core import CCP4Utils
from ..core.CCP4ErrorHandling import CErrorReport, CException, Severity
from ..core.CCP4Modules import JOBCONTROLLER
from ..core.CCP4Modules import PROJECTSMANAGER
from ..core.CCP4Modules import TASKMANAGER
from ..core.CCP4Modules import WORKFLOWMANAGER


class COpenJob(QtCore.QObject):

    jobStatusUpdated = QtCore.Signal(str)
    jobUpdated = QtCore.Signal(str,str)
    workflowJobStatusUpdated = QtCore.Signal(tuple)
    
    ''' Cache for db job info '''
    ERROR_CODES = {101 : {'description' : 'No projectId set in COpenJob'},
                   102 : {'description' : 'No jobId set in COpenJob'},
                   103 : {'description' : 'Error loading parameters from file'},
                   104 : {'description' : 'Failed to find task definition file - is taskName correct?'},
                   105 : {'description' : 'No container in CopenJob.saveParams'},
                   120 : {'description' : 'Failed to open job - no taskname set'},
                   121 : {'description' : 'Failed running job internally'},}

    MODE = ['taskname', 'taskversion', 'jobnumber', 'projectid',
            'projectname', 'status', 'evaluation', 'jobtitle', 'childjobs']

    def __init__(self, jobId=None, projectId=None):
        QtCore.QObject.__init__(self)
        self.__dict__['_projectId'] = projectId
        self.reset()
        if jobId is not None:
            self.set(jobId)
        PROJECTSMANAGER().db().jobStatusUpdated.connect(self.handleJobStatusUpdated)
        PROJECTSMANAGER().db().jobFinished.connect(self.handleJobStatusUpdated)
        PROJECTSMANAGER().db().jobUpdated.connect(self.handleJobUpdated)
        # Signal from the CDbApi.getRecentlyStartedJobs() that is called repeatedly from PROJECTSMANAGER
        # This should pick up jobs created with 'Running' status by CPluginScript
        PROJECTSMANAGER().db().jobStarted.connect(self.handleJobStarted)

    def reset(self):
        self.__dict__['jobId'] = None
        self.__dict__['_projectName'] = None
        self.__dict__['_projectDir'] = None
        self.__dict__['previous'] = None
        self.__dict__['info'] = {}
        self.__dict__['childjobstaskname'] = []
        self.__dict__['container'] = None
        self.__dict__['isWorkflow'] = False
        self.__dict__['clonedFromJobId'] = None

    def set(self,jobId=None):
        # This assumes that project remains the same
#How can this be? This should be a uuid. Python3 throws a wobbler here. (SJM 13/11/2018)
        #if jobId < 0:
#This is my naive replacement.
        if len(jobId) == 0:
            jobId = self.__dict__['previous']
        else:
            self.__dict__['previous'] = copy.deepcopy(self.jobId)
        self.__dict__['jobId'] = jobId
        if jobId is not None:
            self.__dict__['info'] = PROJECTSMANAGER().db().getJobInfo(jobId=jobId, mode=COpenJob.MODE)
            if self.__dict__['info']['taskname'] is not None:
                self.__dict__['isWorkflow'] = (WORKFLOWMANAGER().getDefFile(self.__dict__['info']['taskname']) is not None)
        else:
            self.__dict__['info'] = {}

    def __setattr__(self, name=None, value=NotImplemented):
        if value is NotImplemented:
            # The 'name' is another instance of COpenJob
            self.__dict__['previous'] = copy.deepcopy(self.jobId)
            self.__dict__['jobId'] = copy.deepcopy(other.jobId)
            self.__dict__['_projectId'] = copy.deepcopy(other.projectId)
            self.__dict__['_projectName'] = copy.deepcopy(other.projectName)
            self.__dict__['isWorkflow'] = copy.deepcopy(other.isWorkflow)
            self.__dict__['childjobstaskname'] = copy.deepcopy(other.childjobstaskname)
            self.__dict__['clonedFromJobId'] = copy.deepcopy(other.clonedFromJobId)
            self.__dict__['info'].update(other.info)
        elif name == 'status':
            if isinstance(value,int):
                value = CCP4DbApi.JOB_STATUS_TEXT[value]
            self.__dict__['info']['status'] = value
        elif name == 'clonedFromJobId':
            self.__dict__['clonedFromJobId'] = value

    def __getattr__(self,key):
        if key in ['jobId','info','container','isWorkflow','childjobstaskname','clonedFromJobId']:
            return self.__dict__[key]
        elif key == 'projectId':
            if self.__dict__['_projectId'] is None:
                self.__dict__['_projectId'] =  self.__dict__['info'].get('projectid',None)
            return self.__dict__['_projectId']
        elif key == 'projectName':
            if self.__dict__['_projectName'] is not None:
                return self.__dict__['_projectName']
            if self.__dict__['_projectId'] is None:
                self.__dict__['_projectId'] =  self.__dict__['info'].get('projectid',None)
            if self.__dict__['_projectId'] is not None:
                self.__dict__['_projectName'] = PROJECTSMANAGER().db().getProjectInfo(projectId= self.__dict__['_projectId'],mode='projectname')
                return self.__dict__['_projectName']
        elif key =='projectDir':
            if self.__dict__['_projectDir'] is not None:
                return self.__dict__['_projectDir']
            if self.__dict__['_projectId'] is None:
                self.__dict__['_projectId'] =  self.__dict__['info'].get('projectid',None)
            if self.__dict__['_projectId'] is not None:
                self.__dict__['_projectDir'] = PROJECTSMANAGER().db().getProjectInfo(projectId=self.__dict__['_projectId'],mode='projectdirectory')
                return self.__dict__['_projectDir']
            else:
                return None
        elif key == 'title':
            text = ''
            if self.jobnumber is not None:
                text = 'Job ' + str(self.jobnumber) + ':  '
            if self.taskname is not None:
                text = text + TASKMANAGER().getTitle(self.taskname)
            return text
        elif key == 'jobDir':
            return  PROJECTSMANAGER().db().jobDirectory(jobId=self.__dict__['jobId'])
        else:
            try:
                return self.__dict__['info'].get(key.lower())
            except:
                print('COpenJob.__getattr__ no' + key)
                return None

    def __str__(self):
        retstring = 'jobId: ' + str(self.__dict__['jobId']) + ' projectId: ' + str(self.__dict__['_projectId']) + \
                    ' projectName: ' + str(self.__dict__['_projectName']) + ' info: ' + str(self.__dict__['info'])
        return retstring

    @QtCore.Slot(dict)
    def handleJobStatusUpdated(self, args):
        #print 'CopenJob.handleJobStatusUpdated',args,self.__dict__['jobId']   # KJS . Check here.
        if args['jobId'] != self.__dict__['jobId']:
            return
        if isinstance(args['status'], int):
            args['status'] = CCP4DbApi.JOB_STATUS_TEXT[args['status']]
        self.__dict__['info']['status'] = args['status']
        self.jobStatusUpdated.emit(args['status'])

    @QtCore.Slot(dict)
    def handleJobUpdated(self,args):
        #print 'COpenJob.handleJobUpdated',args,self.__dict__['jobId']
        if args['jobId'] == self.__dict__['jobId']: 
            if args['key'] == 'evaluation':
                self.__dict__['info']['evaluation'] = CCP4DbApi.JOB_EVALUATION_TEXT[args['value']]
            self.jobUpdated.emit(args['key'], args['value'])
        elif self.__dict__['isWorkflow'] and args['jobId'] in self.__dict__['info']['childjobs'] and args['key'] == 'status':
            #print 'COpenJob.handleJobUpdated workflowJobStatusUpdated'
            childIndx = self.__dict__['info']['childjobs'].index(args['jobId'])
            if childIndx < len(self.__dict__['childjobstaskname']):
                taskName = self.__dict__['childjobstaskname'][childIndx]
            else:
                taskName = None
            self.workflowJobStatusUpdated.emit((self.__dict__['jobId'],args['jobId'],taskName,args['value'])) # KJS. Fixed error (value -> current )

    @QtCore.Slot(dict)
    def handleJobStarted(self,args):
        #print 'COpenJob.handleJobStarted (jobId,projectId,status,parentJobId,taskName)',argList,self.__dict__['isWorkflow']
        if not args['parentJobId'] == self.__dict__['jobId'] or ("childjobs" in self.__dict__['info'] and args['jobId'] in self.__dict__['info']['childjobs']): return
        self.__dict__['info']['childjobs'].append( args['jobId']  )
        self.__dict__['childjobstaskname'].append(args['taskName'])
        if self.__dict__['isWorkflow']:
            #parentJobId,jobId,taskName,status
            self.workflowJobStatusUpdated.emit((args['parentJobId'],args['jobId'],args['taskName'],args['status']))

    def createJob(self,taskName=None,taskVersion=None,cloneJobId=None,contextJobId=None,jobNumber=None,copyInputFiles=False):
        if self.__dict__['_projectId'] is None:
            return CErrorReport(self.__class__,101)
        if taskName is None and cloneJobId is not None:
            pass
        rv = self.createContainer(taskName=taskName,taskVersion=taskVersion)
        if rv.maxSeverity()>Severity.WARNING:
            return rv
        if self.__dict__['jobId'] is None:
            jobId, pName, jNumber = PROJECTSMANAGER().newJob(taskName=taskName,projectId=self.projectId,jobNumber=jobNumber)
            self.set(jobId)
        if cloneJobId is not None:
            self.__dict__['clonedFromJobId'] = cloneJobId
            cloneParamsFile = PROJECTSMANAGER().makeFileName(jobId=cloneJobId, mode='JOB_INPUT')
            try:
                self.__dict__['container'].loadDataFromXml(cloneParamsFile)
            except:
                return CErrorReport(self.__class__, 103, cloneParamsFile)
            if copyInputFiles:
                myDir =  PROJECTSMANAGER().jobDirectory(jobId=self.jobId)
                sourceProjectDir = PROJECTSMANAGER().db().getProjectDirectory(jobId=cloneJobId)
                self.copyInputFiles(self.__dict__['container'].inputData,sourceProjectDir=sourceProjectDir)
        PROJECTSMANAGER().setOutputFileNames(container=self.container, projectId=self.projectId, jobNumber=self.jobnumber, force=True)
        if contextJobId is not None:
            self.setInputByContextJob(contextJobId=contextJobId)
        rv = self.saveParams()
        return rv

    def copyInputFiles(self,container, sourceProjectDir):
        from ..core import CCP4File
        errorReport = CErrorReport()
        for key in container.dataOrder():
            obj0 = container.__getattr__(key)
            try:
                objList, xmlText, keyValues = obj0.saveToDb()
                jobParamName = obj0.objectName()
            except:
                print('ERROR in copyInputFiles for',key)
                objList, xmlText, keyValues = [],None,{}
                jobParamName = ''
            importList = []
            for obj in objList:
                if isinstance(obj,CCP4File.CDataFile) and obj.exists():
                    self.copyInputFile(obj,sourceProjectDir)
                elif isinstance(obj,CCP4Data.CList) and isinstance (obj.subItemObject(),CCP4File.CDataFile):
                    for fileObj in obj:
                        if fileObj.exists():
                            self.copyInputFile(obj,sourceProjectDir)

    def copyInputFile(self,fileObj,sourceProjectDir):
        relPath = os.path.relpath(fileObj.__str__(),sourceProjectDir)
        jobDir = os.path.join('CCP4_JOBS','job_' + str(self.jobNumber) )
        if relPath.startswith(('CCP4_IMPORTED_FILES', os.path.join('CCP4_JOBS','job_' + str(self.jobNumber) ) ) ):
            targetFile = os.path.join(self.projectDir,relPath)
            if not os.path.exists(targetFile):
                try:
                    print('COpenJob.copyInputFile copying',fileObj.__str__(),targetFile)
                    shutil.copyfile(fileObj.__str__(),targetFile)
                except:
                    print('ERROR copying file',fileObj.__str__(),targetFile)
                else:
                    try:
                        dbFileId = str(fileObj.dbFileId)
                        fileInfo = PROJECTSMANAGER().db().getFileInfo(dbFileId,mode='sourceFileName')
                    except:
                        print('COpenJob.copyInputFile failed',fileObj.get())
                    else:
                        print('COpenJob.copyInputFile  dbFileId',dbFileId)
                        print('COpenJob.copyInputFile  sourceFileName',fileInfo)
                        fileObj.setFullPath(targetFile)
                        fileObj.dbFileId.unSet()
                        PROJECTSMANAGER().db().createFile(jobId=self.jobId,projectId=self.projectId,fileObject=fileObj,sourceFileName=fileInfo.get('sourcefilename',None))

    def setInputByContextJob(self,contextJobId=None):
        from ..core import CCP4File
        #Loop over inputData files to pull best file of type from database
        db = PROJECTSMANAGER().db()
        for key in self.__dict__['container'].inputData.dataOrder():
            dobj = self.__dict__['container'].inputData.get(key)
            if isinstance(dobj,CCP4File.CDataFile):
                fileIdList = db.getFileByJobContext(contextJobId=contextJobId,fileType=dobj.qualifiers('mimeTypeName'),
                   subType=dobj.qualifiers('requiredSubType'),contentFlag=dobj.qualifiers('requiredContentFlag'),projectId=self.__dict__['_projectId'])
                #print 'COpenJob.setInputByContextJob',contextJobId,key,dobj.qualifiers('mimeTypeName'),dobj.qualifiers('requiredSubType'),fileIdList
                if len(fileIdList)>0:
                    fileInfo = db.getFileInfo(fileId=fileIdList[0],mode=['jobid','filename','relpath','projectid','annotation','filecontent','filesubtype'])
                    dobj.set({'baseName' : fileInfo['filename'], 'relPath' :  fileInfo['relpath'],
                              'project' : self.__dict__['_projectId'],'annotation' : fileInfo['annotation'],
                              'dbFileId' :fileIdList[0],'contentFlag' :fileInfo['filecontent'],
                              'subType' :fileInfo['filesubtype'] } )
            elif isinstance(dobj,CCP4Data.CList) and isinstance(dobj.subItemObject(),CCP4File.CDataFile):
                subObj = dobj.subItemObject()
                #print  'COpenJob.setInputByContextJob',key, subObj 
                fileIdList = db.getFileByJobContext(contextJobId=contextJobId,fileType=subObj.qualifiers('mimeTypeName'),
                     subType=subObj.qualifiers('requiredSubType'),contentFlag=subObj.qualifiers('requiredContentFlag'),projectId=self.__dict__['_projectId'])
                for idx in range(len(fileIdList)):
                    try:
                        fileInfo = db.getFileInfo(fileId=fileIdList[idx],mode=['jobid','filename','relpath','projectid','annotation','filecontent','filesubtype'])
                    #print 'COpenJob.setInputByContextJob',contextJobId,key,subObj.qualifiers('mimeTypeName'),subObj.qualifiers('requiredSubType'),fileIdList
                    except:
                        print('Error in CDbApi.setInputByContextJob unacceptable fileId',fileIdList[idx])
                    else:
                        while len(dobj)<= idx:
                            dobj.addItem()
                            dobj[len(dobj)-1].set( { 'baseName' : fileInfo['filename'], 'relPath' :  fileInfo['relpath'],
                                                    'project' : self.__dict__['_projectId'],'annotation' : fileInfo['annotation'],
                                                    'dbFileId' :fileIdList[idx],'contentFlag' :fileInfo['filecontent'],
                                                    'subType' :fileInfo['filesubtype'] } )

    def createContainer(self, taskName=None, taskVersion=None):
        from ..core import CCP4Container
        defFile = TASKMANAGER().lookupDefFile(taskName,taskVersion)
        if defFile is None:
            return CErrorReport(self.__class__,104,str(taskName)+' '+str(taskVersion))
        else:
            self.__dict__['info'].update({'taskname' : taskName, 'taskversion' : taskVersion})
        # Set up data container
        self.__dict__['container'] = CCP4Container.CContainer(parent=self,definitionFile=defFile,guiAdmin=True)
        return CErrorReport()
    
    def openJob(self):
        if self.__dict__['jobId'] is None:
            return CErrorReport(self.__class__,102)
        if self.__dict__['info'].get('taskname',None) is None:
            return CErrorReport(self.__class__,120)
        rv = self.createContainer(taskName=self.__dict__['info']['taskname'],taskVersion=self.__dict__['info'].get('taskversion',None))
        if rv.maxSeverity()>Severity.WARNING:
            return rv
        rv = self.loadParams()
        return rv

    def runJob(self,runMode='Now',serverParams=None):
        err =  CErrorReport()
        runMode = runMode.lower()
        if self.__dict__['jobId'] is None:
            return CErrorReport(self.__class__,102)
        if self.__dict__['container'] is None:
            return CErrorReport(self.__class__,105)
        PROJECTSMANAGER().importFiles(jobId=self.jobId,container=self.container)
        rv = self.saveParams()
        if rv.maxSeverity()>Severity.WARNING: return rv
        #Record input files in database
        PROJECTSMANAGER().db().gleanJobFiles(jobId=self.__dict__['jobId'],container=self.__dict__['container'],
                                             projectId=self.__dict__['_projectId'],roleList=[CCP4DbApi.FILE_ROLE_IN])
        self.cleanupJobDir()
        if runMode.startswith('remo'):
            JOBCONTROLLER().createServerParams(self.__dict__['jobId'],serverParams)
        elif TASKMANAGER().isInternalPlugin(self.taskName):
            #print 'COpenJob.runJob isInternalPlugin'
            ok = PROJECTSMANAGER().runInternalTask(jobId=self.jobId,projectId=self.projectId,taskName=self.taskName)
            if not ok: err.append(self.__class__,121,self.taskName)
        else:
            PROJECTSMANAGER().updateJobStatus(jobId=self.jobId,status=CCP4DbApi.JOB_STATUS_QUEUED)
        return err

    def cleanupJobDir(self):
        # Delete a pre-existing report - assume we are restarting job
        try:
            os.remove(PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='REPORT'))
        except:
            pass
        return CErrorReport()

    def deleteJob(self):
        if self.__dict__['jobId'] is None:
            return CErrorReport(self.__class__,102)
        PROJECTSMANAGER().deleteJob(jobId=self.jobId)
        self.reset()
        return CErrorReport()

    def stopJob(self,kill=True,delete=False):
        if self.__dict__['jobId'] is None:
            return CErrorReport(self.__class__,102)    
        if kill:
            err = JOBCONTROLLER().killJobProcess(jobId=self.jobId)
            if delete:
                err.extend(self.deleteJob())  # KJS > Fixed. deleteJob takes no args apparently....
            return err
        else:
            jobDir = PROJECTSMANAGER().jobDirectory(jobId=self.jobId,create=False)
            if jobDir is not None:
                CCP4Utils.saveFile( os.path.join(jobDir,'INTERRUPT'),'Signal to interrupt job')
        return CErrorReport()

    def loadParams(self, fileName=None):
        if fileName is None:
            fileName = PROJECTSMANAGER().makeFileName(jobId = self.__dict__['jobId'],mode='JOB_INPUT')
        try:
            self.__dict__['container'].loadDataFromXml(fileName)
        except:
            return CErrorReport(self.__class__,103,fileName)
        PROJECTSMANAGER().setOutputFileNames(container=self.container,projectId=self.projectId,
                                   jobNumber=self.jobnumber,force=True)
        return CErrorReport()

    def saveParams(self, fileName=None):
        from ..core import CCP4File
        if self.__dict__['jobId'] is None:
            return CErrorReport(self.__class__,102)
        if self.__dict__['container'] is None:
            return CErrorReport(self.__class__,105)
        if fileName is None:
            fileName = PROJECTSMANAGER().makeFileName(jobId=self.jobId,mode='JOB_INPUT')
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        cHeader = self.container.getHeader()
        #print 'COpenJob.saveParams',f.header,cHeader
        if cHeader is not None:
            f.header.set(cHeader)
        f.header.setCurrent()
        f.header.function.set('PARAMS')
        f.header.jobId = self.jobId
        if self.__dict__['info']['projectname'] is not None:
            f.header.projectName.set(self.__dict__['info']['projectname'])
        if self.__dict__['info']['projectid'] is not None:
            f.header.projectId.set(self.__dict__['info']['projectid'])
        if self.__dict__['info']['jobnumber'] is not None:
            f.header.jobNumber.set(self.__dict__['info']['jobnumber'])
        bodyEtree = self.__dict__['container'].getEtree()
        f.saveFile(bodyEtree=bodyEtree)
        return CErrorReport()

    def childOpenJob(self,index):
        if self.__dict__['info']['childjobs'] is not None:
            if index < 0:
                index = len(self.__dict__['info']['childjobs']) + index
            if index >= 0 and index < len(self.__dict__['info']['childjobs']):
                return COpenJob(jobId=self.__dict__['info']['childjobs'][index])
        else:
            return None


# Sorting function
def compareJobNumber(jobA,jobB):
    return cmp(int (jobA.split('_')[-1]),int (jobB.split('_')[-1]))
 
class CProjectDirToDb:

    ERROR_CODES = {101 : { 'description' : 'Project directory does not exist/is not directory'},
                   102 : { 'description' : 'Error writing xml file'}}

    def __init__(self, projectId,projectName, projectDirectory):
        if not os.path.exists(projectDirectory) or not os.path.isdir(projectDirectory):
            raise CException(self.__class__, 101, projectDirectory)
        # !! Need to trap for input error here
        self.projectId = projectId
        self.projectDirectory = projectDirectory
        self.projectName= projectName

    def createDb(self, xmlFile=None):
        if xmlFile is None:
            xmlFile = os.path.join(os.path.split(self.projectDirectory)[0],os.path.split(self.projectDirectory)[1]+'.ccp4i2db.xml')
        self.etree = self.createEtree(self.projectDirectory)
        self.jobsEtree = ET.Element('jobs')
        self.etree.append(self.jobsEtree)
        self.globJobs(os.path.join(self.projectDirectory,'CCP4_JOBS'))
        self.saveEtree(xmlFile)

    def createEtree(self,projectDirectory):
        root = ET.Element('root')
        ele =  ET.Element('projectname')
        ele.text = self.projectName
        root.append(ele)
        ele = ET.Element('projectdirectory')
        ele.text = projectDirectory
        root.append(ele)
        ele = ET.Element('projectcreated')
        ele.text = str(time.time())
        root.append(ele)
        ele = ET.Element('username')
        ele.text = CCP4Utils.getUserId()
        root.append(ele)
        return root

    def saveEtree(self,fileName):
        from ..core import CCP4File
        f = CCP4File.CI2XmlDataFile(fullPath=fileName)
        f.header.setCurrent()
        f.header.function.set('PROJECTDATABASE')
        f.header.projectName = self.projectName
        f.header.projectId = self.projectId
        try:
            f.saveFile(bodyEtree=self.etree)
        except:
            raise CException(self.__class__,102)

    def globJobs(self,dir,parentJobNumber=None,parentJobId=None):
        print('Loading from:',dir,'parentJob:',parentJobNumber,parentJobId)
        from ..core import CCP4Container, CCP4File, CCP4TaskManager
        jobDirList = glob.glob(os.path.join(dir,'job_*'))
        jobDirList.sort(cmp=functools.cmp_to_key(compareJobNumber))
        #print 'CProjectDirToDb.globJobs sorted',jobDirList
        for jobDir in jobDirList:
            header = None
            container = None
            jobId = None
            jobNumber = jobDir.split('_')[-1]
            if parentJobNumber is not None: jobNumber = parentJobNumber + '_' + jobNumber
            paramsFile = self.getParamsFile(jobDir,jobNumber)
            if paramsFile is None:
                print('Error finding params file in',jobDir)
            else:
                try:
                    header = CCP4File.CI2XmlHeader()
                    header.loadFromXml(paramsFile)
                    jobId = int(header.jobId)
                except:
                    print('Error loading container from',paramsFile)
                    header= None
                if header is not None:
                    try:
                        defFile = CCP4TaskManager.TASKMANAGER().lookupDefFile(name=header.pluginName.get(),version=header.pluginVersion.get())
                        #print 'CProjectDirToDb.globJobs defFile',defFile
                        container = CCP4Container.CContainer()
                        container.loadContentsFromXml(defFile,guiAdmin=True)
                        container.loadDataFromXml(paramsFile)
                    except:
                        print('Error loading container from',paramsFile)
                        container = None
                jobEtree = self.getJobTree(header,container,parentJobId)
                self.jobsEtree.append(jobEtree)
                self.globJobs(jobDir,parentJobNumber=jobNumber,parentJobId=jobId)


    def getParamsFile(self, jobDir, jobNumber):
        paramsFile = os.path.join(jobDir, self.projectName + '_' + jobNumber + '.params.xml')
        #print 'paramsFile',paramsFile
        if not os.path.exists(paramsFile):
            paramsFile = os.path.join(jobDir, self.projectName + '_' + jobNumber + '.input_params.xml')
            if not os.path.exists(paramsFile):
                paramsFile = None
        return paramsFile

    def getJobTree(self, header, container, parentJobId=None):
        from ..core import CCP4File
        job = ET.Element('job')
        for item,name  in [['jobId', 'jobid'], ['jobNumber', 'jobnumber'], ['pluginName', 'taskname'],
                           ['creationTime', 'creationtime'], ['creationTime', 'finishtime']]:
            ele = ET.Element(name)
            ele.text = str(header.get(item))
            job.append(ele)
        ele = ET.Element('userAgent')
        ele.text = 'CCP4i2'
        job.append(ele)
        if parentJobId is not None:
            ele = ET.Element('parentJobId')
            ele.text = str(parentJobId)
            job.append(ele)
        if container is not None:
            ele = ET.Element('status')
            ele.text = container.guiAdmin.jobStatus.__str__()
            job.append(ele)
            inputFiles = ET.Element('inputFiles')
            job.append(inputFiles)
            keyList = container.inputData.dataOrder()
            for key in keyList:
                obj = container.inputData.__getattr__(key)
                if isinstance(obj,CCP4File.CDataFile) and not isinstance(obj,CCP4File.CXmlDataFile) \
                                                      and obj.isSet() and obj.dbFileId.isSet():
                    ele = obj.getEtree()
                    ele.tag = str(obj.__class__.__name__)[1:]
                    inputFiles.append( ele )
            outputFiles = ET.Element('outputFiles')
            job.append(outputFiles)
            keyList = container.outputData.dataOrder()
            for key in keyList:
                obj = container.outputData.__getattr__(key)
                if isinstance(obj,CCP4File.CDataFile) and not isinstance(obj,CCP4File.CXmlDataFile) \
                                                      and obj.isSet() and obj.dbFileId.isSet():
                    ele = obj.getEtree()
                    ele.tag = str(obj.__class__.__name__)[1:]
                    outputFiles.append(ele)
            return job


class CMakeProjectDbXml(QtCore.QThread):

    jobLoaded = QtCore.Signal(int)

    ERROR_CODES = {101 : { 'severity' : Severity.WARNING,'description' : 'Error interpreting jobnumber from directory' },
                   102 : { 'severity' : Severity.WARNING,'description' : 'Unknown error attempting to find task info for job' },
                   103 : { 'severity' : Severity.WARNING,'description' : 'Error loading def.xml data file for job' },
                   104 : { 'severity' : Severity.WARNING, 'description' : 'File specified in def file does not exist' },
                   105 : { 'severity' : Severity.WARNING, 'description' : 'File specified in def file does not have dbFileId' },
                   106 : { 'severity' : Severity.WARNING,'description' : 'Can not find file type id for file' },
                   107 : { 'severity' : Severity.WARNING,'description' : 'Error interpreting xdata' },
                   108 : { 'severity' : Severity.WARNING,'description' : 'Error loading job.ccp4db.xml for job number' },
                   109 : { 'severity' : Severity.WARNING,'description' : 'No params.xml or input_params.xml file in job directory' },
                   110 : { 'severity' : Severity.WARNING, 'description' : 'Failed to find import info in job backup file for file' },
                   111 : { 'description' : '' }}

    def __init__(self,db=None,projectDir=None,projectName=None):
        QtCore.QThread.__init__(self)
        self.db=db
        self.projectDir = projectDir
        self.projectName = projectName
        self.projectId = None
        self.diagnostic = True
        self.errReport = CErrorReport()
        self.projectTable = ET.Element('projectTable')
        self.jobTable = ET.Element('jobTable')
        self.fileTable = ET.Element('fileTable')
        self.fileuseTable = ET.Element('fileuseTable')
        self.importfileTable = ET.Element('importfileTable')
        self.exportfileTable = ET.Element('exportfileTable')
        self.xdataTable = ET.Element('xdataTable')
        self.commentTable = ET.Element('commentTable')
        self._fileIdList = []

    def run(self):
        self.loadProject()
        return

    def loadProject(self):
        self.lastJobNumber=None
        topJobList = self.makeJobList(os.path.join(self.projectDir,'CCP4_JOBS'))
        nDone=0
        for jobNum in topJobList:
            jobDir = self.jobDirectory(jobNumber=jobNum)
            #print 'loadProject jobDir',jobDir
            self.loadJob(jobNumber=jobNum,jobDirectory=jobDir)
            self.lastJobNumber=jobNum
            nDone += 1
            self.jobLoaded.emit(nDone)
        ele = self.projectEle(projectId=self.projectId)
        self.projectTable.append(ele)
        return self.errReport
    
    def numberOfJobs(self):
        if self.projectDir is None:
            return 0
        jobDirList = glob.glob(os.path.join(self.projectDir,'CCP4_JOBS','job_*'))
        return len(jobDirList)

    def saveXmlFile(self, xmlFile=None):
        from ..core import CCP4File
        if xmlFile is None:
            xmlFile = self.projectDir+'.ccp4db.xml'
        body = ET.Element('ccp4i2_body')
        for key in ['projectTable','jobTable','fileTable','fileuseTable','importfileTable','exportfileTable','xdataTable','commentTable']:
            body.append(getattr(self,key))
        fileObj = CCP4File.CI2XmlDataFile(fullPath=xmlFile)
        if fileObj.exists(): os.remove(fileObj.__str__())
        fileObj.header.setCurrent()
        fileObj.header.function.set('PROJECTDATABASE')
        if self.projectName is not None:
            fileObj.header.projectName = self.projectName
        else:
            fileObj.header.projectName = self.projectTable[0].get('projectname')
        fileObj.header.projectId = self.projectTable[0].get('projectid')
        fileObj.saveFile(bodyEtree=body)
        return xmlFile
        
    
    def makeJobList(self,directory,jobNumber=None):
        jobDirList = glob.glob(os.path.join(directory,'job_*'))
        jobList = []
        jobNumList = []
        for jobDir in jobDirList:
            try:
                jobNumList.append(int(jobDir.split('_')[-1]))
            except:
                self.errReport.append(self.__class__,101,jobDir)
        jobNumList.sort()
        # Prepend parent job number
        if jobNumber is None:
            jobNumber = ''
        else:
            jobNumber = jobNumber+ '.'
        for item in jobNumList:
            jobList.append(jobNumber+str(item))
        return jobList

    def jobDirectory(self,jobNumber):
        directory =  os.path.join(self.projectDir,'CCP4_JOBS')
        for job in jobNumber.split('.'):
            directory = os.path.join(directory,'job_'+job)
        return directory

    def loadJob(self,jobNumber='',jobDirectory='',parentJobId=None):
        backup = None
        if jobNumber.count('.') == 0 and os.path.exists(os.path.join(jobDirectory,'job.ccp4db.xml')):
            try:
                from ..core import CCP4File
                backupObj = CCP4File.CI2XmlDataFile(os.path.join(jobDirectory,'job.ccp4db.xml'))
                backup = backupObj.getBodyEtree()
            except:
                self.errReport.append(self.__class__,108,jobNumber)
        defFile = os.path.join(jobDirectory,'params.xml')
        #print 'CMakeProjectDbXml.loadJob',jobNumber,jobDirectory,defFile
        status = CCP4DbApi.JOB_STATUS_FINISHED
        if not os.path.exists(defFile):
            defFile = os.path.join(jobDirectory,'input_params.xml')
            jobFiles = glob.glob(os.path.join(jobDirectory,'*'))
            if len(jobFiles) <= 1:
                status = CCP4DbApi.JOB_STATUS_PENDING
            else:
                status = CCP4DbApi.JOB_STATUS_FAILED
            if backup is None:
                self.errReport.append(self.__class__,109,jobDirectory)
                return
            else:
                container = None
                header = backupObj.header
        else:
            header, container = self.loadContainer(defFile)
        print('CMakeProjectDbXml.loadJob',header,repr(container),container)
        ctime = os.path.getctime(jobDirectory)
        ele,jobId = self.jobEle(jobNumber=jobNumber,status=status,header=header,parentJobId=parentJobId,ctime=ctime)
        if backup is not None:
            try:
                backupEle = backup.xpath('./jobTable/job')[0]
                #print 'CMakeProjectDbXml.loadJob backupEle',backupEle
                for key,value in list(backupEle.items()):
                    if key == 'evaluation':
                        # Beware that an unset evaluation is usually recorded as '0' in the job backup xml
                        # but '0' is not an allowed value.  This arises from CDbApi.getJobInfo return 'Unknown' for unset evaluation
                        if value != '0':
                            ele.set(key,value)
            except:
                pass
        self.jobTable.append(ele)
        if container is not None:
            self.loadJobParams(jobId=jobId,jobDirectory=jobDirectory,container=container,backup=backup)
        self.loadFromBackup(jobId=jobId,backup=backup)
        subJobList = self.makeJobList(directory=jobDirectory,jobNumber=jobNumber)
        for subJob in subJobList:
            jobDir = self.jobDirectory(jobNumber=subJob)
            self.loadJob(jobNumber=subJob,jobDirectory=jobDir,parentJobId=jobId)

    def loadJobParams(self,jobId=None,jobDirectory=None,container=None,backup=None):
        #This closely parallels CDbApi.gleanJobFiles()
        from ..core import CCP4File
        for subcontainer,file_role in [[container.inputData,CCP4DbApi.FILE_ROLE_IN],[container.outputData,CCP4DbApi.FILE_ROLE_OUT]]:
            keyList = subcontainer.dataOrder()
            print('loadJobParams keyList',repr(subcontainer),keyList)
            for key in keyList:
                obj = subcontainer.__getattr__(key)
                if isinstance(obj,CCP4File.CDataFile) and not isinstance(obj,CCP4File.CXmlDataFile):
                    if obj.isSet():
                        path = os.path.join(self.projectDir,obj.relPath.__str__(),obj.baseName.__str__())
                        print('loadJobParams', obj.objectName(),path, os.path.exists(path))
                        if not os.path.exists(path) :
                            self.errReport.append(self.__class__,104,str(obj))
                        else:
                            if not obj.dbFileId.isSet():
                                fileId = None
                            else:
                                fileId = CCP4DbApi.UUIDTYPE(obj.dbFileId)
                            if file_role==CCP4DbApi.FILE_ROLE_IN:
                                # Its a pre-existing input file
                                if fileId is None:
                                    self.errReport.append(self.__class__,105,str(obj))
                                elif not self._fileIdList.count(fileId):
                                    # File not already registered - implies its an imported file
                                    # Do we have extra info from backup?
                                    sourcefilename = None
                                    creationtime = None
                                    backupEle = None
                                    importid= None
                                    if backup is not None:
                                        try:
                                            backupEle = backup.xpath("./importfileTable/importfile[@fileid='"+fileId+"']")[0]
                                            #print 'loadJobParams backupEle',backupEle
                                            sourcefilename = backupEle.get('sourcefilename',None)
                                            creationtime = backupEle.get('creationtime',None)
                                            importid = backupEle.get('importid',None)
                                        except:
                                            pass
                                            #self.errReport.append(self.__class__,110,str(obj))
                                            #print 'loadJobParams importfile',backup,backupEle,sourcefilename,creationtime
                                    if importid is not None and fileId is not None:
                                        ele = self.importEle(importid=importid,fileid=fileId,sourcefilename=sourcefilename,creationtime=creationtime)
                                        self.importfileTable.append(ele)
                                    ele = self.fileEle(jobId=jobId,fileObject=obj)
                                    if ele is not None:
                                        self.fileTable.append(ele)
                                        self._fileIdList.append(fileId)
                                else:
                                    # Make a FileUse record
                                    ele = self.fileuseEle(jobId=jobId,fileId=fileId)
                                    self.fileuseTable.append(ele)
                            else:
                                # Its an output file but may be replicate of output file in parent job
                                if  fileId is not None and self._fileIdList.count(fileId):
                                    # Make a FileUse record
                                    ele = self.fileuseEle(jobId=jobId,fileId=fileId,roleId=CCP4DbApi.FILE_ROLE_OUT)
                                    self.fileuseTable.append(ele)
                                else:
                                    ele = self.fileEle(jobId=jobId,fileObject=obj)
                                    if ele is not None:
                                        self.fileTable.append(ele)
                                        self._fileIdList.append(obj.dbFileId.__str__())

    def loadFromBackup(self,jobId=None,backup=None):
        if backup is None: return
        exportEleList = backup.xpath('./exportfileTable/exportfile')
        for exEle in exportEleList:
            fileid = exEle.get('fileid',None)
            if fileid is not None:
                ele = self.exportEle(exportid=exEle.get('exportid', None), fileId=fileid,
                                     exportfilename=exEle.get('exportfilename', None), creationtime=exEle.get('creationtime', None))
                self.exportfileTable.append(ele)
        commentEleList = backup.xpath('./commentTable/comment')
        for comEle in commentEleList:
            ele = self.commentEle(jobId=jobId,commentid=comEle.get('commentid',None),username=comEle.get('username',None),
                                  timeofcomment=comEle.get('timeofcomment',None),comment=comEle.get('comment',None))
            self.commentTable.append(ele)

    def commentEle(self,jobId=None,commentid=None,username=None,timeofcomment=None,comment=None):
        ele = ET.Element('comment')
        ele.set('commentid',commentid)
        if jobId is not None:
            ele.set('jobid',str(jobId))
        if username is not None:
            ele.set('username',username)
        if timeofcomment is not None:
            ele.set('timeofcomment',str(timeofcomment))
        if comment is not None:
            ele.set('comment',comment)
        return ele

    def exportEle(self,exportid=None,fileId=None,exportfilename=None,creationtime=None):
        ele = ET.Element('exportfile')
        ele.set('exportid',str(exportid))
        if fileId is not None: ele.set('fileid',str(fileId))
        if exportfilename is not None: ele.set('exportfilename',exportfilename)
        if creationtime is not None: ele.set('creationtime',creationtime)
        return ele

    def xdataEle(self,xdataId=None,jobId=None,dataObject=None):
        try:
            dataEle =dataObject.getEtree()
            dataClass = dataObject.__class__.__name__
        except:
            self.errReport.append(self.__class__,107,str(dataObject))
            return None,None
        self._lastXdataId = self._lastXdataId + 1
        ele = ET.Element('xdata')
        ele.set('xdataid',str(xdataId))
        ele.set('xdataclass',str(dataClass))
        ele.set('jobid',str(jobId))
        dataEle.tag = 'xdataxml'
        ele.append(dataEle)
        #print 'xdataEle',ET.tostring(ele)
        return ele

    def importEle(self,importid=None,fileid=None,sourcefilename=None,creationtime=None):
        # !!! Needs more info
        ele = ET.Element('importfile')
        ele.set('importid',str(importid))
        ele.set('fileid',str(fileid))
        if sourcefilename is not None: ele.set('sourcefilename',sourcefilename)
        if creationtime is not None: ele.set('creationtime',creationtime)
        return ele

    def fileuseEle(self,jobId=None,fileId=None,roleId=CCP4DbApi.FILE_ROLE_IN):
        ele = ET.Element('fileuse')
        ele.set('fileid',str(fileId))
        ele.set('jobid',str(jobId))
        ele.set('roleid',str(roleId))
        return ele

    def fileEle(self,fileId=None,jobId=None,fileObject=None):
        ele = ET.Element('file')
        if fileId is not None:
            ele.set('fileid',str(fileId))
        elif fileObject.dbFileId.isSet():
            ele.set('fileid',fileObject.dbFileId.__str__())
        else:
            return None
        ele.set('jobid',str(jobId))
        ele.set('filename',fileObject.baseName.__str__())
        # relPath now deduced from job
        #ele.set('relpath',fileObject.relPath.__str__())
        if fileObject.annotation.isSet() and len(fileObject.annotation)>0:
            ele.set('annotation',fileObject.annotation.__str__())
        mimeType = fileObject.qualifiers('mimeTypeName')
        if mimeType is not None:
            fileTypeId = 0
            for fileTypeDef in CCP4DbApi.FILETYPELIST:
                if fileTypeDef[1] == mimeType:
                    fileTypeId = fileTypeDef[0]
                    break
            ele.set('filetypeid',str(fileTypeId))
        else:
            self.errReport.append(self.__class__,106,str(fileObject))
        return ele

    def jobEle(self,jobNumber='',status=None,header=None,parentJobId=None,ctime=0):
        ele = ET.Element('job')    
        jobId = header.jobId.__str__()
        ele.set('jobid',jobId)
        ele.set('jobnumber',header.jobNumber.__str__())
        # Try to deal with  sub-job def file headers without projectid
        # by saving self.projectId when first reasonable value found
        projectId = header.projectId.__str__()
        print('CMakeProjectDbXml.jobEle header',header)
        if len(projectId)==0:
            projectId = self.projectId
        elif self.projectId is None:
            self.projectId = projectId
        ele.set('projectid',projectId)
        if header.creationTime.isSet():
            ele.set('creationtime',header.creationTime.__int__().__str__())
        else:
            ele.set('creationtime',str(ctime))
        ele.set('status',str(status))
        if parentJobId is not None: ele.set('parentjobid',str(parentJobId))
        if header.pluginName.isSet():
            ele.set('taskname',header.pluginName.__str__())
        else:
            ele.set('taskname','Unknown')
        if header.pluginTitle.isSet(): ele.set('jobtitle',header.pluginTitle.__str__())
        ele.set('useragent','1')
        if ctime != 0: ele.set('finishtime',str(ctime))
        return ele,jobId

    def projectEle(self,projectId=None):
        ele = ET.Element('project')
        ele.set('projectid',projectId)
        ele.set('projectname',self.projectName)
        ele.set('projectdirectory',self.projectDir)
        ele.set('lastjobnumber',str(self.lastJobNumber))
        userId = CCP4Utils.getUserId()
        if userId is not None: ele.set('username',userId)
        return ele

    def loadContainer(self,fileName):
        from ..core import CCP4Container, CCP4File, CCP4TaskManager
        header = CCP4File.CI2XmlHeader()
        header.loadFromXml(fileName)
        #print 'CMakeProjectDbXml.loadContainer',header
        if self.projectName is None:
            self.projectName=header.project.__str__()
        try:
            defFile = CCP4TaskManager.TASKMANAGER().lookupDefFile(name=header.pluginName.__str__(),version=header.pluginVersion.__str__())
        except CException as e:
            self.errReport.extend(e)
            return header,None
        except Exception as e:
            self.errReport.append(self.__class__,102,'Job number: ' + header.jobNumber.__str__() + '  ' + str(e))   # KJS : Fixed a couple of typos.
            return header,None
        container = CCP4Container.CContainer()
        #print 'CMakeProjectDbXml.loadContainer',defFile
        err = container.loadContentsFromXml(defFile,guiAdmin=True)
        self.errReport.extend(err)
        if err.maxSeverity()>Severity.WARNING: return None
        try:
            container.loadDataFromXml(fileName=fileName,guiAdmin=True)
        except CException as e:
            self.errReport.extend(e)
            return header,None
        except Exception as e:
            self.errReport.append(self.__class__,103,'Job number: ' + header.jobNumber.__str__() + '  ' + str(e))
            return header,None
        return header,container

class CJobDbBackup:
    # Backup file in job directory

    def __init__(self,jobId=None,jobNumber=None,taskName=None,jobDirectory=None,projectName=None,fileName=None,projectId=None):
        from ..core import CCP4File
        self.jobId = jobId
        self._diagnostic = False
        if fileName is not None:
            self.xmlFile = CCP4File.CI2XmlDataFile(fullPath=fileName)
        else:
            self.xmlFile = CCP4File.CI2XmlDataFile(fullPath=os.path.join(jobDirectory,'job.ccp4db.xml'))
        if self.xmlFile.exists():
            root = self.xmlFile.getEtreeRoot()
            #print 'CJobDbBackup root',root.tag
            self.body = copy.deepcopy(root.xpath('./ccp4i2_body')[0])
        else:
            jobInfo = {}
            if jobNumber is None or taskName is None:
                jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode=['jobnumber','taskname'])        
            if projectName is None: projectName = PROJECTSMANAGER().db().getProjectInfo(projectId=projectId,mode='projectname')
            self.xmlFile.header.setCurrent()
            self.xmlFile.header.function.set('PROJECTDATABASE')
            self.xmlFile.header.jobId = jobId
            self.xmlFile.header.jobNumber = jobInfo.get('jobnumber',jobNumber)
            self.xmlFile.header.pluginName = jobInfo.get('taskname',taskName)
            self.xmlFile.header.projectName = projectName
            self.xmlFile.header.projectId = projectId
            self.body = self.buildTree(jobId=jobId,jobNumber=jobNumber)

    def save(self,updateAll=False):
        if updateAll: self.scrapeFromDb(db=PROJECTSMANAGER().db())
        self.xmlFile.saveFile(bodyEtree=self.body)

    def buildTree(self,jobId=None,jobNumber=None):
        top = ET.Element('ccp4i2_body')
        for table in ['jobTable','importfileTable','exportfileTable','commentTable']:
            ele = ET.Element(table)
            top.append(ele)
        ele = ET.Element('job')
        if jobId is not None: ele.set('jobid',str(jobId))
        if jobNumber is not None: ele.set('jobnumber',str(jobNumber))
        top.xpath('./jobTable')[0].append(ele)
        return top

    def updateJob(self,key,value):
        if self._diagnostic: print('CJobDbBackup.updateJob',key,value)
        ele = self.body.xpath('./jobTable/job')[0]
        #print 'CJobDbBackup.updateJob',ele
        ele.set(key,str(value))

    def editImportFile(self,fileId=None,importId=None,sourceFileName=None,fileName=None,annotation=None,creationTime=None):
        if self._diagnostic: print('CJobDbBackup.editImportFile',fileId,importId,sourceFileName,fileName,annotation,creationTime)
        try:
            ele = self.body.xpath("./importfileTable/importfile[@importid='"+str(importId)+"']")[0]
        except:
            ele =  ET.Element('importfile')
            self.body.xpath('./importfileTable')[0].append(ele)
        if importId is not None: ele.set('importid',str(importId))
        if fileId is not None: ele.set('fileid',str(fileId))
        if sourceFileName is not None: ele.set('sourcefilename',str(sourceFileName))
        # This fileName should just be a base name (not full path) to help
        if fileName is not None:  ele.set('filename',str(fileName))
        if annotation is not None:  ele.set('annotation',str(annotation))
        if creationTime is not None:  ele.set('creationtime',str(creationTime))

    def editExportFile(self,fileId=None,exportId=None,exportFileName=None,creationTime=None):
        if self._diagnostic: print('CJobDbBackup.editExportFile',fileId,exportId,exportFileName,creationTime)
        try:
            ele = self.body.xpath("./exportfileTable/exportfile[@exportid='"+str(exportId)+"']")[0]
        except:
            ele =  ET.Element('exportfile')
            self.body.xpath('./exportfileTable')[0].append(ele)
        if exportId is not None: ele.set('exportid',str(exportId))
        if fileId is not None: ele.set('fileid',str(fileId))
        if exportFileName is not None:  ele.set('exportfilename',str(exportFileName))
        if creationTime is not None:  ele.set('creationtime',str(creationTime))

    def editComment(self,commentId=None,userName=None,timeOfComment=None,comment=None):
        if self._diagnostic: print('CJobDbBackup.editComment',commentId,userName,timeOfComment,comment)
        try:
            ele = self.body.xpath("./commentTable/comment[@commentid='"+str(commentId)+"']")[0]
        except:
            ele =  ET.Element('comment')
            self.body.xpath('./commentTable')[0].append(ele)
        if commentId is not None: ele.set('commentid',str(commentId))
        if userName is not None: ele.set('username',str(userName))
        if timeOfComment is not None: ele.set('timeofcomment',str(timeOfComment))
        if comment is not None: ele.set('comment',str(comment))

    def scrapeFromDb(self,db=None):
        saveItems = ['jobnumber','evaluation','preceedingjobid','jobtitle','creationtime','finishtime']
        try:
            ele = self.body.xpath('./jobTable/job')[0]
        except:
            ele = ET.Element('job')
            top.xpath('./jobTable')[0].append(ele)   # KJS : Worse than useless, this is messed up.
        info = db.getJobInfo(jobId=self.jobId,mode =saveItems )
        if info['evaluation'] in CCP4DbApi.JOB_EVALUATION_TEXT:
            info['evaluation'] = CCP4DbApi.JOB_EVALUATION_TEXT.index(info['evaluation'])
        else:
            info['evaluation'] = None
        for item in saveItems:
            if info[item] is not None: ele.set(item,str(info[item]))
        commentList = db.getComments(jobId=self.jobId)
        if self._diagnostic: print('CJobDbBackup.scrapeFromDb commentList',commentList)
        for cid,uName,machTime,comment in commentList:
            self.editComment(commentId=cid,userName=uName,timeOfComment=machTime,comment=comment)
        exportList = db.getExportFileInstances(jobId=self.jobId)
        for info in exportList:
            self.editExportFile(fileId=info['fileid'],exportId=info['exportid'],
                        exportFileName=info['exportfilename'],creationTime=info['creationtime'])
        importList = db.getImportFileInstances(jobId=self.jobId)
        for info in importList:
            self.editImportFile(fileId=info['fileid'],importId=info['importid'],sourceFileName=info['sourcefilename'],
                     creationTime=info['creationtime'],annotation=info['annotation'])


class CCopyJobDirectories:

    ERROR_CODES = {100 : { 'description' : 'Failed copying job directory'},
                   101 : { 'description' : 'Error creating temporary directory for files to send'},
                   102 : { 'description' : 'Error copying file to temporary directory'},
                   103 : { 'description' : 'Failed to find input file to copy'}}

    def __init__(self,projectId=None,jobIdList=None,targetDir=None,copyData=False):
        self.projectId = projectId
        self.jobIdList = jobIdList
        self.targetDir = targetDir
        self.copyData = copyData

    def copy(self):
        if not os.path.exists(self.targetDir): os.mkdir(self.targetDir)
        if not os.path.split(self.targetDir)[1] == 'CCP4_JOBS':
            self.targetDir = os.path.join(self.targetDir,'CCP4_JOBS')
            os.mkdir(self.targetDir)
        for jobId in self.jobIdList:
            self.copyJob(jobId)

    def copyJob(self,jobId):
        err = CErrorReport()
        jobDir = PROJECTSMANAGER().jobDirectory(jobId=jobId,projectId=self.projectId)
        targetJobDir = os.path.join(self.targetDir,os.path.split(jobDir)[1])
        if self.copyData:
            # Can just copy everything!
            try:
                shutil.copytree(jobDir,targetJobDir)
            except:
                err.append(self.__class__,100,details='From '+str(jobDir)+' to '+str(targetJobDir))
        else:
            # Copy log files and recurse through sub-directories
            err.append(self.copyJobLogFiles(jobDir,targetJobDir))
        if self.copyData:
            self.copyDataFiles(jobId)

    def copyDataFiles(self,jobId=None):
        # Need to copy input files from preceeding jobs
        err = CErrorReport()
        inputFileIdList = PROJECTSMANAGER().db().getFilesUsedInJobList(jobList=[jobId])
        if len(inputFileIdList)==0: return
        
        projectDir = PROJECTSMANAGER().db().getProjectDirectory(jobId=inputFileIdList[0][2])
        for inputFileId,fileName,sourceJobId,jobNumber,importId,sourceFilename in inputFileIdList:
            print('copyDataFiles',inputFileId,fileName,sourceJobId,jobNumber,importId,sourceFilename)
            # Replicate job directory structure for source data file
            jobDir = PROJECTSMANAGER().jobDirectory(jobId=sourceJobId,projectId=self.projectId)
            sourceFile = os.path.join( jobDir,fileName)
            if os.path.exists(sourceFile):
                targetJobDir = copy.deepcopy(self.targetDir)
                for jobNum in jobNumber.split('.'):
                    targetJobDir = os.path.join(targetJobDir,'job_'+str(jobNum))
                if not os.path.exists(targetJobDir): os.mkdir(targetJobDir)
                print('Copying input file',sourceFile)
                shutil.copy(sourceFile,targetJobDir)
            else:
                sourceImportFile = os.path.join(projectDir,'CCP4_IMPORTED_FILES',fileName)
                if os.path.exists(sourceImportFile):
                    print('Copying input file',sourceImportFile)
                    shutil.copy(sourceImportFile, os.path.join(self.targetDir,'CCP4_IMPORTED_FILES'))
                else:
                    err.append(self.__class__,103,str(sourceFile)+' or '+str(sourceImportFile))
        return err

    def copyJobLogFiles(self,jobDir,targetJobDir):
        err = CErrorReport()
        if not os.path.exists(targetJobDir):
            try:
                os.mkdir(targetJobDir)
            except:
                err.append(self.__class__,101, targetJobDir)
                return err
        for fileName in ['params.xml','input_params.xml','stderr.txt','stdout.txt','log.txt','com.txt']:
            filePath = os.path.join(jobDir,fileName)
            if os.path.exists(filePath):
                try:
                    shutil.copy(filePath,targetJobDir)
                except:
                    err.append(self.__class__,102,str(fileName)+' to '+str( targetJobDir))
        subJobDirList = glob.glob(os.path.join(jobDir,'job_*[0-9]'))
        for subJobDir in subJobDirList:
            if os.path.isdir(subJobDir):
                targetSubJobDir = os.path.join(targetJobDir,os.path.split(subJobDir)[1])
                err.append(self.copyJobLogFiles(subJobDir,targetSubJobDir))
        return err

def makeJobBackupForProject(projectName):
    db = PROJECTSMANAGER().db()
    projectId = db.getProjectId(projectName=projectName)
    projectDir = db.getProjectDirectory(projectId=projectId)
    jobIdList = db.getProjectJobList(projectId=projectId,topLevelOnly=True)
    #print 'CCP4Dbutils.makeJobBackupForProject',projectName,projectId,projectDir,jobIdList
    for jobId in jobIdList:
        makeJobBackup(jobId=jobId,projectName=projectName,db=db)

def makeJobBackup(jobId=None,projectName=None,db=None):
    if db is None:
        db = PROJECTSMANAGER().db()
    jobInfo=db.getJobInfo(jobId=jobId,mode=['projectname','projectid','jobNumber'])
    jobDirectory = db.jobDirectory(jobId=jobId)
    #print 'CCP4Dbutils.makeJobBackup',jobNumber,jobDirectory
    backup = CJobDbBackup(jobId=jobId,jobNumber=jobInfo['jobnumber'],jobDirectory=jobDirectory,projectName=jobInfo['projectname'],projectId=jobInfo['projectid'])
    backup.scrapeFromDb(db=db)
    backup.save()
