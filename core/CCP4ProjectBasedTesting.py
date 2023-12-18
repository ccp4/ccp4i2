from __future__ import print_function

"""
     CCP4ProjectBasedTesting.py: CCP4 GUI Project
     Copyright (C) 2014 STFC

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
   Liz Potterton Jan 2014 - Test ccp4i2 by rerunning projects and comparing output
"""
import os
import sys
import copy
import shutil
import time
from utils import startup
from dbapi import CCP4DbApi
from dbapi import CCP4DbUtils
from core import CCP4Container
from core import CCP4TaskManager
from core import CCP4Utils
import ccp4mg  # Ensure mmdb/hklfile etc available in testing
from PySide2 import QtCore
from lxml import etree
from xml.etree import ElementTree as ET
from core.CCP4ErrorHandling import *
from core.CCP4Modules import PROJECTSMANAGER,JOBCONTROLLER,QTAPPLICATION,PROCESSMANAGER

DIAGNOSTIC = False
MODES = {'runMode' : {'default' : 3 , 'allowed' : {0 : 'no testing' , 1 : 'run jobs', 2: 'run tests', 3: 'run jobs and tests'}},
         'verbosity' : {'default' : 1 , 'allowed' : {0 : 'fails only' , 1 : 'jobs and parameters tested', 3 : 'jobs, parameters and parameter values'}},
         'testSubJobs' : {'default' : 0 , 'allowed' : {0 : 'no' , 1 : 'for failed jobs', 2: 'for all jobs'}},
         'copyFiles' : {'default' : 0 , 'allowed' : {0 : 'no' , 1 : 'yes'}},
         'resetBaseline' : {'default' : 0 , 'allowed' : {0 : 'no' , 1 : 'yes'}}}


class Logger(object):
    '''
    For running stand-alone write to sys.stdout and to a log file.
    Within CCP4i2 gui should use the PRINTHANDLER instead
    '''
    def __init__(self, filename="ProjectBasedTesting.log"):
        self.terminal = sys.stdout
        self.log = open(filename, "a",1)

    def write(self, message):
        self.terminal.write(message)
        self.log.write(message)

    def close(self):
        self.log.close()


def run(sourceProjectList=[], outputDirectory=None, dbFile=None):
    '''Run a CProjectBasedTesting from a Python command line
    Set up Qt application (for event loop) and job controller so can run
    jobs as in graphical mode'''
    app = QTAPPLICATION(graphical=False)
    jc = startup.startJobController()
    jc.setDiagnostic(False)
    if dbFile is not None:
        jc.setDbFile(dbFile)
    PROJECTSMANAGER().startCheckForFinishedJobs()
    startTime = time.strftime('%y-%m-%d-%H-%M',time.localtime())
    print('Test results will be saved to: ' + os.path.join(outputDirectory, 'CCP4_test-' + startTime + '.log'))
    log = Logger(os.path.join(outputDirectory, 'test' + startTime + '.log'))
    t = CProjectBasedTesting(sourceProjectList=sourceProjectList, outputDirectory=outputDirectory,
                             dbFile=dbFile, parent=app, log=log, startTime=startTime)
    t.finished.connect(app.quit)
    t.runTests()
    sys.exit(app.exec_())


class CProjectBasedTesting(QtCore.QObject):
    reportUpdated = QtCore.Signal(str)
    finished = QtCore.Signal()

    ERROR_CODES = {101 : {'description' : 'Failed to open database'},
                   102 : {'description' : 'Project not found in database'},
                   103 : {'description' : 'Failed extracting database xml file from compressed file'},
                   104 : {'description' : 'Failed deleting old test project directory'},
                   201 : {'description' : 'Program command file not found for job'},
                   202 : {'description' : 'Differences in command file for job'},
                   203 : {'description' : 'Error occured testing for differences in command file for job'},
                   210 : {'description' : 'Job contains different number of subjobs'},
                   211 : {'description' : 'Failed creating CPluginScript object'}}

    def __init__(self, sourceProjectList=[], outputDirectory=None, dbFile=None, useCurrentDb=False, parent=None, log=None, startTime=None, selectedJobs=None, logXmlRoot=None, **kw):
        QtCore.QObject.__init__(self, parent)
        self.report = ''
        self.sourceProjectList = sourceProjectList
        self.outputDirectory = outputDirectory
        self.dbFile = dbFile
        if selectedJobs is None:
            self.selectedJobs = None
        else:
            self.selectedJobs = []
            for item in selectedJobs.split(','):
                rnge = item.split('-')
                if len(rnge) == 1:
                    self.selectedJobs.append(rnge[0])
                elif len(rnge) == 2:
                    try:
                        for i in range(int(rnge[0]), int(rnge[1])):
                            self.selectedJobs.append(str(i))
                    except:
                        pass
        #print 'selectedJobs', self.selectedJobs
        if startTime is None:
            self.startTime = time.strftime('%y-%m-%d-%H-%M', time.localtime())
        else:
            self.startTime = startTime
        self.sourceProjectId = None
        self.sourceProjectName = None
        self.testProjectName = None
        self.testProjectDir = None
        self.testProjectId = None
        self.sourceProjectJobs = []
        self.runningJobId = None
        self.sourceJobInfo = None
        if log is None:
            self.log = None
            self.logByProject = True
        else:
            self.log = log
            self.logByProject = False
        self.logFiles = []
        
        #MN added a structured element to the report logging - allows better interface to pytest or unittest
        if logXmlRoot is None: self.logXmlRoot = None
        else: self.logXmlRoot = logXmlRoot
        
        if self.outputDirectory is None:
            self.outputDirectory = os.getcwd()
        if self.dbFile is None:
            self.dbFile = os.path.join(self.outputDirectory, 'db.sqlite')
        self.workDir = os.path.join(self.outputDirectory, 'tmp_dir')
        if not os.path.exists(self.workDir):
            os.mkdir(self.workDir)
        for name,defn in list(MODES.items()):
            setattr(self, name, kw.get(name, defn['default']))
        if self.resetBaseline:
            self.copyFiles = 1
        if not useCurrentDb:
            try:
                PROJECTSMANAGER().db().openDb(fileName=self.dbFile, createDb=True, diagnostic=False)
            except:
                raise CException(self.__class__, 101, self.dbFile)
            else:
                JOBCONTROLLER().setDbFile(self.dbFile)
        PROJECTSMANAGER().db().jobFinished.connect(self.handleJobFinished)
        print('Project based testing set up using database:', self.dbFile)

    def runTests(self):
        #JOBCONTROLLER().setDiagnostic(False)
        self.nextJob()

    @QtCore.Slot(dict)
    def handleJobFinished(self, args):
        if DIAGNOSTIC:
            print('CProjectBasedTesting.handleJobFinished', args)
        jobId = args.get('jobId', None)
        if jobId is None or jobId != self.runningJobId:
            return
        self.analyseJob(self.sourceJobInfo['jobid'], jobId)
        if len(self.logFiles) > 0:
            self.reportUpdated.emit(self.logFiles[0])
        self.nextJob()

    @QtCore.Slot()
    def nextJob(self):
        if DIAGNOSTIC:
            print('nextJob', len(self.sourceProjectJobs))
        while len(self.sourceProjectJobs) == 0:
            #self.finishProject()
            #start on the next project
            if len(self.sourceProjectList) == 0:
                self.finish()
                return
            sourceProject0 = self.sourceProjectList.pop(0)
            if self.logXmlRoot is not None:
                self.currentProjectNode = ET.SubElement(self.logXmlRoot,'Project',projectPath=str(sourceProject0))
            if sourceProject0.endswith('zip'):
                sourceProject = self.unpackProject(sourceProject0)
            else:
                sourceProject = sourceProject0
            self.setupProject(sourceProject)
            if self.logByProject:
                if self.log is not None:
                    self.log.close()
                self.logFiles.append(os.path.normpath(os.path.join(self.testProjectDir, 'CCP4_test-' + self.startTime + '.log')))
                self.log = open(self.logFiles[-1], 'a', 1)
            self.log.write('\n\nPROJECT: ' + self.testProjectName + ' saved to ' + self.testProjectDir + '\n')
        self.sourceJobInfo = self.sourceProjectJobs.pop(0)
        if DIAGNOSTIC:
            print('CProjectBasedTesting.nextJob sourceJobInfo', self.sourceJobInfo)
        if self.selectedJobs is not None:
            while self.sourceJobInfo['jobnumber'] not in self.selectedJobs:
                #print 'CProjectBasedTesting.nextJob jobnumber', self.sourceJobInfo['jobnumber']
                if len(self.sourceProjectJobs) > 0:
                    self.sourceJobInfo = self.sourceProjectJobs.pop(0)
                else:
                    self.finish()
                    return
        if hasattr(self, 'currentProjectNode') and self.currentProjectNode is not None:
            self.currentJobNode = ET.SubElement(self.currentProjectNode,'Job',number=self.sourceJobInfo['jobnumber'])
        job = CCP4DbUtils.COpenJob(projectId=self.testProjectId)
        ret = job.createJob(taskName=self.sourceJobInfo['taskname'], jobNumber=self.sourceJobInfo['jobnumber'], cloneJobId=self.sourceJobInfo['jobid'], copyInputFiles=self.copyFiles)
        if self.copyFiles:
            self.copyImportedFiles(fromJobId=self.sourceJobInfo['jobid'], toJobId=job.jobId)

        sourceProjectId = PROJECTSMANAGER().db().getJobInfo(self.sourceJobInfo['jobid'],'projectid')
        targetProjectId = PROJECTSMANAGER().db().getJobInfo(job.jobId,'projectid')
        sourceDir = PROJECTSMANAGER().db().getProjectDirectory(projectId=sourceProjectId)
        targetDir = PROJECTSMANAGER().db().getProjectDirectory(projectId=targetProjectId)
        if os.path.isdir(os.path.join(sourceDir,"CCP4_TEST_FILES")):
            shutil.copytree(os.path.join(sourceDir,"CCP4_TEST_FILES"), os.path.join(targetDir,"CCP4_TEST_FILES"))

        #Try running the job
        self.log.write('\nJob: ' + str(job.jobNumber) + ' ' + str(job.taskName) + '\n')
        self.runningJobId=copy.deepcopy(job.jobId)
        PROJECTSMANAGER().db().resetLastJobNumber(self.testProjectId)
        if CCP4TaskManager.TASKMANAGER().isInternalPlugin(job.taskName):
          PROJECTSMANAGER().runInternalTask(jobId = self.runningJobId,projectId=self.testProjectId,
                               taskName = job.taskName)
        else:
          PROJECTSMANAGER().updateJobStatus(jobId=self.runningJobId, status=CCP4DbApi.JOB_STATUS_QUEUED)
          QTAPPLICATION().doCheckForFinishedJobs.emit()

    def unpackProject(self, compressedFile):
        from qtcore import CCP4Export
        # Extract database xml file from compressedFile
        try:
            xmlFile = PROJECTSMANAGER().extractDatabaseXml(compressedFile)
        except CException as e:
            raise e
        except Exception as e:
            raise CException(self.__class__, 103, compressedFile)
        # Load xml into CDbXml object and thence to db temporary tables
        dbImport = CCP4DbApi.CDbXml(db=PROJECTSMANAGER().db(), xmlFile=xmlFile)
        header = dbImport.headerInfo(load=True)
        #timeTag = '_'+str(int(time.time()))
        timeTag = ''
        dirName = os.path.join(self.outputDirectory, dbImport.projectName + timeTag)
        dbImport.projectDirectory = dirName
        dbImport.createProject()
        dbImport.createTempTables()
        dbImport.loadTempTable()
        # Unpack project files from the tar file (possibly in separate thread) 
        # Pass import thread dbImport to enable query database and flagging loaded jobs/files
        try:
            if os.path.exists(dirName):
                shutil.rmtree(dirName)
        except:
            raise CException(self.__class__, 104, dirName)
        os.mkdir(dirName)
        self.importThread = CCP4Export.ImportProjectThread(self, projectDir=dirName, compressedFile=compressedFile, dbImport=dbImport)
        errReport = self.importThread.run()
        #dbImport.listTempFiles()
        dbImport.cleanupTempTables()
        # Seem to need to disable the exclude mechanism 
        dbImport.importTempTables(testExl=False)
        dbImport.removeTempTables()
        if DIAGNOSTIC:
            print('CProjectBasedTesting.unpackProject db file', PROJECTSMANAGER().db()._fileName)
        projectName = copy.deepcopy(dbImport.projectName)
        self.lastCompressedFile = copy.deepcopy(compressedFile)
        return projectName

    def setupProject(self,sourceProject):
        if DIAGNOSTIC: print('CProjectBasedTesting,setupProject',sourceProject)
        try:
            self.sourceProjectId = PROJECTSMANAGER().db().getProjectId(projectName=sourceProject)
            self.sourceProjectName=sourceProject
        except:
            try:
                self.sourceProjectName = PROJECTSMANAGER().db().getProjectInfo(projectId=sourceProject, mode='projectname', checkPermission=False)
                self.sourceProjectId = sourceProject
            except:
                raise CException(self.__class__, 102, sourceProject)
        ifExists = True
        while ifExists:
            self.testProjectName = self.sourceProjectName + '_rerun_' + self.startTime
            self.testProjectDir =os.path.normpath(os.path.join(self.outputDirectory, self.testProjectName))
            ifExists = os.path.exists(self.testProjectDir)
        #Delete any old test project
        #PROJECTSMANAGER().deleteProject(projectName=self.testProjectName,deleteDirectory=True)
        # NB use PROJECTSMANAGER method that create directory
        self.testProjectId = PROJECTSMANAGER().createProject(projectName=self.testProjectName, projectPath=self.testProjectDir)
        if DIAGNOSTIC:
            print('CProjectBasedTesting,setupProject testProjectId', self.testProjectId)
            print('CProjectBasedTesting db file', PROJECTSMANAGER().db()._fileName)
        self.sourceProjectJobs = PROJECTSMANAGER().db().getProjectJobListInfo(projectId=self.sourceProjectId, mode=['jobid', 'jobnumber', 'taskname'],
                                                                              topLevelOnly=True, jobStatus=[CCP4DbApi.JOB_STATUS_FINISHED, CCP4DbApi.JOB_STATUS_UNSATISFACTORY] )
        #print 'CProjectBasedTesting.setupProject',self.sourceProjectId,self.sourceProjectJobs
        return len(self.sourceProjectJobs)

    def copyImportedFiles(self, fromJobId, toJobId):
        importedFileInfoList = PROJECTSMANAGER().db().getImportFileInstances(jobId=fromJobId, brief=False)
        for info in importedFileInfoList:
            targetFileId = PROJECTSMANAGER().db().createFile(jobId=toJobId,fileTypeId=info['filetypeid'], sourceFileName=info['sourcefilename'],
                                                             sourceFileId=info['sourcefileid'], sourceFileAnnotation=info['annotation'],
                                                             importNumber=info['importnumber'], reference=info['reference'],
                                                             jobParamName=info['jobparamname'], baseName=info['filename'], pathFlag=info['pathflag'],
                                                             annotation=info['annotation'], fileContent=info['filecontent'], subType=info['filesubtype'])
            targetFile = PROJECTSMANAGER().db().getFullPath(fileId=targetFileId)
            sourceFile = PROJECTSMANAGER().db().getFullPath(fileId=info['fileid'])
            if DIAGNOSTIC:
                print('Copying imported file:',sourceFile,'to',targetFile)
            try:
                shutil.copyfile(sourceFile, targetFile)
            except:
                print('ERROR copying imported file:', sourceFile, 'to', targetFile)

    def analyseJob(self,sourceJobId, testJobId):
        self.reportJobStatus(sourceJobId, expected=True)
        self.reportJobStatus(testJobId)
        self.listPluginErrorReport(testJobId)
        severity = self.runJobTest(sourceJobId, testJobId)
        if self.testSubJobs == 2 or (self.testSubJobs == 1 and severity >= SEVERITY_WARNING):
            self.runSubJobTests(sourceJobId, testJobId)

    def runSubJobTests(self, sourceJobId, testJobId, testSubJobs=False):
        sourceSubJobList = PROJECTSMANAGER().db().getJobInfo(sourceJobId, 'childjobs')
        testSubJobList = PROJECTSMANAGER().db().getJobInfo(testJobId, 'childjobs')
        nSubJobs = min(len(sourceSubJobList),len(testSubJobList))
        if nSubJobs > 0:
            self.log.write('Testing sub-jobs..\n')
            for idx in range(nSubJobs):
                info = PROJECTSMANAGER().db().getJobInfo(testSubJobList[idx], ['jobnumber','taskname'])
                self.log.write('\nJob: ' + str(info['jobnumber']) + ' ' + str(info['taskname']) + '\n')
                self.runJobTest(sourceSubJobList[idx], testSubJobList[idx], testSubJobs)

    def runJobTest(self, sourceJobId, testJobId, testSubJobs=False):
        report = CErrorReport()
        taskName = PROJECTSMANAGER().db().getJobInfo(sourceJobId, 'taskname')
        cls = CCP4TaskManager.TASKMANAGER().getPluginScriptClass(taskName)
        if DIAGNOSTIC:
            print('CCP4ProjectBasedTesting.runJobTest cls', taskName, cls)
        try:
            sourcePlugin = cls(parent=self, name='source', workDirectory=self.workDir)
        except:
            report.append(self.__class__, 211, 'Taskname: ' + str(taskName) + ' class: ' + str(cls))
        else:
            #Beware dummy plugins such as crank2 sub-tasks may not have a def file and so do not have a container
            if sourcePlugin.container is not None:
                sourcePlugin.container.loadDataFromXml(fileName=PROJECTSMANAGER().db()._makeJobFileName(jobId=sourceJobId, mode='PARAMS'))
                testContainer = CCP4Container.CContainer()
                testContainer.loadDataFromXml(fileName=PROJECTSMANAGER().db()._makeJobFileName(jobId=testJobId, mode='PARAMS'))
                report.extend(sourcePlugin.assertSame(testContainer, diagnostic=DIAGNOSTIC))
            else:
                from core import CCP4PluginScript
                report.append(CErrorReport(CCP4PluginScript.CPluginScript, 312, details='For taskname:' + taskName))
        if self.verbosity == 0:
            minSeverity = SEVERITY_ERROR
        else:
            minSeverity=SEVERITY_OK
        self.log.write(report.report(ifStack=False, mode=1, minSeverity=minSeverity))
        self.log.write('\n')
        if hasattr(self, 'currentJobNode'):
            self.currentJobNode.append(report.getEtree())
        if testSubJobs:
            self.runSubJobTests(sourceJobId, testJobId)
        return report.maxSeverity()

    def reportJobStatus(self, jobId, expected=False):
        status = PROJECTSMANAGER().db().getJobInfo(jobId, 'status')
        if not expected:
            self.log.write('Job database status: ' + status)
            if hasattr(self, 'currentJobNode'):
                ET.SubElement(self.currentJobNode,'Status').text = str(status)
        else:
            self.log.write('Job database expected status: ' + status)
            if hasattr(self, 'currentJobNode'):
                ET.SubElement(self.currentJobNode,'ExpectedStatus').text = str(status)

    def listPluginErrorReport(self, jobId):
        try:
            tree = CCP4Utils.openFileToEtree(PROJECTSMANAGER().makeFileName(jobId=jobId, mode='DIAGNOSTIC'))
            report = CErrorReport()
            report.setEtree(tree)
            #print 'listPluginErrorReport',tree
        except:
            self.log.write('Failed loading diagnostic file for job')
        else:
            #print 'listPluginErrorReport', PROJECTSMANAGER().makeFileName(jobId=jobId, mode='DIAGNOSTIC'), len(report)
            if len(report) > 0:
                self.log.write(report.report(ifStack=True, mode=0, minSeverity=SEVERITY_WARNING))
                if hasattr(self, 'currentJobNode'):
                    ET.SubElement(self.currentJobNode,'Warnings').text = report.report(ifStack=True, mode=0, minSeverity=SEVERITY_WARNING)
                    if report.maxSeverity() >= SEVERITY_ERROR:
                        ET.SubElement(self.currentJobNode,'Errors').text = report.report(ifStack=True, mode=0, minSeverity=SEVERITY_ERROR)
                        

    def finishProject(self):
        from qtcore import CCP4Export
        #print 'finishProject',self.resetBaseline,self.testProjectDir,getattr(self,'lastCompressedFile',None)
        if (not self.resetBaseline) or (self.testProjectDir is None) or (getattr(self, 'lastCompressedFile', None) is None):
            return
        #print 'finishProject resetBaseline',self.testProjectDir
        # Want to export the test project to a zip file and copy to the the original test zip
        # but backup the original compressed file
        os.rename(self.lastCompressedFile, self.lastCompressedFile + '.backup_' + time.strftime('%y-%m-%d-%H-%M', time.localtime()))
        dbxml = os.path.join(self.testProjectDir, 'CCP4_TMP', 'DATABASE' + str(int(time.time())) + '.db.xml')
        print('Creating temporary database xml file in:',dbxml)
        errReport = PROJECTSMANAGER().db().exportProjectXml(self.testProjectId, fileName=dbxml)
        self.exportThread = CCP4Export.ExportProjectThread(self, projectDir=self.testProjectDir, dbxml=dbxml, target=self.lastCompressedFile)
        self.lastCompressedFile = None
        self.exportThread.finished.connect(self.nextJob)
        self.exportThread.wait()
        print('after exportThread.wait')

    def finish(self):
        if self.log is not None:
            self.log.write('\nTESTS FINISHED\n')
            self.log.close()
            if len(self.logFiles) > 0:
                self.reportUpdated.emit(self.logFiles[0])
        self.finished.emit()

    def compareComFiles(self, sourceJobId, testJobId):
        report = CErrorReport()
        sourceComFile = PROJECTSMANAGER().db()._makeJobFileName(jobId=sourceJobId, mode='COM')
        testComFile = PROJECTSMANAGER().db()._makeJobFileName(jobId=testJobId, mode='COM')
        if os.path.exists(sourceComFile):
            if not os.path.exists(testComFile):
                jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=testJobId, mode=['jobnumber','taskname'])
                report.append(self.__class__, 201, 'Job ' + jobInfo.get('jobnumber')+' ' + jobInfo.get('taskname'))
            else:
                sourceText = CCP4Utils.readFile(sourceComFile)
                testText = CCP4Utils.readFile(sourceComFile)
                dmp =  diff_match_patch.diff_match_patch()
                try:
                    diffs = dmp.diff_main(text1, text2) # KJS : Small problem here... no text1, no text2
                except:
                    report.append(self.__class__, 203, name = 'Job ' + jobInfo.get('jobnumber') + ' ' + jobInfo.get('taskname'))
                else:
                    if len(diffs) > 0:
                        report.append(self.__class__,202, name = 'Job ' + jobInfo.get('jobnumber')+' '+jobInfo.get('taskname'), details =str(diffs))
        sourceChildJobs = PROJECTSMANAGER().db().getJobInfo(sourceJobId, mode='childjobs')
        testChildJobs = PROJECTSMANAGER().db().getJobInfo(testJobId, mode='childjobs')
        if len(sourceChildJobs) == 0 and len(testChildJobs) == 0:
            return report
        # Might as well test the sub-jobs run here 
        if len(sourceChildJobs) != len(testChildJobs):
            report.append(self.__class__, 210, name = 'Job ' + jobInfo.get('jobnumber') + ' ' + jobInfo.get('taskname'))
            for idx in range(min(len(sourceChildJobs), len(testChildJobs))):
                sourceSubtask = PROJECTSMANAGER().db().getJobInfo(sourceChildJobs[idx], mode=['taskname','jobnumber'])
                testSubtask = PROJECTSMANAGER().db().getJobInfo(testChildJobs[idx], mode=['taskname','jobnumber'])
                if sourceSubtask.get('jobnumber') != testSubtask.get('jobnumber') or sourceSubtask.get('taskname') != testSubtask.get('taskname'):
                    report.append(self.__class__, 211, name = 'Job ' + jobInfo.get('jobnumber') + ' ' + jobInfo.get('taskname'),
                                  details='Source:' + sourceSubtask.get('jobnumber') + ' '+sourceSubtask.get('taskname') + ' test:'+
                                  testSubtask.get('jobnumber') + ' ' + testSubtask.get('taskname') )
                else:   
                    self.compareComFiles(sourceChildJobs[idx],testChildJobs[idx])

    def evaluate(self):
        if not os.path.exists(os.path.join(self.sourceProjectDir(), 'CCP4_TEST_SYSTEM', 'xml')):
            self.runTestSys('import')
        if not os.path.exists(os.path.join(self.sourceProjectDir(), 'CCP4_TEST_SYSTEM', 'runs')):
            self.runTestSys('run --simulate')
        self.runTestSys('eval')

    def sourceProjectDir(self):
        return PROJECTSMANAGER().db().getProjectDirectory(projectId=self.sourceProjectId)

    def handleRunTestSysFinished(self, pid):
        print('handleRunTestSysFinished', pid)

    def runTestSys(self,argList):
        #exePath = os.path.join(os.environ('TESTSYSSRC'),'testsys.py')
        exePath = '/Users/lizp/Desktop/pavols_test_sys/testsys.py'
        projectDir = self.sourceProjectDir()
        PROCESSMANAGER().startProcess(exePath, arg=argList, editEnv=[['TESTSYS', os.path.join(projectDir, 'CCP4_TEST_SYSTEM')]], handler=[self.handleRunTestSysFinished])


class BuildTestSuite:

    ERROR_CODES = {210 : {'description' : 'Cannot create test suite - not all jobs are successfully completed'},
                   211 : {'description' : 'Cannot create test suite - failed creating sub-directory in project directory'},
                   212 : {'description' : 'Overwriting existing job in test directory', 'severity' : SEVERITY_WARNING},
                   213 : {'description' : 'Failed to create job directory in test directory'}}

    def __init__(self, projectId):
        self.projectId = projectId

    def run(self):
        self.errorReport = CErrorReport()
        projectInfo = PROJECTSMANAGER().db().getProjectInfo(projectId=self.projectId)
        #print 'PROJECTSMANAGER().makeTestSuite projectId',projectId,projectInfo
        unfinishedJobs = PROJECTSMANAGER().db().getProjectUnfinishedJobs(projectId=self.projectId)
        if len(unfinishedJobs) > 0:
            return CErrorReport(self.__class__, 210, details='for project: ' + str(projectInfo.get('projectname', 'UNKNOWN')), stack=False)
        self.unittestDir = os.path.normpath(os.path.join(projectInfo['projectdirectory'], 'CCP4_TEST_SYSTEM'))
        if not os.path.exists(self.unittestDir):
            try:
                os.mkdir(self.unittestDir)
            except:
                return CErrorReport(self.__class__, 211,details='in project directory: ' + str(projectInfo.get('projectdirectory', 'UNKNOWN')), stack=False)
        jobList = PROJECTSMANAGER().db().getProjectJobListInfo(mode=['jobid', 'jobnumber', 'taskname'], projectId=self.projectId, topLevelOnly=True)
        # Create the template test for each job
        #print 'makeTestSuite jobList', jobList
        for jobInfo in jobList:
            programName = CCP4TaskManager.TASKMANAGER().getTaskAttribute(jobInfo['taskname'], 'TASKCOMMAND', script=True)
            if programName is not None:
                self.addTest(jobInfo['jobnumber'], programName)
        return self.errorReport

    def addTest(self, jobNumber, programName):
        sourceDir = os.path.join(CCP4Utils.getTestSysDir(), 'import', programName)
        #print 'BuildTestSuite.addTest', sourceDir
        if not os.path.exists(sourceDir):
            return
        targetDir = os.path.join(self.unittestDir, 'import', programName)
        if not os.path.exists(targetDir):
            shutil.copytree(sourceDir, targetDir)
        path = os.path.join(targetDir, 'default_test', 'job_' + jobNumber)
        if os.path.exists(path):
            self.errorReport.append(self.__class__, 212, path)
            os.rmdir(path)
        try:
            os.mkdir(path)   # KJS : Not clear what this function is supposed to be doing ...check.
        except:
            self.errorReport.append(self.__class__, 213, path)

