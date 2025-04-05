"""
Copyright (C) 2010 University of York
Liz Potterton - May 2010 - Quick version with hard-coded project - pending database
"""

# NB  load methods from ccp4mg and inconsistent with data types in the init

import copy
import datetime
import functools
import glob
import os
import re
import shutil
import sys
import tarfile
import tempfile
import time
import zipfile

from lxml import etree
from PySide2 import QtCore

from . import CCP4File
from . import CCP4Utils
from .CCP4ErrorHandling import CErrorReport, CException, Severity
from .CCP4Modules import JOBCONTROLLER, DEMODATAMANAGER, PREFERENCES
from .CCP4QtObject import CObject
from .CCP4Utils import getCCP4I2Dir, getDotDirectory, getHOME
from .CCP4WarningMessage import warningMessage


    if CProjectsManager.insts is None:
        CProjectsManager.insts = CProjectsManager()
    return CProjectsManager.insts


#----------------------------------------------------------------------
class CProjectsManager(CObject):
#----------------------------------------------------------------------

    projectsListChanged = QtCore.Signal()
    doCheckForFinishedJobs = QtCore.Signal()
    jobUpdated = QtCore.Signal(str)

    #beware handling of directories does not allow for different context
    # ? but do we have different CProjectsManager for different context?

    insts = None
    CHECKINTERVAL = 1000
    ERROR_CODES = {101 : {'severity' : Severity.ERROR, 'description' : 'Failed creating CCP4I2_TEST project in'},
                   102 : {'severity' : Severity.ERROR, 'description' : 'CProjectsManager.addProject project of that name exists'},
                   103 : {'severity' : Severity.ERROR, 'description' : 'CProjectsManager.addProject imput project wrong type'},
                   104 : {'severity' : Severity.ERROR,'description' : 'Creating project, no name given for project'},
                   105 : {'severity' : Severity.ERROR,'description' : 'CProjectsManager.createProject Error creating project'},
                   106 : {'severity' : Severity.ERROR,'description' : 'CProjectsManager.createProject Error setting project directory'},
                   107 : {'severity' : Severity.ERROR,'description' : 'CProjectsManager.createProject Error making project directory'},
                   108 : {'severity' : Severity.ERROR,'description' : 'CProjectsManager.openDirectoriesDef  No HOME/USERPROFILE directory'},
                   109 : {'severity' : Severity.ERROR,'description' : 'CProjectsManager.deleteProject no project of that name'},
                   110 : {'severity' : Severity.ERROR,'description' : 'CProjectsManager.openDirectoriesDef Error opening directories.def file'},
                   111 : {'severity' : Severity.ERROR,'description' : 'CProjectsManager.write_ccp4i_directories No HOME/USERPROFILE directory, Error saving directories.def '},
                   112 : {'severity' : Severity.ERROR,'description' : 'Attempting to create project in directory that is used by another project.'},
                   113 : {'description' : 'No directory name provided for making directory'},
                   114 : {'description' : 'Error creating directory'},
                   115 : {'description' : 'Error creating sub-directory in project directory'},
                   116 : {'description' : 'Error creating database config file in project directory'},
                   117 : {'description' : 'Unknown error loading from XML database'},
                   118 : {'description' : 'Failed attempting to write job control file'},
                   120 : {'description' : 'Error reading CCP4I2 projects file'},
                   130 : {'description' : 'Error calling makeFileName without job number'},
                   140 : {'description' : 'Error calling deleteJob without jobId'},
                   141 : {'description' : 'Error deleting job - failed to delete directories - possibly files still in use'},
                   142 : {'description' : 'CProjectsManager.jobDirectory - directory does not exist'},
                   143 : {'description' : 'CProjectsManager failed saving imported files from a deleted job'},
                   144 : {'description' : 'CProjectsManager failed interpreting job directory'},
                   151 : {'description' : 'CProjectsManager.setOutputFileNames error accessing container'},
                   152 : {'description' : 'CProjectsManager.setOutputFileNames error for parameter'},
                   153 : {'description' : 'CProjectsManager.setInputFileNames error accessing container'},
                   154 : {'description' : 'CProjectsManager.setInputFileNames error for parameter'},
                   160 : {'description' : 'Importing file failed to copy'},
                   161 : {'description' : 'Importing file failed to record in database'},
                   162 : {'description' : 'Unknown error loading from params file for CProjectsManager.importData'},
                   171 : {'description' : 'Error searching database for project jobs'},
                   172 : {'description' : 'Error saving to project compressed file'},
                   173 : {'description' : 'Error deleting temporary database xml file'},
                   174 : {'description' : 'Error creating temporary database xml file'},
                   175 : {'description' : 'Error creating project compressed file'},
                   176 : {'description' : 'Error opening project compressed file'},
                   177 : {'description' : 'Error reading project compressed file'},
                   178 : {'description' : 'Compressed file format not recognised (expecting tar.gz or zip)'},
                   180 : {'description' : 'Error reading input_params.xml file to find list of input files to job'},
                   181 : {'description' : 'Error extracting an input file path from input_params.xml'},
                   182 : {'description' : 'Input file for job does not exist'},
                   200 : {'severity' : Severity.OK,'description' : 'Successfully deleted project directory'},
                   201 : {'description' : 'Error attempting to delete project directory'},
                   202 : {'severity' : Severity.OK,'description' : 'Successfully deleted project from database'},
                   203 : {'description' : 'Error attempting to delete project from database'},
                   204 : {'description' : 'No projectId provided for delete project from database'},
                   210 : {'description' : 'Failed to delete temporary file created in exporting MTZ'}}

    def __init__(self, parent=None, database=None, **kw):
        CObject.__init__(self, parent)
        self._db = None
        self._dbTaskTitlesLoaded = False
        if database is not None:
            self.setDatabase(database)
        self.blockExit = False

    def setDatabase(self, database):
        if self._db is not None:
            self._db.close()
            self._db = None
            self._dbTaskTitlesLoaded = False
        self._db = database
        self.lastJobFinishCheckTime = time.time()
        self.lastJobStartedCheckTime = time.time()
        #self.initialiseDirectories()
        #print 'PROJECTSMANAGER.setDatabase',self._db
        self._db.jobFinished.connect(self.backupDBIfTopLevel)

    def loadDbTaskTitles(self):
        from . import CCP4TaskManager
        # Load task titles on startup to enable text searching
        if self._dbTaskTitlesLoaded:
            return
        taskTitlesList = CCP4TaskManager.TASKMANAGER().getAllShortTitles()
        self.db().loadTaskTitleTable(taskTitlesList)
        self._dbTaskTitlesLoaded = True

    def startCheckForFinishedJobs(self):
        self.db().jobToDelete.connect(self.handleJobToDelete)
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self._checkForFinishedJobs)
        self.timer.start(CProjectsManager.CHECKINTERVAL)

    def db(self, label=None):
        return self._db

    @QtCore.Slot(dict)
    def backupDBIfTopLevel(self,args):
        #We do not back up if we can verify that this is a sub job
        if hasattr(args,"has_key") and "parentJobId" in args and args["parentJobId"] is not None:
            print("Not backing up after sub-job")
            return
        if hasattr(args,"has_key") and "parentJobId" not in args:
            print("Not backing up after possible sub-job")
            return
        projectInfo = PROJECTSMANAGER().db().getProjectInfo(projectId=args["projectId"])
        dbxml = os.path.join(projectInfo['projectdirectory'],'DATABASE.db.xml')
        print("Exporting project meta-data as XML file:")
        print(dbxml)
        PROJECTSMANAGER().db().exportProjectXml(args["projectId"],fileName=dbxml)
        print("Backing up after top-level job finished")
        self.backupDB()

    def backupDB(self):
        print("----------------------------------------backupDB----------------------------------------")
        def datetimesort(k1,k2):
            t1 = datetime.datetime.strptime(k1.lstrip("database_sqlite_backup-"),"%d%m%Y-%H%M%S")
            t2 = datetime.datetime.strptime(k2.lstrip("database_sqlite_backup-"),"%d%m%Y-%H%M%S")
            if t1 < t2:
                return -1
            if t1 > t2:
                return 1
            if t1 == t2:
                return 0
            return None
        n = datetime.datetime.now()
        basefname = datetime.datetime.strftime(n,"database_sqlite_backup-%d%m%Y-%H%M%S")
        fname = os.path.join(CCP4Utils.getDotDirectory(),'db',basefname)
        self.dbOrigDir =  os.path.join(CCP4Utils.getDotDirectory(),'db')
        self.db().backupDB(fname)
        #Now clean up the old ones.
        backups = glob.glob(os.path.join(self.dbOrigDir,"database_sqlite_backup-*[0-9]"))
        baseBackups = [ os.path.basename(f) for f in backups ]
        sortedBaseBackups = sorted(baseBackups,key=functools.cmp_to_key(datetimesort),reverse=True)
        sortedFullBaseBackups = [ os.path.join(self.dbOrigDir,f) for f in sortedBaseBackups ]
        print(len(sortedFullBaseBackups))
        backupsToDelete = sortedFullBaseBackups[30:]
        print(backupsToDelete)
        for b in backupsToDelete:
            print(b)
            os.remove(b)

    def backupDBXML(self):
        proj_dir_list0=PROJECTSMANAGER().db().getProjectDirectoryList()
        root = etree.Element("ProjectRoots")
        for projectInfo in proj_dir_list0:
            projectRoot = etree.Element("project")
            projectRoot.text =  projectInfo[2]
            root.append(projectRoot)
        dbListBackupName = os.path.join(CCP4Utils.getDotDirectory(),'projectList-backup.xml')
        dbListBackupFile = open(dbListBackupName,"w+")
        CCP4Utils.writeXML(dbListBackupFile,etree.tostring(root,pretty_print=True))
        dbListBackupFile.close()
        print("Backed up list of projects to",dbListBackupName)

    def initialiseDirectories(self):
        self._db.setDirectoryAlias('CCP4I2_TOP', getCCP4I2Dir())

    def Exit(self):
        self.timer.stop()
        #sys.__stdout__.write('CProjectsManager.Exit blockExit '+str(self.blockExit)+'\n');sys.__stdout__.flush()

    def makeDefaultProject(self, projectName='CCP4I2_TEST', projectPath=None):
        try:
            defaultProject = self.db().getProjectId(projectName=projectName)
        except:
            if projectPath is None:
                projectPath = os.path.join(getHOME(), projectName)
            defaultProject = self.createProject(projectName, projectPath)
        return defaultProject

    def projectStatus(self,projectName=None):
        try:
            pid = self._db.getProjectId(projectName=projectName)
        except:
            return 1  # project name not recognised
        '''
        # This is probably not appropriate here
        userName = CCP4Utils.getUserId()
        self._db.getProjectPermission(projectId=pid,userName=userName)
        '''
        info = self._db.getProjectInfo(projectId=pid)
        if not os.path.exists(info['projectdirectory']):
            return 2 # Project directory not found
        return 0

    def createProject(self, projectName=None, projectPath=None):
        self.makeProjectDirectory(projectPath)
        projectID = self._db.createProject(projectName=projectName, projectDirectory=projectPath)
        self.projectsListChanged.emit()
        self.backupDB()
        return projectID

    def deleteProject(self, projectId=None, projectName=None, deleteDirectory=False):
        e = CException()
        if projectId is None and projectName is not None:
            try:
                projectId = self.db().getProjectId(projectName=projectName)
            except:
                pass
        if projectId is None:
            return  CException(self.__class__, 204)
        projectInfo = self.db().getProjectInfo(projectId=projectId)
        try:
            self._db.deleteProject(projectId=projectId)
            e = CException(self.__class__, 202, projectInfo['projectname'])
        except:
            e = CException(self.__class__, 203, projectInfo['projectname'])
        if deleteDirectory:
            try:
                subDirList = glob.glob(os.path.join(projectInfo['projectdirectory'], 'CCP4_*'))
                for subDir in subDirList:
                    if os.path.isfile(subDir):
                        os.remove(subDir)
                    else:
                        shutil.rmtree(subDir, True)
                if os.path.exists(os.path.join(projectInfo['projectdirectory'], 'DATABASE.db.xml')):
                    os.remove(os.path.join(projectInfo['projectdirectory'], 'DATABASE.db.xml'))
                if len(glob.glob(os.path.join(projectInfo['projectdirectory'], '*'))) == 0:
                    shutil.rmtree(projectInfo['projectdirectory'], True)
                e = CException(self.__class__, 200, projectInfo['projectdirectory'])
            except:
                e = CException(self.__class__, 201, projectInfo['projectdirectory'])
        self.projectsListChanged.emit()
        self.backupDB()
        return e

    def openDirectoriesDef(self, mode, **kw):
        ccp4_dir = CCP4Utils.getHOME()
        if  not os.path.exists(ccp4_dir):
            raise CException(self.__class__, 108, ccp4_dir)
        if sys.platform.find("win32") >= 0:
            subd = 'CCP4'
            opsys = 'windows'
        else:
            subd = '.CCP4'
            opsys = 'unix'
        directories_def = os.path.join(ccp4_dir, subd, opsys, 'directories.def')
        #print 'CCP4ProjectsManager.openDirectoriesDef',directories_def
        if mode == 'w' and os.path.exists(directories_def):
            backup = directories_def + '.' + time.strftime('%y_%m_%d_%H-%M')
            shutil.copyfile(directories_def, backup)
        f = CCP4Utils.openFile(fileName=directories_def, mode=mode, overwrite = 1)
        if not f:
            raise CException(self.__class__, 110, directories_def)
        return [f, opsys]

    def loadDirectoriesDef(self, time=''):
        '''Load the ccp4i directories.def file into a dictionary'''
        rv = self.openDirectoriesDef('r', time='')
        f, opsys = rv
        text = f.read()
        f.close()
        lines = text.split('\n')
        array = {}
        directories = []
        projects = []
        self.header = []
        for line in lines:
            if len(line.strip()) <= 0:
                pass
            elif line[0] == '#':
                self.header.append(line)
            else:
                varname,line0 = line.split(' ', 1)
                line0 = line0.strip()
                try:
                    if line[0] != '"':
                        dtype, line0 = line0.split(' ', 1)
                except:
                    # Fails if no type in the definition
                    pass
                linetmp = line0.strip()
                if linetmp:
                    if linetmp[0] == '"' and linetmp[-1] == '"':
                        linetmp = linetmp.strip('"')
                    other = 1
                    for item in ['N_PROJECTS', 'PROJECT_ALIAS', 'PROJECT_PATH', 'N_DEF_DIRS', 'DEF_DIR_PATH', 'DEF_DIR_ALIAS']:
                        if varname.count(item):
                            other = 0
                else:
                    other = 1
                if other:
                    self.header.append(line)
                else:
                    array[varname] = (linetmp)
        if 'N_PROJECTS' in array:
            for n in range(1,int(array['N_PROJECTS']) + 1):
                if 'PROJECT_ALIAS,' + str(n) in array and 'PROJECT_PATH,' + str(n) in array:
                    path = array['PROJECT_PATH,' + str(n)]
                    if opsys == 'windows' and not re.search(r'/$', path):
                        path = path + '/'
                    projects.append([array['PROJECT_ALIAS,' + str(n)], path])
        if 'N_DEF_DIRS' in array:
            for n in range(1, int(array['N_DEF_DIRS']) + 1):
                if 'DEF_DIR_ALIAS,' + str(n) in array and 'DEF_DIR_PATH,' + str(n) in array:
                    path = array['DEF_DIR_PATH,' + str(n)]
                    if opsys == 'windows' and not re.search(r'/$', path):
                        path = path + '/'
                    directories.append([array['DEF_DIR_ALIAS,' + str(n)], path])
        if directories:
            self.directories = directories
        if projects:
            self.projects = projects
        #print 'CCP4ProjectsManager.loadDirectoriesDef',self.directories,self.projects
        return 0

    def getProjectDirectory(self, projectName=None, projectId=None, testAlias=False):
        if projectName in ['CCP4I2', 'CCP4I2_TOP'] or projectId in ['CCP4I2', 'CCP4I2_TOP'] :
            return CCP4Utils.getCCP4I2Dir()
        elif projectName == 'CCP4I2_TEST' or projectId == 'CCP4I2_TEST':
            return CCP4Utils.getTestTmpDir()
        try:
            pDir = self._db.getProjectDirectory(projectName=projectName, projectId=projectId)
        except:
            pDir = None
        if pDir is None and testAlias and projectName is not None:
            try:
                pDir = self._db.getAliasDirectory(projectName)
            except:
                pDir = None
        return pDir

    def projectNameStatus(self, projectName):
        proj_list = self._db.listProjectNames()
        if projectName in proj_list:
            return 1
        else:
            return 0

    def getProjectsList(self):
        proj_list = self._db.listProjectNames()
        #alias_list = self.directories.keys()
        alias_list = self._db.listDirectoryAliases()
        proj_list.sort()
        alias_list.sort()
        return [proj_list, alias_list]

    def aliasForDirectory(self, directory='', testAliases=True):
        abs_dir = os.path.abspath(directory)
        projectId, projectName, relpath = self._db.matchProjectDirectory(abs_dir, testAliases=testAliases)
        return projectName

    def interpretDirectory(self, directory='', ifAlias=False):
        '''Is the input directory one of the project or alias directories
        or in a sub-directory of one
        Return the project name and a relPath from project dir'''
        bestResult = [None, None, None]
        if self._db is None:
            return bestResult
        dirIn = os.path.abspath(directory)
        for  ProjectID, ProjectName, ProjectDirectory, ParentID in self._db.listProjects():
            try:
                relPath = os.path.relpath(dirIn, ProjectDirectory)
            except:
                pass
            else:
                if relPath == '.':
                    return ProjectName, '', ProjectID
                elif relPath[0:2] != '..':
                    if bestResult[0] is None or len(os.path.split(relPath)) < len(os.path.split(bestResult[1])):
                        bestResult = ProjectName, relPath, ProjectID
        if ifAlias:
            #for name,path in self.directories.items():
            for name, path in self._db.listDirectoryAliases():
                try:
                    relPath = os.path.relpath(dirIn, str(path))
                except:
                    pass
                else:
                    if relPath == '.':
                        return name, '', None
                    elif relPath[0:2] != '..':
                        if bestResult[0] is None or len(os.path.split(relPath)) < len(os.path.split(bestResult[1])):
                            bestResult = name, relPath, None
        return bestResult

    def add_quotes(self, alias):
        if alias.count(' '):
            return '"' + alias + '"'
        else:
            return alias

    def write_ccp4i_directories(self):
        text = ''
        for line in self.header:
            text = text + line +'\n'
        text = text + '%-30s_positiveint1      %s'%('N_PROJECTS',str(len(self.projects)))  + '\n'
                # '%-30s%s'%('PROJECT_ALIAS,0','\"\"') + '\n' + \
                # '%-30s%s'%('PROJECT_PATH,0','\"\"') + '\n'
        n = 0
        for alias,path in self.projects:
            n = n + 1
            text = text + '%-30s_text      %s'%('PROJECT_ALIAS,'+str(n),self.add_quotes(alias)) + '\n' + \
                        '%-30s_dir      %s'%('PROJECT_PATH,'+str(n),self.add_quotes(path)) + '\n'
        text = text + '%-30s_positiveint1      %s'%('N_DEF_DIRS',str(len(self.directories))) + '\n'
                # '%-30s_text     %s'%('DEF_DIR_ALIAS,0','\"\"') + '\n' + \
                # '%-30s_dir     %s'%('DEF_DIR_PATH,0','\"\"') + '\n'
        n = 0
        for alias,path in self.directories:
            n = n + 1
            text = text + '%-30s_text      %s'%('DEF_DIR_ALIAS,' + str(n), self.add_quotes(alias)) + '\n' + \
                        '%-30s_dir      %s'%('DEF_DIR_PATH,' + str(n), self.add_quotes(path)) + '\n'
        #print 'write_ccp4i_directories',text
        rv = self.open_directories_file(mode='w')
        if not rv[0]:
            f = rv[1]
            try:
                f.write(text)
                f.close()
            except:
                raise CException(self.__class__, 111)
        return 0

    def makeProjectDirectory(self, directory=None):
        '''Ensure that the project directory and CCP4_JOBS sub-directory exist
        If not directory was set then this will make CCP4_JOBS in current
        working directory'''
        if directory is None:
            raise CException(self.__class__, 113)
        if not os.path.exists(directory):
            try:
                os.mkdir(directory)
            except:
                raise CException(self.__class__, 114, 'Directory: ' + directory)
        for sub_dir in ['CCP4_JOBS', 'CCP4_TMP', 'CCP4_IMPORTED_FILES', 'CCP4_PROJECT_FILES', 'CCP4_COOT']:
            sub_dir_path = os.path.join(directory, sub_dir)
            if not os.path.exists(sub_dir_path):
                try:
                    os.mkdir(sub_dir_path)
                except:
                    raise CException(self.__class__, 115, 'Directory: ' + directory)

    def save(self, fileName=None):
        errorReport = CErrorReport()
        if fileName is None:
            fileName = os.path.join(getDotDirectory(),'projects.xml')
        from .CCP4File import CI2XmlDataFile
        f = CI2XmlDataFile(fullPath=fileName)
        f.header.function.set('PROJECTDIRECTORIES')
        f.header.setCurrent()
        bodyEtree = self.getEtree()
        try:
            f.saveFile(bodyEtree=bodyEtree)
        except CException as e:
            errorReport.append(e)
        except:
            errorReport.append(self.__class__, 130, fileName)
        return errorReport

    def restore(self, fileName=None):
        errorReport = CErrorReport()
        if fileName is None:
            fileName = os.path.join(getDotDirectory(), 'projects.xml')
        from .CCP4File import CI2XmlDataFile
        # Beware CDataFile will attempt to interpret any file path and call PROJECTSMAMANGER
        # --> recursion and crash
        f = CI2XmlDataFile(baseName='projects.xml', relPath=getDotDirectory())
        #print 'PM.restore before loadFile'
        try:
            f.loadFile()
        except CException as e:
            errorReport.extend(e)
        except Exception as e:
            errorReport.append(self.__class__, 117, str(fileName))
        else:
            #print 'PM.restore done loadFile'
            if f.header.function != 'PROJECTDIRECTORIES':
                raise CException(self.__class__, 107, 'Filename: '  +fileName)
            e = self.setEtree(f.getBodyEtree())
            if len(e) > 0:
                errorReport.extend(e)
        if len(errorReport) > 0:
            print('PROJECTSMANAGER restore errrors:\n', errorReport.report())
        return errorReport

    def newJob(self, taskName=None, jobTitle=None, projectId=None, jobNumber=None):
        #if jobTitle is None:
        #  jobTitle = CCP4TaskManager.TASKMANAGER().getTitle(taskName)
        jobId = self._db.createJob(projectId, taskName, jobTitle=jobTitle, jobNumber=jobNumber)
        jobInfo = self._db.getJobInfo(jobId, ['jobnumber', 'projectid'])
        projectName = self._db.getProjectInfo(jobInfo['projectid'], 'projectName')
        #print 'PROJECTMANAGER.newJob', projectId,projectName,jobInfo['jobnumber']
        #self.emitSignal('jobCreated',jobId)
        #self.backupDB()
        return jobId, projectName, jobInfo['jobnumber']

    def deleteJob(self, jobId=None, jobNumber=None, projectName=None, projectId=None, importFiles=[], deleteImportFiles=True):
        from ..dbapi import CCP4DbApi
        from ..dbapi import CCP4DbUtils
        print('PROJECTSMANAGER.deleteJob importFiles', jobId, importFiles, deleteImportFiles)
        if jobId is None:
            jobId = self._db.getJobId(jobNumber=jobNumber, projectName=projectName, projectId=projectId)
            if jobId is None:
                raise CException(self.__class__, 140)
        # Delete files from CCP4_IMPORT_FILES
        # Beware may need to query database so do this before zapping job from db
        if deleteImportFiles and len(importFiles) > 0:
            importDir = os.path.join(self.getProjectDirectory(projectId=projectId), 'CCP4_IMPORTED_FILES')
            for  importId,filename in importFiles:
                filePath = os.path.join(importDir, filename)
                #print 'PROJECTSMANAGER.deleteJob attempting to delete',filePath
                if os.path.exists(filePath):
                    try:
                        os.remove(filePath)
                    except:
                        pass
                else:
                    importInfo = self._db.getImportFileInfo(importId=importId, mode=['importnumber', 'sourcefilename'])
                    #print 'PROJECTSMANAGER.deleteJob importInfo',importInfo
                    sourceFileName, ext = os.path.splitext(os.path.split(importInfo['sourcefilename'])[-1])
                    filePath = os.path.join(importDir,sourceFileName + '_' + str(importInfo['importnumber']) + ext)
                    #print 'PROJECTSMANAGER.deleteJob attempting to delete',filePath
                    if os.path.exists(filePath):
                        try:
                            os.remove(filePath)
                        except:
                            pass
        try:
            status = self._db.getJobInfo(jobId, 'status')
        except:
            # Assume job no longer exisits - this could happend since job may appear in the tree of
            # dependent jobs (from CDbApi.getFollowOnJobs() as used in project viewer code) more than once
            print('PROJECTMANAGER.deleteJob job does not exist', jobId)
            return
        if status in (CCP4DbApi.JOB_STATUS_RUNNING, CCP4DbApi.JOB_STATUS_REMOTE):
            err = JOBCONTROLLER().killJobProcess(jobId=jobId)
            if len(err) > 0:
                raise err
        jobDir = self.jobDirectory(jobId=jobId)
        if len(importFiles) == 0 or deleteImportFiles:
            try:
                shutil.rmtree(jobDir)
            except:
                raise CException(self.__class__, 141, str(jobDir))
            self._db.deleteServerJob(jobId=jobId)
            self._db.deleteJob(jobId=jobId, deleteChildren=True)
        else:
            try:
            #if 1:
                copyJobDir = jobDir + '_copy'
                shutil.move(jobDir, copyJobDir)
                os.mkdir(jobDir)
                for importId, filename in importFiles:
                    src =os.path.normpath(os.path.join(copyJobDir, filename))
                    if os.path.exists(src):
                        dst = os.path.normpath(os.path.join(jobDir, filename))
                        shutil.move(src, dst)
                # Find the params or input_params file and save
                # This should help if need to recover the database from project directory
                # Dont change pluginName - can not find right def file to enable loading
                for name in ['input_params.xml', 'params.xml']:
                    if os.path.exists(os.path.join(copyJobDir, name)):
                        shutil.move(os.path.join(copyJobDir, name), os.path.join(jobDir, name))
                        break
                try:
                    shutil.rmtree(copyJobDir)
                except:
                    pass
            except:
                raise CException(self.__class__, 143, str(jobDir))
            # Add a dummy program.xml
            dummy = etree.Element('dummy')
            CCP4Utils.saveEtreeToFile(dummy, os.path.join(jobDir, 'program.xml'))
            self._db.setJobToImport(jobId=jobId, projectId=projectId)
            # Save the job db backup
            backup = CCP4DbUtils.CJobDbBackup(jobId=jobId, jobDirectory=jobDir, projectId=projectId)
            backup.save(updateAll=True)
        #self.backupDB()

    def runInternalTask(self, jobId, projectId, taskName):
        print('runInternalTask')
        from . import CCP4PluginScript
        from . import CCP4TaskManager
        from ..dbapi import CCP4DbApi
        cls = CCP4TaskManager.TASKMANAGER().getPluginScriptClass(taskName)
        if cls is not None:
            db = self.db()
            db.updateJobStatus(jobId=jobId, status='Running')
            jobInfo = db.getJobInfo(jobId=jobId, mode=['projectid', 'projectname', 'jobnumber'])
            workDirectory = os.path.join(db.getProjectDirectory(projectId=jobInfo['projectid']), 'CCP4_JOBS', 'job_' + str(jobInfo['jobnumber']))
            if not os.path.exists(workDirectory):
                os.mkdir(workDirectory)
            name = str(jobInfo['projectname']) + '_' + str(jobInfo['jobnumber'])
            pluginObj = cls(parent=self, name=name, jobId=jobId, projectId=projectId)
            pluginObj.container.header.projectId = projectId
            pluginObj.container.header.jobId = jobId
            comFile = self.makeFileName(jobId=jobId, mode='JOB_INPUT')
            pluginObj.container.loadDataFromXml(comFile)
            dbHandler = CCP4PluginScript.CDatabaseHandler(projectId=jobInfo['projectid'], jobNumber=jobInfo['jobnumber'], projectName=jobInfo['projectname'])
            dbHandler.openDb()
            pluginObj.setDbData(handler=dbHandler, projectId= jobInfo['projectid'], jobNumber=jobInfo['jobnumber'], projectName=jobInfo['projectname'], jobId=jobId)
            ret = pluginObj.process()
            if ret == CCP4PluginScript.CPluginScript.SUCCEEDED:
                self.db().gleanJobFiles(jobId=jobId, container=pluginObj.container, projectId=projectId, roleList=[CCP4DbApi.FILE_ROLE_OUT])
                self.updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_FINISHED)
            else:
                self.updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_FAILED)
            return True
        else:
            return False

    def cleanupJob(self,jobDirectory):
        pass

    @QtCore.Slot()
    def _checkForFinishedJobs(self):
        self.doCheckForFinishedJobs.emit()

    @QtCore.Slot()
    def checkForFinishedJobs(self):
        t = time.time()
        self.blockExit = True
        finishedJobs = self._db.getRecentlyFinishedJobs(after=self.lastJobFinishCheckTime)
        self.lastJobFinishCheckTime = t
        t = time.time()
        startedJobs = self._db.getRecentlyStartedJobs(after=self.lastJobStartedCheckTime)
        self.lastJobStartedCheckTime = t
        self.blockExit = False

    @QtCore.Slot(dict)
    def handleJobToDelete(self, args={}):
        from ..dbapi import CCP4DbApi
        if not bool(PREFERENCES().DELETE_INTERACTIVE_JOBS):
            self.db().updateJobStatus(jobId=args.get('jobId'), status=CCP4DbApi.JOB_STATUS_FINISHED)
            return
        if bool(PREFERENCES().SHOW_DELETE_INTERACTIVE_JOBS):
            jobInfo = self.db().getJobInfo(args.get('jobId'), ['taskname', 'projectid', 'jobnumber'])
            from ..qtgui import CCP4ProjectViewer
            pV = CCP4ProjectViewer.PROJECTVIEWER(jobInfo['projectid'])
            print('handleJobToDelete', pV, jobInfo)
            if pV is not None:
                pV.queryDeleteJob(args.get('jobId'), jobInfo)
                return
        # Silently delete job
        try:
            self.db().deleteJob(jobId=args.get('jobId'), deleteChildren=True)  
        except:
            print('PROJECTSMANAGER.handleJobToDelete failed to delete job')

    def makeFileName(self, jobId=None, mode='PARAMS', jobInfo=None, qualifier=None):
        '''Generate suitable name for job output file
        This should give same result as CPluginScript.makeFileName()'''
        myDir = self.jobDirectory(jobId=jobId)
        taskName = self._db.getJobInfo(jobId=jobId, mode='taskname')
        #exceptions = { 'molrep_mr' : { 'LOG' : 'foo.bar' } }
        defNames = {'ROOT' : '', 'PARAMS' : 'params.xml', 'JOB_INPUT' : 'input_params.xml', 'PROGRAMXML' : 'program.xml',
                    'LOG' : 'log.txt', 'STDOUT' : 'stdout.txt', 'STDERR' : 'stderr.txt',
                    'INTERRUPT' : 'interrupt_status.xml', 'DIAGNOSTIC' : 'diagnostic.xml', 'REPORT' : 'report.html',
                    'DIAGNOSTIC_REPORT' : 'diagnostic_report.html', 'TABLE_RTF' : 'tables.rtf',
                    'TABLES_DIR' : 'tables_as_csv_files', 'XML_TABLES_DIR' : 'tables_as_xml_files',
                    'LOG' : 'log.txt', 'COM' : 'com.txt', 'MGPICDEF' : 'report.mgpic.py', 'PIC' : 'report.png', 'RVAPIXML' : 'i2.xml'}
        #fileName = exceptions.get(taskName,{}).get(mode,None)
        #if fileName is None:
        fileName = defNames.get(mode, 'unknown.unk')
        if qualifier is not None:
            base, ext = fileName.split('.', 1)
            fileName = base + '_' + str(qualifier) + '.' + ext
        return os.path.join(myDir, fileName)

    def jobDirectory(self,jobId=None, jobNumber=None, projectDirectory=None, projectId=None, create=True, subDir=None):
        try:
            directory = self._db.jobDirectory(jobId=jobId, jobNumber=jobNumber, projectId=projectId, projectDirectory=projectDirectory)
        except:
            excptTxt = 'jobId:' + str(jobId) + ' jobNumber:' + str(jobNumber) + ' projectDirectory:' \
                        + str(projectDirectory) + ' projectId:' + str(projectId)
            raise CException(self.__class__, 144, excptTxt)
        if not os.path.exists(directory):
            if create:
                try:
                    os.mkdir(directory)
                except:
                    raise CException(self.__class__, 142, directory)
            else:
                return None
        if subDir is not None:
            directory = os.path.join(directory, subDir)
            if not os.path.exists(directory):
                try:
                    os.mkdir(directory)
                except:
                    raise CException(self.__class__, 142, directory)
        return directory

    def updateJobStatus(self, jobId=None, status=None, runMode=None):
        self._db.updateJobStatus(jobId=jobId, status=status)
        self._db.commit()
        self.jobUpdated.emit(jobId)

    def setInputFileNames(self, container, contextJobId=None, projectId=None, force=False):
        from ..dbapi import CCP4DbApi
        myErrorReport = CErrorReport()
        try:
            dataList = container.inputData.dataOrder()
        except:
            myErrorReport.append(self.__class__, 153)
            return myErrorReport
        #print 'setInputFileNames',contextJobId,dataList
        for objectName in dataList:
            #try:
                dobj = container.inputData.find(objectName)
                #print 'PROJECTSMANAGER.setInputFileNames',objectName,dobj,dobj.isSet()
                if isinstance(dobj,CCP4File.CDataFile) and (force or not dobj.isSet()):
                    if contextJobId < 0:
                        #dobj.setDbFileId(None)
                        dobj.unSet()
                    else:
                        if dobj.qualifiers('fromPreviousJob'):
                            mimeType = dobj.qualifiers('mimeTypeName')
                            if CCP4DbApi.FILETYPES_TEXT.count(mimeType):
                                fileType = CCP4DbApi.FILETYPES_TEXT.index(mimeType)
                            else:
                                fileType = CCP4DbApi.FILETYPES_TEXT[0]
                            subType =  dobj.qualifiers('requiredSubType')
                            contentFlag = dobj.qualifiers('requiredContentFlag')
                            fileIdList = self.db().getFileByJobContext(contextJobId=contextJobId, fileType=fileType,
                                                                       subType=subType, contentFlag=contentFlag, projectId=projectId)
                            #print 'PROJECTSMANAGER.setInputFileNames',objectName,dobj,mimeType,fileType,fileIdList
                            if len(fileIdList) > 0:
                                fileInfo = self.db().getFileInfo(fileId=fileIdList[0], mode=['jobid', 'filename', 'relpath', 'projectid', 'annotation', 'filecontent', 'filesubtype'])
                                if projectId is None:
                                    projectId = fileInfo['projectid']
                                #print 'PROJECTSMANAGER.setInputFileNames',objectName,contextJobId,fileInfo
                                #dobj.setDbFileId(fileIdList[0])
                                dobj.set({'baseName' : fileInfo['filename'], 'relPath' : fileInfo['relpath'], 'project' : projectId,
                                          'annotation' : fileInfo['annotation'], 'dbFileId' : fileIdList[0],
                                          'contentFlag' : fileInfo['filecontent'], 'subType' : fileInfo['filesubtype']})
                            else:
                                #dobj.setDbFileId(None)
                                dobj.unSet()
            #except:
            #    myErrorReport.append(self.__class__,154,objectName)
            #print 'setInputData',objectName,dobj
        return myErrorReport

    def setOutputFileNames(self, container=None, projectId=None, jobNumber=None, force=True):
        '''Where output data is a file name set suitable name'''
        #print 'PROJECTSMANAGER.setOutputFileNames',container,projectId,jobNumber
        #pName = self.db().getProjectInfo(projectId=projectId,mode='projectname')
        relPath = 'CCP4_JOBS'
        for job in jobNumber.split('.'):
            relPath = os.path.join(relPath,'job_' + job)
        #relPath = relPath[0:-1]
        #jobName = taskName+'_'+re.sub('\.','_',jobNumber)
        jobName = ''
        myErrorReport = CErrorReport()
        if container is None:
            myErrorReport.append(self.__class__, 151)
        else:
            try:
                dataList = container.outputData.dataOrder()
            except:
                myErrorReport.append(self.__class__, 151)
        if len(myErrorReport) > 0:
            return myErrorReport
        for objectName in dataList:
            try:
                dobj = container.outputData.find(objectName)
                #print 'setOutputData get',objectName,dobj.get(),dobj.isSet()
                if isinstance(dobj,CCP4File.CDataFile) and (force or not dobj.isSet()):
                    dobj.setOutputPath(jobName=jobName, projectId=projectId, relPath=relPath)
            except:
                myErrorReport.append(self.__class__, 152, objectName)
        return myErrorReport

    def getJobParams(self, jobId=None, inputParams=False):
        from . import CCP4Container
        from . import CCP4TaskManager
        if inputParams:
            filename = self.makeFileName(jobId=jobId, mode='JOB_INPUT')
        else:
            filename = self.makeFileName(jobId=jobId, mode='PARAMS')
        if not os.path.exists(filename):
            return None
        header = CCP4File.xmlFileHeader(filename)
        defFile = CCP4TaskManager.TASKMANAGER().lookupDefFile(str(header.pluginName), str(header.pluginVersion))
        if defFile is None:
            return None
        container = CCP4Container.CContainer(definitionFile=defFile)
        container.loadDataFromXml(filename)
        return container

    def importFilePath(self, jobId=None, projectId=None, baseName=None, fileExtensions=None,
                       use_CCP4_IMPORTED_FILES=True, importNumber=1, newFile=True, importId=None):
        # Get the file baseName from the sourcefilename saved in the import table
        if baseName is None and importId is not None:
            importInfo = self.db().getImportFileInfo(importId=importId, mode=['sourcefilename', 'importnumber'])
            baseName = os.path.split(importInfo['sourcefilename'])[1]
            base,ext = os.path.splitext(baseName)
            if not ext[1:] in fileExtensions:
                baseName = base + '.' + fileExtensions[0]
            importNumber = importInfo['importnumber']
        # Convert file extension to standard
        if fileExtensions is not None and len(fileExtensions) > 0:
            base, ext = os.path.splitext(baseName)
            if not ext[1:] in fileExtensions:
                baseName = base + '.' + fileExtensions[0]
        if not use_CCP4_IMPORTED_FILES:
            jobDir = self.db().jobDirectory(jobId=jobId)
            filePath = os.path.join(jobDir,baseName)
        else:
            filePath = os.path.join(self.db().getProjectDirectory(projectId=projectId), 'CCP4_IMPORTED_FILES', baseName)
            root,ext = os.path.splitext(filePath)
            root = root + '_'
            filePath = root + str(importNumber) + ext
            if not newFile:
                return filePath, importNumber
            else:
                while os.path.exists(filePath):
                    importNumber += 1
                    filePath = root + str(importNumber) + ext
        #print 'PROJECTSMANAGER.importFilePath filePath',jobId,'projectId',projectId,use_CCP4_IMPORTED_FILES,filePath
        return filePath, importNumber

    def alreadyImportedId(self, sourceFileName=None, projectId=None, contentFlag=None, sourceFileReference=None, fileType=None):
        # is there already an imported file with same sourceFileName and same checksum?
        # Return the matching importId OR the checksum for sourceFileName neded to create new import record
        importList=self.db().getImportedFile(sourceFileName=sourceFileName, projectId=projectId, fileContent=contentFlag,
                                             reference=sourceFileReference, fileType=fileType)
        if len(importList) == 0:
            importId = None
            checksum = None
            dbFileId = None
            annotation = None
        else:
            checksum = importList[0][2]
            dbFileId = importList[0][1]
            importId = importList[0][0]
            annotation = importList[0][3]
        sourceChecksum = None
        if dbFileId is None or checksum is not None:
            try:
                fileObj = CCP4File.CDataFile(fullPath=sourceFileName)
                sourceChecksum = fileObj.checksum()
                # If checksums dont match cancel the importId
                if checksum is not None and sourceChecksum != checksum:
                    dbFileId = None
                    importId = None
            except:
                pass
        #print 'CProjectmaanger.alreadyImportedId',sourceFileName,dbFileId,checksum,sourceChecksum
        return dbFileId, importId, sourceChecksum, annotation

    def importFiles(self, jobId=None, container=None, projectId=None):
        #print 'PROJECTSMANAGER.importFiles',jobId,repr(container),projectId
        from ..dbapi import CCP4DbApi
        err= CErrorReport()
        if container is None:
            try:
                container = self.db().getParamsContainer(jobId=jobId)
            except CException as e:
                err.extend(e)
            except Exception as e:
                err.append(self.__class__, 161, 'Error loading container for job')
        ifImportedFiles = False
        projectId = self.db().getJobInfo(jobId=jobId, mode='projectId')
        for key in container.inputData.dataOrder():
            obj0 = container.inputData.__getattr__(key)
            objList, xmlText, keyValues = obj0.saveToDb()
            #print 'PROJECTSMANAGER.importFiles',key,obj0,objList,xmlText
            for obj in objList:
                if not isinstance(obj, CCP4File.CDataFile):
                    #print 'PROJECTSMANAGER.importFiles not file',obj,obj.objectName()
                    pass
                else:
                    #print 'PROJECTSMANAGER.importFiles COPY',key,container.guiParameters.__getattr__('COPY_'+key,True)
                    if obj.exists() and not obj.dbFileId.isSet() and obj.qualifiers('saveToDb') and not obj.qualifiers('isDirectory') \
                                and (container.__getattr__('guiParameters', None) is None \
                                     or container.guiParameters.__getattr__('COPY_' + key, True)):
                        obj.blockSignals(True)
                        #print 'PROJECTSMANAGER.importFiles',obj.objectName(),obj.__dict__.get('sourceFileName','NO sourceFileName'),obj
                        if 'sourceFileName' in obj.__dict__ and obj.__dict__['sourceFileName'] is not None:
                            # The file has already being imported (i.e. for CMiniMtzDataFile)
                            sourceFileName = obj.__dict__['sourceFileName']
                            #print 'PROJECTSMANAGER.importFiles previously imported from sourceFileName',sourceFileName
                            sourceFileReference = obj.__dict__.get('sourceFileReference', None)
                            sourceFileAnnotation = obj.__dict__.get('sourceFileAnnotation', None)
                            targetObj = obj
                            targetFile, importNumber=self.importFilePath(jobId=jobId, projectId=projectId, baseName=os.path.split(sourceFileName)[-1],
                                                                         fileExtensions=obj.qualifiers('fileExtensions'))
                            #Has the monster-MTZ been imported before for a different set of columns
                            dbFileId, importId, checksum, dbAnnotation = self.alreadyImportedId(sourceFileName=sourceFileName, projectId=projectId)
                            if importId is not None:
                                #it has been imported before so just use the same 'targetFile'
                                targetFile, importNumber = self.importFilePath(projectId=projectId, importId=importId,
                                                                               fileExtensions=obj.qualifiers('fileExtensions'), newFile=False)
                                #print 'targetFile from importId',targetFile,importNumber
                            else:
                                targetFile, importNumber=self.importFilePath(jobId=jobId, projectId=projectId, baseName=os.path.split(sourceFileName)[-1],
                                                                             fileExtensions=obj.qualifiers('fileExtensions'))
                        else:
                            # Test if user attempting to import an already imported file
                            if ('sourceFileName' in obj.__dict__) == 'COPYDONE':
                                sourceFileName = None
                            else:
                                sourceFileName = obj.__str__()
                            sourceFileReference = obj.__dict__.get('sourceFileReference', None)
                            sourceFileAnnotation = obj.__dict__.get('sourceFileAnnotation', None)
                            mimeTypeName = obj.qualifiers('mimeTypeName')
                            if CCP4DbApi.FILETYPES_TEXT.count(mimeTypeName):
                                sourceFileType = CCP4DbApi.FILETYPES_TEXT.index(mimeTypeName)
                            else:
                                sourceFileType = None
                            dbFileId, importId, checksum, dbAnnotation = self.alreadyImportedId(sourceFileName=sourceFileName,projectId=projectId,
                                                                                                sourceFileReference=sourceFileReference,fileType=sourceFileType)
                            #print 'File has been previously imported?',sourceFileName,'previous file id',dbFileId,'checksum comparison',checksum
                            if dbFileId is None:
                                targetFile,importNumber=self.importFilePath(jobId=jobId, projectId=projectId, baseName=os.path.split(sourceFileName)[-1],
                                                                            fileExtensions=obj.qualifiers('fileExtensions'))
                                value = {}
                                if obj.contentFlag.isSet():
                                    value['contentFlag'] = int(obj.contentFlag)
                                if obj.subType.isSet():
                                    value['subType'] = int(obj.subType)
                                targetObj = obj.__class__(value=value, name=obj.objectName())
                                if obj.__dict__.get('sourceFileAnnotation', None) is not None:
                                    targetObj.__dict__['sourceFileAnnotation'] = copy.deepcopy(obj.__dict__['sourceFileAnnotation'])
                            else:
                                # This external file is already in the database - update the obj with the appropriate data
                                # and skip copying or creating an import file record
                                # setDbFileId should force a dataChanged signal and an update to the widget
                                # which will not then overwrite the data object with old values
                                obj.setDbFileId(dbFileId)
                                print('PROJECTSMANAGER.importFiles handling previous import', obj.get(), mimeTypeName, sourceFileType)
                                sourceFileName = None
                                sourceFileReference = None
                                targetObj = None
                        #print 'PROJECTSMANAGER.importFiles',sourceFileName,sourceFileReference,repr(targetObj)
                        # Copy a sourceFileName to the CCP4_IMPORT_FILES directory
                        if sourceFileName is not None and not os.path.exists(targetFile):
                            print('Copying file', sourceFileName, 'to', targetFile)
                            try:
                                if obj.qualifiers('isDirectory'):
                                    shutil.copytree(sourceFileName, targetFile)
                                else:
                                    shutil.copyfile(sourceFileName, targetFile)
                                # Set the targetObj path now that the file exisits to get around validity failures
                                if not targetObj.isSet():
                                    targetObj.setFullPath(targetFile)
                            except:
                                err.append(self.__class__, 160, sourceFileName + ' to ' + targetFile)
                                targetObj = None
                        if targetObj is not None:
                            # Set annotation on the imported targetObj and create an imported file record in the db
                            # Beware downloaded file already has a more use annotatio set!
                            if not targetObj.annotation.isSet():
                                if sourceFileAnnotation is not None:
                                    targetObj.annotation = sourceFileAnnotation
                                else:
                                    anno = obj.qualifiers('guiLabel')
                                    if anno is NotImplemented or anno is None:
                                        anno = 'Imported'
                                    else:
                                        anno = anno + ' imported'
                                    targetObj.annotation = anno + ' from ' + os.path.split(sourceFileName)[1] + ' by job ' + \
                                                                str(self.db().getJobInfo(jobId=jobId, mode='jobnumber'))
                            #print 'PROJECTSMANAGER.importFiles targetObj sourceFileAnnotation',targetObj,targetObj.__dict__.get('sourceFileAnnotation',None)
                            try:
                                fileId = self.db().createFile(fileObject=targetObj, projectId=projectId, jobId=jobId, sourceFileName=sourceFileName,
                                                              reference=sourceFileReference, importNumber=importNumber, sourceFileAnnotation=sourceFileAnnotation)
                                print('Database recording file', targetObj.fullPath.__str__(), 'imported from', sourceFileName, 'file id:', fileId)
                            except CException as e:
                                err.extend(e)
                            except Exception as e:
                                #print 'PROJECTSMANAGER.importFiles',str(e)
                                err.append(self.__class__, 161, sourceFileName)
                            else:
                                print('PROJECTSMANAGER.importFiles resetting path', obj.objectName(), targetFile)
                                ifImportedFiles = True
                                if targetObj != obj:
                                    obj.setFullPath(targetFile)
                                    obj.annotation = targetObj.annotation.__str__()
                                obj.dbFileId = fileId
                        obj.blockSignals(False)
        return ifImportedFiles, err

    def getSceneFiles(self, jobId=None):
        scenefileList = glob.glob(os.path.join(self.jobDirectory(jobId=jobId), 'scene_*.scene.xml'))
        if len(scenefileList) == 0:
            scenefileList = glob.glob(os.path.join(self.jobDirectory(jobId=jobId), '*.scene.xml'))
            if len(scenefileList) == 0:
                scenefileList = glob.glob(os.path.join(self.jobDirectory(jobId=jobId), '*', '*.scene.xml'))
        scenefileList.sort()
        #print 'PROJECTSMANGER.getSceneFiles', scenefileList
        return scenefileList

    def cleanupAllProjects(self, context='temporary'):
        projectList = self.db().getRecentProjects()
        print('cleanupAllProjects',projectList)
        for pid, pname, access in projectList:
            cleanup = CPurgeProject(projectId=pid,db=self.db())
            cleanup.purgeProject(context=context)

    def makeDbXml(self, projectId=None, after=None):
        projectInfo = self.db().getProjectInfo(projectId=projectId)
        dbxml = os.path.join( projectInfo['projectdirectory'], 'CCP4_TMP', 'DATABASE' + str(int(time.time())) + '.db.xml')
        print('Creating temporary database xml file in:', dbxml)
        errReport = self.db().exportProjectXml(projectId, fileName=dbxml, after=after)
        return errReport

    def importProject(self, compressedFile):
        dbfile = self.extractDatabaseXml(compressedFile)
        projectInfo = self.getProjectInfoFromXml(dbfile)
        #print 'importProject projectInfo',projectInfo

    def projectInfo(self,compressedFile):
        dbfile = self.extractDatabaseXml(compressedFile)
        projectInfo = self.getProjectInfoFromXml(dbfile)
        #print 'importProject projectInfo',projectInfo

    def extractDatabaseXml(self, compressedFile, tempDir=None):
        if compressedFile.endswith('tar.gz'):
            return self.extractDatabaseXmlFromTar(compressedFile, tempDir=tempDir)
        elif compressedFile.endswith('zip'):
            return self.extractDatabaseXmlFromZip(compressedFile, tempDir=tempDir)
        else:
            raise CErrorReport(self.__class__, 178, compressedFile)

    def extractDatabaseXmlFromTar(self, compressedFile, tempDir=None):
        if tempDir is None:
            importDir = tempfile.gettempdir()
        else:
            importDir = tempDir
        try:
            tf = tarfile.open(compressedFile, mode='r:gz')
        except:
            raise CException(self.__class__, 176, compressedFile)
        try:
            member = tf.getmember('DATABASE.db.xml')
            tf.extractall(path=importDir, members=[member])
        except:
            tf.close()
            raise CException(self.__class__, 177, 'Extracting from ' + compressedFile + ' to ' + importDir)
        tf.close()
        dbfile = os.path.join(importDir, 'DATABASE.db.xml')
        #print 'extractDatabaseXml dbfile', dbfile
        return dbfile

    def extractDatabaseXmlFromZip(self, compressedFile, tempDir=None, diagnostic=False):
        from ..qtcore import CCP4Export
        if tempDir is None:
            importDir = tempfile.gettempdir()
        else:
            importDir = tempDir
        try:
            zip = zipfile.ZipFile(compressedFile, mode='r', allowZip64=CCP4Export.ALLOWZIP64)
        except Exception as e:
            err = CException(self.__class__, 176, compressedFile + '\n', str(e))
            if diagnostic:
                print(err.report())
            raise err
        if diagnostic:
            for zinfo in zip.infolist():
                print('zinfo', zinfo.filename, zinfo.date_time)
        try:
            zip.extractall(path=importDir, members=['DATABASE.db.xml'])
        except Exception as e:
            zip.close()
            raise CException(self.__class__, 177, 'Extracting from ' + compressedFile + ' to ' + importDir + '\n' + str(e))
        zip.close()
        dbfile = os.path.join(importDir, 'DATABASE.db.xml')
        #print 'extractDatabaseXml dbfile', dbfile
        return dbfile

    def getJobInputFiles(self, projectDir='', jobIdList=[], jobNumberList=[], useDb=False, excludeI2files=False):
        # Find the input files that would also need to be exported with a job
        # ***** Should have checked jobs have valid input before we get this far ************
        from . import CCP4Container
        from . import CCP4Utils
        errReport = CErrorReport()
        fileList = []
        fileIdList = []
        fromJobIdList = []
        if useDb:
            demoDataDirList = DEMODATAMANAGER().getTestDatasets()
            demoData = []
            for f1, f2 in demoDataDirList:
                demoData.append(f1)
            demoData = tuple(demoData)
            #print 'getJobInputFiles demoDataDirs', demoData, type(demoData), excludeI2files
            dbList = self.db().getFilesUsedInJobList(jobList=jobIdList)
            #print 'PROJECTSMANGER.getJobInputFiles dbList', jobIdList, dbList
            for fileId,baseName,jobId,jobNumber,importId,sourceFilename in dbList:
                if excludeI2files:
                    print('PROJECTSMANGER.getJobInputFiles', sourceFilename, os.path.relpath(sourceFilename, CCP4Utils.getCCP4I2Dir()), os.path.relpath(sourceFilename, projectDir))
                    if os.path.relpath(sourceFilename,CCP4Utils.getCCP4I2Dir()).startswith(('wrappers', 'pipelines')):
                        print('PROJECTSMANGER.getJobInputFiles sourceFilename exclude from wrapper', sourceFilename)
                    elif  os.path.relpath(sourceFilename,projectDir).startswith(demoData):
                        print('PROJECTSMANGER.getJobInputFiles sourceFilename exclude from demoData', sourceFilename)
                # Beware the mini-MTZ converted from imported file is in the job directory
                if importId is not None and  os.path.exists(os.path.join(projectDir, 'CCP4_IMPORTED_FILES', baseName)):
                    fileList.append( ['CCP4_IMPORTED_FILES', baseName])
                    fileIdList.append(fileId)
                    if not jobId in jobIdList and not jobId in fromJobIdList:
                        fromJobIdList.append([jobId, jobNumber])
                else:
                    relPath =  os.path.join('CCP4_JOBS')
                    for jN in jobNumber.split('.'):
                        relPath = os.path.join(relPath,'job_' + jN)
                    if os.path.exists(os.path.join(projectDir, relPath, baseName)):
                        fileList.append([relPath, baseName])
                        fileIdList.append(fileId)
                        if not jobId in jobIdList and not jobId in fromJobIdList:
                            fromJobIdList.append([jobId, jobNumber])
                    else:
                        errReport.append(self.__class__, 182, 'For job number:' + str(jobNumber))
        else:
            for jobNumber in jobNumberList:
                try:
                    directory = os.path.join(projectDir, 'CCP4_JOBS')
                    for jN in jobNumber.split('.'):
                        directory = os.path.join(directory, 'job_' + jN)
                    container = CCP4Container.CContainer(parent=self)
                    container.loadDataFromXml(os.path.join(directory, 'input_params.xml'))
                except:
                    errReport.append(self.__class__, 180, 'For job number:' + str(jobNumber))
                else:
                    #print 'PROJECTSMANGER.getJobInputFiles dataOrder',container.inputData.dataOrder()
                    for key in container.inputData.dataOrder():
                        obj0 = container.inputData.__getattr__(key)
                        try:
                            objList, xmlText, keyValues = obj0.saveToDb()
                            for obj in objList:
                                if isinstance(obj,CCP4File.CDataFile) and obj.isSet() and obj.exists():
                                    item = [obj.relPath.__str__(), obj.baseName.__str__()]
                                    #print 'getJobInputFiles',key,item
                                    if not item in fileList:
                                        fileList.append(item)
                                        fileIdList.append(obj.dbFileId.__str__() )
                                        fromJobId = self.db().getFileInfo(fileId=obj.dbFileId.__str__(), mode=['jobid', 'jobnumber'])
                                        #print 'PROJECTSMANGER.getJobInputFiles',key,obj.dbFileId.__str__(),fromJobId
                                        if not fromJobId in fromJobIdList:
                                            fromJobIdList.append(fromJobId)
                        except:
                            print('getJobInputFiles error', key)
                            errReport.append(self.__class__, 181, 'For job number:' + str(jobNumber))
        return fileList, fileIdList, fromJobIdList, errReport

    def compressProject(self, projectId, after=None, jobList=None, excludeI2files=False, fileName=None, blocking=False, parentWidget=None):
        from ..qtcore import CCP4Export
        #print 'CProjectManager.compressProject',after,jobList,excludeI2files,fileName
        projectInfo = self.db().getProjectInfo(projectId=projectId)
        # If there is a limited set of jobs then find the input jobs that are not output by jobs on that list
        inputFilesList, inputFileIdList, fromJobList, errReport = self.getJobInputFiles(projectDir=projectInfo['projectdirectory'], jobIdList=jobList, useDb=True, excludeI2files=excludeI2files)
        #print 'CProjectManager.compressProject inputFilesList',inputFilesList
        #print 'CProjectManager.compressProject fromJobList',fromJobList
        fromJobIdList = []
        fromJobNumberList = []
        for item in fromJobList:
            fromJobIdList.append(item[0])
            fromJobNumberList.append(item[1])
        dbxml = os.path.join(projectInfo['projectdirectory'], 'CCP4_TMP', 'DATABASE' + str(int(time.time())) + '.db.xml')
        print('Creating temporary database xml file in:',dbxml)
        # exportProjectXml returns list of TOP-LEVEL jobNumbers for the export
        jobNumberList, errReport = self.db().exportProjectXml(projectId, fileName=dbxml, recordExport=True, status='exportable', after=after, jobList=jobList, inputFileList=inputFileIdList, inputFileFromJobList=fromJobIdList)
        #print 'CProjectManagerDialog.compressProject jobNumberList',jobNumberList,
        if errReport.maxSeverity() > Severity.WARNING:
            warningMessage(errReport, 'Export project','Error creating XML database file', parent=self)
            return
        if jobList is not None:
            directoriesList = []
        else:
            directoriesList = ['CCP4_IMPORTED_FILES','CCP4_PROJECT_FILES']
        self.exportThread = CCP4Export.ExportProjectThread(self, projectDir=projectInfo['projectdirectory'], dbxml=dbxml, target=fileName, jobList=jobNumberList, inputFilesList=inputFilesList, directoriesList=directoriesList, extraJobList=fromJobNumberList)
        self.exportThread.finished.connect(functools.partial(self.doneSavingJobData, projectInfo['projectname'], fileName, parentWidget))
        if not blocking:
            self.exportThread.start()
        else:
            self.exportThread.run()

    @QtCore.Slot(str,str,'QWidget')
    def doneSavingJobData(self, projectName, fileName, parentWidget):
        if self.exportThread.errorReport.maxSeverity() > Severity.WARNING:
            warningMessage(self.exportThread.errorReport, 'Saving job data', 'Error saving data files for export to\n' + str(fileName), parent=parentWidget)
        self.exportThread.deleteLater()
        self.exportThread = None

    def exportJobMtzFile(self, jobId):
        from . import CCP4TaskManager
        from . import CCP4XtalData
        # Devise name for the merged file and check if it has already been created
        jobDir = self.jobDirectory(jobId=jobId, create=False)
        jobInfo = self.db().getJobInfo(jobId, ['taskname', 'jobnumber'])
        exportFile = os.path.join(jobDir, 'exportMtz.mtz')
        if os.path.exists(exportFile):
            return exportFile
        exportParams = CCP4TaskManager.TASKMANAGER().getTaskAttribute(jobInfo['taskname'], 'EXPORTMTZPARAMS', default=[])
        if len(exportParams) == 0:
            return None
        # Get the source reflection data either from aimless or an imported file
        # getSourceReflectionFile() returns a dict with elements: fileName, source, jobNumber
        reflnInfo = self.getSourceReflectionFile(jobId = jobId, jobParamNameList=exportParams[0])
        #print 'CProjectViewer.exportJobFile getSourceReflectionFile',jobInfo,reflnInfo
        if reflnInfo.get('fileName', None) is None:
            return None
        # Query database for filenames and job info for the input and ouptput objects
        fileInfo = []
        paramNameList = []
        for ep in exportParams[1:]:
            if isinstance(ep, list):
                param, lab = ep
            else:
                param = ep
                lab = None
            fInfo = self.db().getJobFilesInfo(jobId=jobId, jobParamName=param)
            if len(fInfo) == 0:
                fInfo = self.db().getJobFilesInfo(jobId=jobId, jobParamName=param, input=True)
            if len(fInfo) > 0:
                fileInfo.append(fInfo[0])
                paramNameList.append(param)
                fileInfo[-1]['label'] = fileInfo[-1]['jobnumber'] + '_' + CCP4TaskManager.TASKMANAGER().getTaskLabel(fileInfo[-1]['taskname'])
                if fileInfo[-1]['importId'] is not None:
                    fileInfo[-1]['label'] = fileInfo[-1]['label'] + '_import'
                if lab is not None:
                    fileInfo[-1]['label'] = fileInfo[-1]['label'] + '_' + lab
        colTagList = CCP4TaskManager.TASKMANAGER().exportMtzColumnLabels(taskName=jobInfo['taskname'], jobId=jobId,
                                                                         paramNameList=paramNameList, sourceInfoList=fileInfo)
        if len(colTagList) > 0:
            for ii in range(len(colTagList)):
                fileInfo[ii]['label'] = colTagList[ii]
        #print 'CProjectViewer.exportJobMtzFile fileInfo',fileInfo
        if len(fileInfo) == 0:
            return None
        # Create CAD inputFiles and command lines
        fileNo = 1
        comLines = []
        # Dont need reflnInfo['fileName'] in inputFiles as it is used for the CMtzDataFile object that runs CAD
        inputFiles = []
        for fInfo in fileInfo:
            if fInfo['fileTypeId'] in [11, 12]:
                # Beware alternative fileContent
                if fInfo.get('fileContent',None) is not None:
                    fC = fInfo['fileContent']
                else:
                    if fInfo['fileTypeId'] == 11:
                        p = CCP4XtalData.CObsDataFile(fInfo['fullPath'])
                    else:
                        p = CCP4XtalData.CPhsDataFile(fInfo['fullPath'])
                    p.setContentFlag()
                    fC = int(p.contentFlag)
            # Use CAD LABOUT line to set column labels in export file
            fileNo += 1
            inputFiles.append(fInfo['fullPath'])
            if  fInfo['fileTypeId'] == 10:
                comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=FREER_' + fInfo['label'])
            elif fInfo['fileTypeId'] == 11:
                if fC == 1:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=I(+)_' + fInfo['label'] + ' E2=SIGI(+)_' + fInfo['label'] + ' E3=I(-)_' + fInfo['label'] + ' E4=SIGI(-)_' + fInfo['label'])
                elif fC == 2:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=F(+)_' + fInfo['label'] + ' E2=SIGF(+)_' + fInfo['label'] + ' E3=F(-)_' + fInfo['label'] + ' E4=SIGF(-)_' + fInfo['label'])
                elif fC == 3:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=I_' + fInfo['label'] + ' E2=SIGI_' + fInfo['label'])
                else:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=F_' + fInfo['label'] + ' E2=SIGF_' + fInfo['label'])
            elif fInfo['fileTypeId'] == 12:
                if fC == 1:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=HLA_' + fInfo['label'] + ' E2=HLB_' + fInfo['label'] + ' E3=HLC_'+fInfo['label']+' E4=HLD_'+fInfo['label'] )
                else:
                    comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=PHI_' + fInfo['label'] + '  E2=FOM_' + fInfo['label'])
            elif fInfo['fileTypeId'] == 13:
                comLines.append('LABOUT FILENUMBER ' + str(fileNo) + ' E1=F_' + fInfo['label'] + '  E2=PHI_' + fInfo['label'])
        #print 'CProjectViewer.exportJobMtzFile',inputFiles
        #print 'CProjectViewer.exportJobMtzFile',comLines
        #  Create an CMtzDataFile object and initialise with the refln data file
        m = CCP4XtalData.CMtzDataFile(reflnInfo['fileName'])
        #print m.runCad.__doc__   #Print out docs for the function
        outfile, err = m.runCad(exportFile, inputFiles ,comLines)
        #print 'CProjectViewer.exportJobMtzFile',outfile,err.report()
        return outfile

    def getSourceReflectionFile(self, jobId=None, jobParamNameList=None):
        from . import CCP4TaskManager
        from . import CCP4XtalData
        exportTaskName = self.db().getJobInfo(jobId=jobId, mode='taskname')
        #print 'into getSourceReflectionFile',exportTaskName
        if not isinstance(jobParamNameList,list):
            jobParamNameList = [jobParamNameList]
        reflnFileList = []
        pList = []
        fList= []
        iList = []
        for jobParam in jobParamNameList:
            reflnFile = None
            fileInfoList = self.db().getJobFilesInfo(jobId=jobId, jobParamName=jobParam, input=True)
            #print 'getSourceReflectionFile fileInfoList',jobParam,fileInfoList
            if len(fileInfoList) > 0:
                taskName = self.db().getJobInfo(jobId=fileInfoList[0]['jobId'], mode='taskname')
                jobNumber =  self.db().getJobInfo(jobId=fileInfoList[0]['jobId'], mode='jobnumber')
                if taskName == 'aimless_pipe':
                    mode = 'Data reduction'
                    from ..pipelines.aimless_pipe.script import aimless_pipe
                    reflnFile = aimless_pipe.exportJobFile(jobId=fileInfoList[0]['jobId'], mode='complete_mtz')
                elif taskName == 'import_merged':
                    mode = 'Imported by import merged'
                    from ..pipelines.import_merged.script import import_merged
                    reflnFile = import_merged.exportJobFile(jobId=fileInfoList[0]['jobId'], mode='complete_mtz')
                elif taskName == 'AlternativeImportXIA2':
                    mode = 'Imported from XIA2'
                    from ..wrappers.AlternativeImportXIA2.script import AlternativeImportXIA2
                    reflnFile = AlternativeImportXIA2.exportJobFile(jobId=fileInfoList[0]['jobId'], mode='complete_mtz', fileInfo=fileInfoList[0])
                else:
                    mode = 'Imported file'
                    importInfo = self.db().getImportFileInfo(fileId=fileInfoList[0]['fileId'])
                    #print 'getSourceReflectionFile importInfo',importInfo
                    reflnFile = importInfo['sourcefilename']
                if reflnFile is not None:
                    reflnFileList.append([jobParam, reflnFile, fileInfoList[0]])
                    pList.append(jobParam)
                    fList.append(reflnFile)
                    iList.append(fileInfoList[0])
        colTagList = CCP4TaskManager.TASKMANAGER().exportMtzColumnLabels(taskName=exportTaskName, jobId=jobId, paramNameList=pList, sourceInfoList=iList)
        if len(colTagList) == 0:
            colTagList = pList
        if len(fList) == 0:
            return {}
        if len(fList) == 1:
            return {'fileName': fList[0], 'source': mode, 'jobNumber' : jobNumber}
        else:
            # Multiple reflection input need to be cadded and colmn labels sorted
            # Likely to have redundant freer?
            #print 'CProjectManager.getSourceReflectionFile reflnFileList',reflnFileList
            jobDir = self.jobDirectory(jobId=jobId, create=False)
            exportFile = os.path.join(jobDir, 'exportMtz_reflns.mtz')
            fileNo = 0
            comLines = []
            inputFiles = []
            for ii in range(len(fList)):
                if ii == 0:
                    exObj = CCP4XtalData.CMtzDataFile(fList[ii])
                else:
                    inputFiles.append(fList[ii])
                refnObj = CCP4XtalData.CMtzDataFile(fList[ii])
                comLine = 'LABOUT FILENUMBER ' + str(ii + 1)
                jj = 0
                for column in refnObj.fileContent.listOfColumns:
                    jj += 1
                    comLine += ' E' + str(jj)+'=' + str(column.columnLabel) + '_' + colTagList[ii]
                comLines.append(comLine)
            #print 'CProjectManager.getSourceReflectionFile input files',inputFiles
            outfile, err = exObj.runCad(exportFile, inputFiles ,comLines )
            #print 'CProjectManager.getSourceReflectionFile',outfile,err.report()
            return {'fileName': outfile, 'source': mode, 'jobNumber' : jobNumber}

    @QtCore.Slot(str)
    def cleanupExportFiles(self, jobId):
        err = CErrorReport()
        jobDir = self.jobDirectory(jobId=jobId)
        fileList = glob.glob(os.path.join(jobDir, 'exportMtz*.*'))
        #print 'cleanupExportFiles',fileList
        for fName in fileList:
            try:
                os.remove(fName)
            except:
                err.append(self.__class__, 999, details=fName)
        return err

    def autoCleanup(self):
        '''For projects accessed more recently than their last cleanup (-1day):
          Remove files from CCP4_TMP directory older than 1day
          Do purge of temporary files'''
        projectList = self.db().getProjectsToCleanup()
        #print 'autoCleanup',projectList
        self.cleanupProjectTMPDir(projectList)
        for projectId, projectDir, projectName in projectList:
            cleanup = CPurgeProject(projectId)
            cleanup.purgeProject(context='temporary')
            self.db().updateProject(projectId=projectId, key='LastCleanupTime', value=time.time())
            print('Finished temporary file cleanup of ', projectName)

    def cleanupProjectTMPDir(self, projectList):
        oneday = 86400.0  # seconds
        currtime = time.time()
        for projectId ,projectDir, projectName in projectList:
            tmpFiles = glob.glob(os.path.join(projectDir, 'CCP4_TMP', '*'))
            tmpFiles.extend(glob.glob(os.path.join(projectDir, 'CCP4_COOT', '*')))
            # If one of these files hasn't been accessed in over a day, then remove it
            for tmpFile in tmpFiles:
                try:
                    time_since_lfaccess = currtime - os.stat(tmpFile).st_atime
                    if time_since_lfaccess > oneday:
                        if os.path.isfile(tmpFile):
                            os.remove(tmpFile)
                        else:
                            os.rmdir(tmpFile)
                except:
                    print('FAILED deleting temporary file:', tmpFile)
                    pass


class CPurgeProject(CObject):
    #      0        1         2        3           4           5        6          7                   8              9            10
    #  'Unknown','Pending','Queued','Running','Interrupted','Failed','Finished','Running remotely','File holder','To delete','Unsatisfactory'
    FINISHED_JOB_STATUS = [4, 5, 6, 8, 10]
    PURGECODES = {1 : 'Scratch file', 2 : 'Scratch file potentially useful diagnotic', 3 : 'Diagnostic file',
                  4 : 'File redundant after report created', 5 : 'Intermediate data file', 6 : 'Redundant on project completion',
                  7 : 'Large file that could be deleted', 0 : 'Retained file - used to override default' }
    SEARCHLIST = [['report.previous_*.html', 1], ['report_tmp.previous_*.html', 1],
                  ['params.previous_*.xml', 1], ['diagnostic.previous_*.xml', 1], ['hklin.mtz' , 1],
                  ['*mtzsplit.log', 2], ['log_mtz*.txt', 2], ['sftools', 2],
                  ['ctruncate', 2 ], ['scratch', 2 ], ['stderr.txt', 3], ['stdout.txt', 3],
                  ['diagnostic.xml', 3], ['report_tmp.html', 4], ['program.xml', 4], ['XMLOUT.xml', 4],
                  ['log.txt', 4 ], ['*.scene.xml' , 6], ['*observed_data_as*', 6], ['*phases_as*', 6],
                  ['*%*/report.html', 6], ['*%*/tables_as_csv_files', 6 ], ['*%*/tables_as_xml_files', 6],
                  ['*%*/params.xml', 6 ], ['com.txt', 6], ['ABCDOUT.mtz',7,'ABCDOUT'], ['DIFFPHIOUT.mtz',7,'DIFFPHIOUT'], ['FPHIOUT.mtz',7,'FPHIOUT']]
    CONTEXTLOOKUP = {'script_finish' : [1, 2], 'script_finish_fail' : [1], 'script_finish_remote' : [1, 2],
                     'temporary' : [1, 2, 4], 'intermediate' : [1, 2, 4, 5], 'extended_intermediate' : [1, 2, 5, 7], 'project_complete' : [1, 2, 3, 4, 5, 6]}
    ERROR_CODES = {101 : {'description' : 'Failed to delete'}}

    def __init__(self, projectId=None, jobId=None, db=None):
        CObject.__init__(self)
        self.error = CErrorReport()
        if db is None:
            self.db = PROJECTSMANAGER().db()
        else:
            self.db = db
        self.projectId = projectId
        if self.projectId is not None:
            self.projectInfo = self.db.getProjectInfo(projectId=projectId)
        else:
            self.projectInfo = {}
        # Get jobId,taskname list for entire project or subjobs of a job
        if jobId is not None:
            self.taskLookup = self.db.getTaskNameLookup(jobId=jobId, extras=True)
        elif projectId is not None:
            self.taskLookup = self.db.getTaskNameLookup(projectId=self.projectId, extras=True)
        else:
            self.taskLookup = {}
        #print 'CPurgeProject.init',self.projectInfo
        #print 'CPurgeProject.init', self.taskLookup

    def getTaskname(self, jobNumber):
        try:
            idx = [x[0] for x in self.taskLookup].index(jobNumber)
        except:
            return None
        else:
            return self.taskLookup[idx][1]

    def getJobId(self, jobNumber):
        try:
            idx = [x[0] for x in self.taskLookup].index(jobNumber)
        except:
            return None
        else:
            return self.taskLookup[idx][2]

    def getChildJobs(self, jobNumber):
        dotCount = jobNumber.count('.') + 1
        idxList = [i for i, v in enumerate(self.taskLookup) if v[0].startswith(jobNumber + '.') and v[0].count('.') == dotCount]
        ret = []
        for idx in idxList:
            ret.append([self.taskLookup[idx][0], self.taskLookup[idx][1], self.taskLookup[idx][2]])
        #print 'getChildJobs',ret
        return ret

    def topJobs(self, lastCleanupTime=None):
        if lastCleanupTime is not None:
            idxList = [i for i, v in enumerate(self.taskLookup) if not v[0].count('.') and v[4] > lastCleanupTime]
        else:
            idxList = [i for i, v in enumerate(self.taskLookup) if not v[0].count('.')]
        ret = []
        for idx in idxList:
            ret.append([self.taskLookup[idx][0], self.taskLookup[idx][1], self.taskLookup[idx][2], self.taskLookup[idx][3]])
        return ret

    def jobNumberMatch(self, subJobNumber, jobNumberSelection):
        if jobNumberSelection == '*':
            return True
        if subJobNumber in jobNumberSelection.split(','):
            return True
        else:
            return False

    def jobNumberFromPath(self, path):
        num = ''
        path,frag = os.path.split(path)
        while frag != 'CCP4_JOBS':
            if frag[0:4] == 'job_':
                num = frag[4:] + '.' + num
            path,frag = os.path.split(path)
        return num[0:-1]

    def subJobSearchList(self, parentSearchList, subJobTaskname, subJobNumber):
        searchList = []
        for p in parentSearchList:
            if p[0].count('/'):
                taskName0,jobNoSele = p[0].split('/')[0].split('%')
                if (taskName0 == '*' or taskName0 == subJobTaskname) and self.jobNumberMatch(subJobNumber, jobNoSele):
                    searchList.append([p[0].split('/',1)[1], p[1], p[2]])
        return searchList

    def searchListIndex(self, searchList, fName):
        for idx in range(len(searchList)):
            if searchList[idx][0] == fName:
                return idx
        return -1

    def getJobDir(self, jobNumber):
        dL = [self.projectInfo['projectdirectory'], 'CCP4_JOBS']
        for j in jobNumber.split('.'):
            dL.append('job_' + j)
        #print 'CPurgeProject.getJobDir',dL
        return os.path.join(*dL)

    def purgeProject(self, severity=None, context='temporary'):
        ''' # Concept of LastCleanupTime complicated by possibility of differnt severities of cleanup
        if self.db is not None and self.projectId is not None:
          lastCleanupTime = self.db.getProjectInfo(projectId=self.projectId,mode='lastcleanuptime')
          self.db.updateProject(projectId=self.projectId,key='LastCleanupTime',value=time.time())
        else:
          lastCleanupTime = None '''
        if severity is None:
            severity = CPurgeProject.CONTEXTLOOKUP.get(context)
        for jN, tN, iN, jS in self.topJobs(lastCleanupTime=None):
            if jS in self.FINISHED_JOB_STATUS:
                self.purgeJob(jobNumber=jN, taskName=tN, jobId=iN, jobStatus=jS, purgeSubJobs=True, severity=severity, signal=False)
        if self.db is not None and self.projectId is not None:
            self.db.projectReset.emit({'projectId' : self.projectId})
        PROJECTSMANAGER().backupDB()

    def getTaskSearchList(self, taskName, severity):
        '''Default search list for this task including the generic search list'''
        # If search item does not have a third element then add one
        from . import CCP4TaskManager
        mySearchList = []
        for idx in range(len(self.SEARCHLIST)):
            try:
                if self.SEARCHLIST[idx][1] in severity:
                    if len(self.SEARCHLIST[idx]) == 2:
                        mySearchList.append([self.SEARCHLIST[idx][0], self.SEARCHLIST[idx][1], None])
                    else:
                        mySearchList.append(self.SEARCHLIST[idx])
            except:
                print('FAILED handling file purge search item', self.SEARCHLIST[idx])
        taskList = CCP4TaskManager.TASKMANAGER().getTaskAttribute(taskName, 'PURGESEARCHLIST', script=True)
        #print 'getTaskSearchList',taskName,taskList
        if taskList is not None:
            for ii in range(len(taskList)):
                if len(taskList[ii]) == 2:
                    search = [taskList[ii][0] , taskList[ii][1], None]
                else:
                    search = taskList[ii]
                idx = self.searchListIndex(mySearchList, search[0])
                if idx >= 0:
                    mySearchList[idx] = search
                else:
                    mySearchList.append(search)
        return mySearchList

    def createReport(self, jobId=None, jobStatus=None, jobNumber=None, func=None):
        from ..report import CCP4ReportGenerator
        self.generator = CCP4ReportGenerator.CReportGenerator(jobId=jobId, jobStatus=jobStatus, jobNumber=jobNumber)
        if jobStatus == 5:
            try:
                reportFile = self.generator.makeFailedReportFile(redo=True)
            except:
                pass
            func()
        else:
            #print 'createReport to makeReportFile'
            try:
                reportFile, newPageOrNewData = self.generator.makeReportFile()
                print('createReport containsPictures', self.generator.report.containsPictures())
                if False and self.generator.report.containsPictures():
                    self.generator.FinishedPictures.connect(func)
                else:
                    func()
            except CException as e:
                if e.maxSeverity() > Severity.WARNING:
                    pass
        return True

    @QtCore.Slot(str,str,str,str,bool,str,str,list,str,bool)
    def purgeJob(self, jobId=None, jobNumber=None, jobStatus=None, taskName=None, purgeSubJobs=True, severity=None, context=None,parentSearchList=[], reportMode='create', signal=True):
        #MN database access will not work if we have no db attached, e.g. if we are running script from the command line
        if (not hasattr(self, "db" )) or self.db is None: return
        if severity is None:
            severity = CPurgeProject.CONTEXTLOOKUP[context]
        if jobNumber is None or taskName is None:
            jobInfo = self.db.getJobInfo(jobId=jobId, mode=['jobnumber', 'taskname'])
            jobNumber = jobInfo['jobnumber']
            taskName = jobInfo['taskname']
        #print 'File cleanup for:',jobNumber,taskName, 'with severity',context
        # For top level jobs handle dependency on report file existing
        if jobNumber.count('.') == 0 and severity.count(4):
            reportFile = os.path.join(self.getJobDir(jobNumber), 'report.html')
            if not os.path.exists(reportFile):
                if reportMode == 'skip':
                    pass
                elif reportMode == 'create':
                    rv = self.createReport(jobId=jobId,jobStatus=jobStatus,jobNumber=jobNumber,func=functools.partial(self.purgeJob,jobId=jobId,jobNumber=jobNumber,jobStatus=jobStatus,taskName=taskName,purgeSubJobs=purgeSubJobs,severity=severity,parentSearchList=parentSearchList,reportMode='skip'))
                    return
        delList,dbDelList = self.purgeJob0(jobId=jobId,jobNumber=jobNumber,jobStatus=jobStatus,taskName=taskName,
            purgeSubJobs=purgeSubJobs,severity=severity,parentSearchList=parentSearchList)
        if len(delList) > 0:
            print('Deleting temporary files for', jobNumber, taskName, delList)
        for d in delList:
            try:
                if os.path.isdir(d):
                    shutil.rmtree(d)
                else:
                    os.remove(d)
            except:
                self.error.append(self.__class__, 101, d)
                print('ERROR deleting', d)
            else:
                pass
                #print 'deleted',d
        # Remove from database
        if len(dbDelList) > 0:
            self.db.deleteFilesOnJobNumberAndParamName(self.projectId, dbDelList)
            if signal:
                #FIXME I think that we have a real problem here. I am not sure we are in the correct process for this to have any effect.
                print("Signalling",'projectReset', {'projectId' : self.projectId})
                self.db.projectReset.emit({'projectId' : self.projectId})

    def purgeJob0(self, jobId=None, jobNumber=None, taskName=None, jobStatus=None, purgeSubJobs=True, severity=[], parentSearchList=[]):
        delList = []
        dbDelList = []
        jobDir = self.getJobDir(jobNumber)
        if taskName is None:
            self.getTaskname(jobNumber)
        #  get the default purge list for this task
        mySearchList = self.getTaskSearchList(taskName, severity)
        #print 'purgeJob0 mySearchList',taskName,mySearchList
        #print 'purgeJob0 parentSearchList',parentSearchList
        # If we have a supplementary purge list from a parent job then the items in that
        # list should override the same in the task list
        for search in parentSearchList:
            idx = self.searchListIndex(mySearchList, search[0])
            #print 'searchListIndex',search,idx
            if idx >= 0:
                mySearchList[idx] = search
            else:
                mySearchList.append(search)
        # if purgeSubJobs is True then we are doing this recursively for subjobs - loop over all child jobs
        if purgeSubJobs:
            for subJobNumber, subJobTaskname,subJobId in self.getChildJobs(jobNumber):
                # Extract any fileDefn from current purgeList that apply to the child job
                subJobSearchList = self.subJobSearchList(mySearchList, subJobTaskname, subJobNumber)
                #print 'purgeJob0 subJobSearchList',subJobNumber,subJobTaskname,subJobSearchList
                subDelList, subDbDelList = self.purgeJob0(jobId=subJobId, jobNumber=subJobNumber, taskName=subJobTaskname, jobStatus=jobStatus,
                                                          purgeSubJobs=purgeSubJobs, severity=severity, parentSearchList=subJobSearchList)
                delList.extend(subDelList)
                dbDelList.extend(subDbDelList)
        #print 'CPurgeProject.purgeJob0',jobNumber,taskName,mySearchList
        # For the items in purge list with appropriate status find the matching files and delete
        for fSearch,status,dbParam in mySearchList:
            if status in severity:
                lastLen = len(delList)
                if fSearch.count('*') + fSearch.count('?') + fSearch.count('['):
                    delList.extend(glob.glob(os.path.join(jobDir, fSearch)))
                else:
                    if os.path.exists(os.path.join(jobDir, fSearch)):
                        delList.append(os.path.join(jobDir, fSearch))
                if dbParam is not None:
                    for ii in range(lastLen, len(delList)):
                        #print 'appending dbDelList',ii,delList[ii]
                        jN = self.jobNumberFromPath(delList[ii])
                        dbDelList.append([jN, dbParam])
        #print 'delList',jobNumber,taskName,delList
        #print 'dbDelList',dbDelList
        return delList, dbDelList
