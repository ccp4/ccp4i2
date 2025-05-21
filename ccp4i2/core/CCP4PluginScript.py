"""
Copyright (C) 2010 University of York
Liz Potterton - May 2010 - plugin script base class
"""

    # All classes are subclassed from CObject which is sub-class of QtCore.QObject that supports
    # signal-slot mechanism or a CCP4 minimal alternative

import copy
import functools
import glob
import os
import re
import shutil
import sys
import traceback as traceback_module

from lxml import etree
from PySide2 import QtCore, QtWidgets

from . import CCP4Utils
from ..utils import startup
from .CCP4ErrorHandling import CErrorReport, CException, Severity
from .CCP4Modules import COMFILEPATCHMANAGER
from .CCP4Modules import PREFERENCES
from .CCP4Modules import PROCESSMANAGER
from .CCP4Modules import PROJECTSMANAGER
from .CCP4QtObject import CObject


class CPluginScript(CObject):

    finished = QtCore.Signal(dict)

    SUCCEEDED = 0
    FAILED = 1
    INTERRUPTED = 2
    MARK_TO_DELETE = 3
    UNSATISFACTORY = 4
    ERROR_CODES = {1 : {'description' : 'Data definition file not found'},
                   2 : {'description' : 'Invalid data for command line'},
                   3 : {'description' : 'No command to run external process has been provided'},
                   4 : {'description' : 'Error reading log file'},
                   5 : {'description' : 'Invalid data for command file'},
                   6 : {'description' : 'Could not load data file'},
                   7 : {'description' : 'Error writing command file'},
                   8 : {'description' : 'Attempting to set invalid working directory'},
                   9 : {'description' : 'Failed starting external process - this can be due to a number of things, but' \
                        'usually is due to the command used by subprocess/QProcess not working for some reason. Missing input' \
                        'files, bad commands, non-functional programs etc. Check log files and stdout.'},
                   10 : {'description' : 'Failed starting external process - no process id returned'},
                   11 : {'description' : 'Error running external process'},
                   12 : {'description' : 'Attempting to set invalid waitForFinish'},
                   13 : {'description' : 'Wrapper class has not reimplemented MakeCommandAndScript() method'},
                   14 : {'description' : 'Can not find specified command template'},
                   15 : {'severity' : Severity.WARNING, 'description' : 'No command line defined in MakeCommandAndScript()'},
                   16 : {'description' : 'Script file does not exist'},
                   17 : {'severity' : Severity.WARNING, 'description' : 'Error attempting to set output file name'},
                   18 : {'description' : 'Failed importing plugin module'},
                   19 : {'description' : 'Failed instantiating plugin object'},
                   20 : {'severity' : Severity.WARNING, 'description' : 'Failed finding output data to check'},
                   21 : {'severity' : Severity.WARNING, 'description' : 'Failed to register new job with database'},
                   22 : {'description' : 'Failed to pass database info to new job'},
                   23 : {'description' : 'Error saving outputData file'},
                   24 : {'description' : 'Failed to create job sub-directory'},
                   25 : {'severity' : Severity.WARNING,'description' : 'External process return code'},
                   26 : {'severity' : Severity.WARNING,'description' : 'Error saving status to param file'},
                   27 : {'description' : 'Error interpreting output data for split MTZ'},
                   28 : {'description' : 'Error interpreting input data for join MTZ'},
                   30 : {'severity' : Severity.WARNING,'description' : 'Warning converting miniMTZ to HKLIN - data not set'},
                   31 : {'description' : 'Error converting miniMTZ to HKLIN - data name not recognised'},
                   32 : {'description' : 'Error converting miniMTZ to HKLIN - failed running mtzjoin'},
                   33 : {'description' : 'Error converting HKLOUT to miniMTZ - data name not recognised'},
                   34 : {'description' : 'Error converting HKLOUT to miniMTZ - failed running mtzsplit'},
                   35 : {'description' : 'Error converting miniMTZ to HKLIN - data type conversion not possible'},
                   36 : {'description' : 'Error merging in FreeR flags - insufficient data'},
                   39 : {'description' : 'Unknown error in script'},
                   40 : {'description' : 'Process timed out after (secs)'},
                   41 : {'description' : 'Error checking input data'},
                   42 : {'description' : 'Error checking output data'},
                   43 : {'description' : 'Error processing input data'},
                   44 : {'description' : 'Error creating command line and script'},
                   45 : {'description' : 'Error in processing output files'},
                   46 : {'description' : 'Error in finish handler'},
                   47 : {'description' : 'Error in checking external process after completion'},
                   48 : {'description' : 'Error in the plugin script startProcess'},
                   49 : {'severity' : Severity.WARNING,'description' : 'Error reading sub-process log file'},
                   50 : {'severity' : Severity.WARNING,'description' : 'Failure while looking for project defaults file'},
                   51 : {'severity' : Severity.WARNING,'description' : 'Failure loading project defaults file'},
                   52 : {'description' : 'Failure creating project defaults file - no definition for parameter'},
                   53 : {'severity' : Severity.WARNING,'description' : 'Warning creating project defaults file - parameter already in file'},
                   54 : {'severity' : Severity.WARNING,'description' : 'Warning creating project defaults file - failed reading existing file'},
                   55 : {'description' : 'Error inserting information on hklin content into log file'},
                   56 : {'description' : 'External process exited with exit code != 0'},
                   57 : {'description' : 'Failed to find the command to run this program.' \
                          'Either it set incorrectly for the script or is not setup at all (ie. missing exe file or script)'},
                   310 : {'description' : 'Output object contains different number of sub-objects'},
                   311 : {'description' : 'Failure in running comparison code'},
                   312 : {'description' : "No evaluation possible with a 'dummy' task"}}
    TASKNAME = None
    TASKVERSION= None
    TASKMODULE= None
    TASKCOMMAND = None
    SUBTASKS = []
    COMTEMPLATE = None
    COMTEMPLATEFILE = None
    COMLINETEMPLATE = None
    DBOUTPUTDATA = None
    INTERRUPTABLE = False
    RESTARTABLE = False
    ASYNCHRONOUS = False
    PERFORMANCECLASS = None
    RUNEXTERNALPROCESS = True
    DESCRIPTION = None
    PURGEFILELIST= []

    def renameFile(self, src, dst):
        if os.path.exists(dst):
            os.unlink(dst)
        os.rename(src, dst)

    def __init__(self, parent=None, name=None, workDirectory=None, dummy=False, taskName=None,**kw):
        CObject.__init__(self, parent=parent, name=name)
        self._command = None
        self._interpreter = None
        self._workDirectory = workDirectory
        self._dbChildJobNumber = 0
        self._dbHandler = None
        self._dbProjectId = None
        self._dbProjectName = None
        self._dbJobNumber = None
        self._dbJobId = None
        self.commandLine = []
        self.commandScript = []
        self.messageList = []
        self.container = None
        self._runningProcessId = None
        self._errorReport = CErrorReport()
        self._timeout = -1
        self._makeHklinInput = ['# Mini-MTZ input to HKLIN:','#   {:<15} {:<10} {:<10} {:<50}'.format('Data type','parameter','job','annotation')]  #List of input to cmtzjoin to be output in the log file
        self._myProgramVersion = None
        self._programVersions = {}
        self.fileSystemWatcher = None
        self._ifAsync = copy.deepcopy(self.ASYNCHRONOUS)
        #print 'CPluginScript.init _ifAsync', self._ifAsync
        self._finishHandler = None
        self._timerList = []
        self.mainWindow = None
        self.editComFile = kw.get('editComFile', False)
        self._readyReadStandardOutputHandler = None
        if taskName is not None:
            self.TASKNAME = taskName
        if not dummy:
            e = self.loadContentsFromXml(name=self.TASKNAME, version=self.TASKVERSION)
            if len(e) > 0:
                self.extendErrorReport(e)
            if e.maxSeverity() > 2:
                raise e
        else:
            self.makeContainer()
        subDirList = self.getSubDirectories()
        #print 'CPluginScript.init subDirList',subDirList
        if len(subDirList) > 0:
            try:
                lastDir = os.path.split(subDirList[-1])[1]
                #print 'CPluginScript.__init__ lastDir',lastDir
                self._dbChildJobNumber = int(lastDir[4:].strip('0'))
            except:
                pass

    def getSubDirectories(self):
        # Beware workflow dir may have job_n_input_params.xml files
        subDirList = glob.glob(os.path.join(self.getWorkDirectory(),'job_*'))
        retList = []
        for sub in subDirList:
            if os.path.isdir(sub): retList.append(sub)
        if len(retList) > 0:
            retList.sort(key=functools.cmp_to_key(self.compareSubDir))
        return retList

    def compareSubDir(self, d1, d2):
        n1 = int(d1.split('_')[-1])
        n2 = int(d2.split('_')[-1])
        return n1-n2

    def setCommand(self,command):
        self._command = command

    def modifyCootBat(self, cootBat):
        if not os.path.splitext(cootBat)[1] == '.bat' or not os.path.exists(cootBat):
            return None
        text = CCP4Utils.readFile(cootBat)
        text0 = re.sub('start /affinity 1 coot-bin.exe %*', 'start /wait /affinity 1 coot-bin.exe %*', text)
        modFile = os.path.join(CCP4Utils.getDotDirectory(), 'runwincoot.bat')
        CCP4Utils.saveFile(modFile, text0, overwrite=True)
        return modFile

    def getCommand(self, exeName=None):
        if not self.RUNEXTERNALPROCESS:
            return None
        if exeName is None:
            if self._command is not None:
                exeName = self._command
            else:
                exeName = self.TASKCOMMAND
        #print 'CPluginScript.getCommand exeName',exeName
        if exeName is None:
            return None
        if exeName == 'coot':
            exePath = str(PREFERENCES().COOT_EXECUTABLE)
            if sys.platform == 'win32':
                altpath = self.modifyCootBat(exePath)
                if altpath is not None and os.path.exists(altpath):
                    exePath = altpath
                    print("Using ",exePath)
            if exePath is not None and not os.path.exists(exePath):
                exePath = None
        elif exeName == 'ccp4mg':
            exePath = str(PREFERENCES().CCP4MG_EXECUTABLE)
            if exePath is None or not  os.path.exists(exePath):
                exePath = os.path.join(os.environ['CCP4'], 'bin', 'ccp4mg')
                if not os.path.exists(exePath):
                    if sys.platform == "win32" and os.path.exists(os.path.join(os.environ["CCP4"], "bin", "ccp4mg.bat")):
                        exePath = os.path.join(os.environ["CCP4"], "bin", "ccp4mg.bat")
        elif exeName[0:5] == 'shelx':
            exePath = None
            #print 'CPluginScript.getCommand SHELXDIR',PREFERENCES().SHELXDIR
            if hasattr(PREFERENCES(),'SHELXDIR') and PREFERENCES().SHELXDIR.exists():
                exePath = os.path.join(str(PREFERENCES().SHELXDIR), exeName)
                if not os.path.exists(exePath): exePath = None
            #print 'CPluginScript.getCommand handling shelx exePath=',exePath
        else:
            exePath = PREFERENCES().EXEPATHLIST.getExecutable(exeName)
            #print 'CPluginScript.getCommand for',exeName,'using executable defined in preferences',exePath
        if exePath is None:
            return exeName
        else:
            return exePath

    command = property(getCommand, setCommand)

    def getWorkDirectory(self, ifRelPath=False):
        if self._workDirectory is not None and os.path.exists(self._workDirectory):
            if ifRelPath and self._workDirectory.count('CCP4_JOBS'):
                return self._workDirectory[self._workDirectory.index('CCP4_JOBS'):]
            else:
                return self._workDirectory
        else:
            tmpDir = CCP4Utils.getTestTmpDir()
            jobDirList = glob.glob(os.path.join(tmpDir, 'CCP4_JOBS', 'job_*'))
            num = 0
            for job in jobDirList:
                num = max(num, int(job[job.rindex('_') + 1:]))
            wrkDir = os.path.join(tmpDir, 'CCP4_JOBS', 'job_' + str(num + 1))
            try:
                print('Making temporary directory for test job:', wrkDir)
                os.mkdir(wrkDir)
            except:
                wrkDir = tmpDir
                print('Failed making temporary directory using:', wrkDir)
            self._workDirectory = wrkDir
            return wrkDir

    def projectId(self):
        return self._dbProjectId

    def projectDirectory(self):
        if self.workDirectory is not None:
            d = self.workDirectory
            b = ''
            level = 0
            while b != 'CCP4_JOBS' and level < 10:
                d,b = os.path.split(d)
                level += 1
            return d
        else:
            print('CPluginScript.projectDirectory no workDirectory set - using current working directory')
            return os.path.normpath(os.getcwd())

    def setDbData(self, handler=None, projectId=None, projectName=None, jobId=None, jobNumber=None):
        if handler is not None:
            self._dbHandler = handler
        if projectName is not None:
            self._dbProjectName = projectName
        if projectId is not None:
            self._dbProjectId = projectId
        if jobId is not None:
            self._dbJobId = jobId
        if jobNumber is not None:
            self._dbJobNumber= jobNumber

    def setWorkDirectory(self, path=None):
        if path is None:
            self.appendErrorReport(8, 'Directory: none')
            return CException(self.__class__, 8, 'Directory: none')
        try:
            absPath = os.path.abspath(path)
        except:
            self.appendErrorReport(8, 'Directory: ' + path)
            return CException(self.__class__, 8, 'Directory: ' + path)
        if not os.path.exists(absPath):
            self.appendErrorReport(9, 'Directory: ' + path)
            return CException(self.__class__, 9, 'Directory: ' + path)
        self._workDirectory = absPath

    workDirectory = property(getWorkDirectory, setWorkDirectory)

    def relPath(self,jobNumber=None):
        if jobNumber is None:
            jobNumber = self._dbJobNumber
        numList = jobNumber.split('.')
        path = os.path.join('CCP4_JOBS', 'job_' + numList[0])
        for num in numList[1:]:
            path = os.path.join(path, 'job_' + num)
        return path

    def getProgramVersions(self):
        #print 'into getProgramVersions'
        #self.setProgramVersion()
        #self._programVersions.update({ self.TASKCOMMAND : self._myProgramVersion} )
        #print 'getProgramVersions',self.TASKCOMMAND,self._myProgramVersion,self._programVersions
        return self._programVersions

    def setProgramVersion(self, searchString=None):
        if searchString is None and self.TASKCOMMAND is not None:
            searchString = self.TASKCOMMAND
        if self._myProgramVersion is None and searchString is not None:
            try:
                logFile = self.makeFileName('LOG')
                if os.path.exists(logFile):
                    text = CCP4Utils.readFile(logFile)
                    self._myProgramVersion = CCP4Utils.searchVersion(text, searchString)
                    self._programVersions.update({self.TASKCOMMAND : self._myProgramVersion})
            except:
                print('Failed finding program version for', self.TASKCOMMAND, 'using searchString', searchString)
        try:
            self.parent()._programVersions.update(self._programVersions)
        except:
            print('Failed adding program version to parent job', self.TASKCOMMAND, self._myProgramVersion)
        return self._myProgramVersion

    def getProcessId(self):
        return self._runningProcessId

    def getJobId(self):
        return self._dbJobId

    def getJobNumber(self):
        return self._dbJobNumber

    def jobNumberString(self):
        if  self._dbJobNumber is None:
            return ''
        else:
            return str(self._dbJobNumber)

    jobId = property(getJobId)
    jobNumber = property(getJobNumber)
    processId = property(getProcessId)

    def getErrorReport(self):
        return self._errorReport

    errorReport = property(getErrorReport)

    def getLastChildJobNumber(self):
        return self._dbChildJobNumber

    lastChildJobNumber = property(getLastChildJobNumber)

    def getLastChildJobDirectory(self):
        wrkDir = os.path.join(self.workDirectory,'job_' + str(self._dbChildJobNumber))
        return wrkDir

    lastChildJobDirectory = property(getLastChildJobDirectory)

    def getDbOutputData(self):
        return self.DBOUTPUTDATA

    def setWaitForFinished(self, value):
        if isinstance(value, int):
            if value <= 0:
                self._ifAsync = True
            else:
                self._timeout = value
                self._ifAsync = False
        else:
            self.appendErrorReport(12)
            return CException(self.__class__, 12)

    def getWaitForFinished(self):
        if self._ifAsync:
            return -1
        else:
            return self._timout

    waitForFinish = property(getWaitForFinished, setWaitForFinished)

    def setAsync(self,mode=True):
        self._ifAsync = mode

    def getAsync(self):
        return self._ifAsync

    doAsync =  property(getAsync, setAsync)

    def setTimeout(self, value=999999):
        self._timeout = value

    def getTimeout(self):
        return self._timeout

    timeout = property(getTimeout,setTimeout)

    def setFinishHandler(self, command=None, **kw):
        self._finishHandler = [command, kw]

    def getFinishHandler(self):
        return self._finishHandler

    finshHandler = property(getFinishHandler, setFinishHandler)

    def appendErrorReport(self, code=0, details='', name = None, label=None, cls=None, recordTime=False, stack=True, exc_info=None):
        if cls is None:
            cls = self.__class__
        if name is None:
            name = 'Error in wrapper '+ str(self.TASKNAME)
            if self.TASKVERSION is not None:
                name = name + ' ' + str(self.TASKVERSION)
        else:
            name = 'Error in wrapper ' + name
        #print 'CPluginScript.appendErrorReport',cls,'code',code
        try:
            jobId = str(self._dbJobId)
            jobDirectory = PROJECTSMANAGER().makeFileName(jobId=jobId,mode='ROOT')
            logfiles = []
            for root, subFolders, files in os.walk(jobDirectory):
                for fn in files:
                    if fn == "log.txt" and os.path.exists(os.path.join(root,fn)):
                        fileName = os.path.join(root,fn)
                        if os.path.exists(fileName):
                            logfiles.append((fileName, os.path.getmtime(fileName)))
            logfiles.sort(key=lambda tup: tup[1]) 
            text = ""
            for f in logfiles:
                theseErrors = ""
                with open(f[0]) as fh:
                    lines = fh.readlines()
                    il = 1
                    for l in lines:
                        errors = re.findall("ERROR",l,re.I)
                        if len(errors)>0:
                            theseErrors += l.rstrip("\n") + " (line "+str(il)+")\n"
                        il += 1
                if len(theseErrors)>0:
                    details += "\n\nPossible causes of job failure in log file "+f[0]+" :\n"
                    details += theseErrors
                    details += "\n"
        except:
           pass # Failed to get and analyze log files.

        self._errorReport.append(cls=cls, code=code, details=details, name=name + ':', label=None, recordTime=recordTime, stack=stack, exc_info=exc_info)

    def extendErrorReport(self, otherReport=None):
        if otherReport is not None and isinstance(otherReport, CErrorReport) and len(otherReport) > 0:
            #print 'CPluginScript.extendErrorReport',len(otherReport),otherReport.report()
            self._errorReport.extend(otherReport)

    def makeContainer(self):
        # Make a skeleton container
        from . import CCP4Container
        self.container = CCP4Container.CContainer(parent=self)
        for item in ['inputData', 'controlParameters', 'outputData']:
            subContainer = CCP4Container.CContainer(name=item)
            self.container.addObject(object=subContainer, name=item)
        self.container.addHeader()

    def loadContentsFromXml(self, name=None, version=None):
        from . import CCP4Container
        myErrorReport = CException()
        #try:
        if 1:
            from . import CCP4TaskManager
            defFile = CCP4TaskManager.TASKMANAGER().lookupDefFile(name=name, version=version)
        #except:
        #  defFile = None
        #print 'CPluginScript.loadContentsFromXml',name,defFile
        if defFile is None and self.TASKNAME is not None:
            defFile = os.path.join(self.path(), self.TASKNAME + '.def.xml')
        #print 'CPluginScript.loadContentsFromXml defFile',defFile
        if defFile is None or not os.path.exists(defFile):
            myErrorReport.append(self.__class__, 1, defFile)
            return myErrorReport
        self.container = CCP4Container.CContainer(parent=self)
        myErrorReport.extend(self.container.loadContentsFromXml(defFile, guiAdmin=True))
        return myErrorReport

    def getProjectDefaultFile(self, taskName=None, next=False):
        # Get the last (or the next) defaults file for task taskName
        paramsFileList = glob.glob(os.path.join(self.projectDirectory(), 'CCP4_PROJECT_FILES', str(taskName) + '*.params.xml'))
        if len(paramsFileList)==0:
            if next:
                return os.path.join(self.projectDirectory(), 'CCP4_PROJECT_FILES', str(taskName) + '_001.params.xml')
            else:
                return None
        paramsFileList.sort()
        if next:
            idx = ('00' + str(int(paramsFileList[-1][-14:-11])+1))[-3:]
            return os.path.join(self.projectDirectory(), 'CCP4_PROJECT_FILES', str(taskName) + '_' + idx + '.params.xml')
        else:
            return paramsFileList[-1]

    def loadProjectDefaults(self):
        err= CErrorReport()
        try:
            defaultFile = self.getProjectDefaultFile(taskName=self.TASKNAME)
        except:
            err.append(self.__class__, 50)
            return err
        if defaultFile is None:
            return err
        try:
            self.container.loadDataFromXml(defaultFile)
            print('Loading project default parameters for', self.TASKNAME, 'from', defaultFile)
        except:
            err.append(self.__class__, 51, defaultFile)
        return err

    def path(self):
        module = __import__(self.__class__.__module__)
        p = os.path.split(module.__file__)[0]
        return p

    def checkInputData(self):
        ''' Check that all items in *inputData* container are set
          This is particularly a check that files are named and exist
          Returns a list of unset items '''
        nonExFiles = self.container.inputData.nonexistantFiles()
        #print 'checkInputData',nonExFiles
        if len(nonExFiles) > 0:
            for f in nonExFiles:
                self.appendErrorReport(6, 'File: ' + f + ' ' + str(getattr(self.container.inputData, f)))
        return nonExFiles

    def checkOutputData(self, container=None):
        '''
          Ensure that where output data is a file name
          there is a sensible name - or create one
        '''
        from . import CCP4File
        myErrorReport = CErrorReport()
        if container is None:
            container = self.container
        if container is None:
            myErrorReport.append(self.__class__, 20)
        else:
            try:
                dataList = container.outputData.dataOrder()
            except:
                myErrorReport.append(self.__class__, 20)
            else:
                #print 'CCP4PluginScript.checkOutputData dataList',dataList
                #jobName = self.objectName()
                #if jobName is None:
                #  jobName = ''
                #else:
                #  jobName = str(jobName)+'_'
                jobName = ''
                for objectName in dataList:
                    #print 'CCP4PluginScript.checkOutputData objectName',objectName
                    try:
                        dobj = container.outputData.find(objectName)
                        #print 'CCP4PluginScript.checkOutputData get',objectName,
                        #print dobj.isSet()
                        if isinstance(dobj,CCP4File.CDataFile) and not dobj.isSet():
                            fullPath = os.path.join(self.getWorkDirectory(),jobName + objectName + '.' + dobj.fileExtensions()[0])
                            #print 'CPluginScript.checkOutputData path',objectName,self.getWorkDirectory(),dobj.fileExtensions(),
                            #print fullPath
                            dobj.setFullPath(fullPath)
                    except:
                        myErrorReport.append(self.__class__, 17, objectName)
                    #print 'CCP4PluginScript.checkOutputData',objectName,str(dobj)
        return myErrorReport

    def appendCommandScript(self, text=None, fileName=None, oneLine=False, clear=False):
        '''
        Add a text string or a list of text strings to the command script
        '''
        from . import CCP4Data
        myErrorReport = CErrorReport()
        if clear:
            self.commandScript = []
        if fileName is not None:
            # Deal with possible CDataFile input
            fileName = str(fileName)
            if not os.path.exists(fileName):
                myErrorReport.append(self.__class__, 16, fileName)
                self.appendErrorReport(16, fileName)
                return myErrorReport
            try:
                text = CCP4Utils.readFile(fileName)
            except CException as e:
                myErrorReport.extend(e)
                self.extendErrorReport(e)
                return myErrorReport
        if text is not None:
            if isinstance(text,(list, CCP4Data.CList)):
                if not oneLine:
                    for item in text:
                        e = self.appendCommandScript(item)
                        if len(e) > 0:
                            myErrorReport.append(e)
                    return
                else:
                    textIn = ''
                    for item in text:
                        try:
                            strText = str(item)
                            textIn = textIn + strText + ' '
                        except:
                            myErrorReport.append(self.__class__, 5)
                    if len(textIn) > 0:
                        textIn = textIn[0:-1]
            else:
                try:
                    textIn = str(text)
                except:
                    #print 'CPluginScript.appendCommandScript: Error trying to str',text
                    myErrorReport.append(self.__class__, 5)
                    self.appendErrorReport(5)
            #print 'CPluginScript.appendCommandScript',textIn
            if len(textIn) > 0 and textIn[-1] != '\n':
                textIn = textIn + '\n'
            self.commandScript.append(textIn)
        return myErrorReport

    def clearCommandLine(self):
        self.commandLine = []

    def appendCommandLine(self, wordList=[], clear=False):
        ''' Add a text string or a list of text strings to the command line '''
        from . import CCP4Data
        myErrorReport = CErrorReport()
        if clear:
            self.clearCommandLine()
        if not isinstance(wordList, list):
            wordList = [wordList]
        for item in wordList:
            if isinstance(item, (list, CCP4Data.CList)):
                for subItem in item:
                    try:
                        myText = str(subItem)
                        self.commandLine.append(myText)
                    except:
                        myErrorReport.append(self.__class__, 2)
            else:
                try:
                    myText = str(item)
                    myText = re.sub(r'\n', ' ', myText)
                    self.commandLine.append(myText)
                except:
                    myErrorReport.append(self.__class__, 2)
        if len(myErrorReport) > 0:
            self.extendErrorReport(myErrorReport)
        #print 'CPluginScript.appendCommandLine',self.commandLine
        return myErrorReport

    def applyComFilePatches(self):
        patchSele = self.container.guiAdmin.find('patchSelection')
        #print 'CPluginScript.applyComFilePatches',patchSele,self.commandScript
        if patchSele is None or not patchSele.patch.isSet():
            return
        script = ''
        for line in self.commandScript:
            script = script + line
        try:
            script, results = COMFILEPATCHMANAGER().applyPatches(str(patchSele.patch), script)
            #print 'CPluginScript.applyComFilePatches',script,results
        except CException as e:
            self._errorReport.extend(e)
        else:
            self.commandScript = script


    def writeCommandFile(self, qualifier=None):
        ''' Convert the list of lines from self.commandScript to a command file '''
        if len(self._makeHklinInput) > 1:
            for idx in range(len(self._makeHklinInput) - 1, -1, -1):
                self.commandScript.insert(0, self._makeHklinInput[idx] + '\n')
        self.commandScript.insert(0, '# Task ' + str(self._dbJobNumber) + ' ' + self.TASKNAME + ' running ' + self.getCommand() + '\n')
        fileName = self.makeFileName('COM', qualifier=qualifier)
        try:
            CCP4Utils.saveFile(fileName=fileName, text_list=self.commandScript)
            return fileName
        except:
            self.appendErrorReport(7, 'Command file name: ' + fileName)
            return None

    def makeFileName(self, format='COM', ext='', qualifier=None):
        ''' Generate consistent names for output files
          Should give same result as PROJECTMANAGER().makeFileName() but here we
          do not call the database as we already have all necessary info '''
        defNames = {'ROOT' : '', 'PARAMS' : 'params.xml', 'JOB_INPUT' : 'input_params.xml', 'PROGRAMXML' : 'program.xml',
                    'LOG' : 'log.txt', 'STDOUT' : 'stdout.txt', 'STDERR' : 'stderr.txt', 'INTERRUPT' : 'interrupt_status.xml',
                    'DIAGNOSTIC' : 'diagnostic.xml', 'REPORT' : 'report.html', 'LOG' : 'log.txt', 'COM' : 'com.txt',
                    'MGPICDEF' : 'report.mgpic.py', 'PIC' : 'report.png', 'RVAPIXML' : 'i2.xml' }

        fileName = defNames.get(format, 'unknown.unk')
        if qualifier is not None:
            base, ext = fileName.split('.', 1)
            fileName = base + '_' + str(qualifier) + '.' + ext
        return os.path.join(self.workDirectory, fileName)

    def processInputFiles(self):
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        return CPluginScript.SUCCEEDED

    def process(self, **kw):
        ''' Check input data is set, create program command script (by calling makeCommandAndScript
        which should be implemented in sub-class and call startProcess '''
        #print 'CPluginScript.process',self.objectName()
        #self.loadProjectDefaults()
        try:
            unsetData = self.checkInputData()
        except:
            self.appendErrorReport(41)
            return self.reportStatus(CPluginScript.FAILED)
        #print 'CPluginScript.process unsetData',unsetData
        if len(unsetData) > 0:
            return self.reportStatus(CPluginScript.FAILED)
        try:
            rv = self.checkOutputData(self.container)
            #print 'CPluginScript.process unsetOutputData',e
        except Exception as e:
            self.appendErrorReport(42, exc_info=sys.exc_info())
        else:
            if len(rv) > 0:
                self.extendErrorReport(rv)
        try:
            status = self.processInputFiles()
        except CException as e:
            return self.reportStatus(CPluginScript.FAILED)
        except Exception as e:
            self.appendErrorReport(43, exc_info=sys.exc_info())
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED
        else:
            #print 'CPluginScript.process processInputFiles',status
            if status == CPluginScript.FAILED:
                return self.reportStatus(CPluginScript.FAILED)
        try:
            rv = self.makeCommandAndScript()
            #print 'CPluginScript.process commandLine',self.command,self.commandLine,self.commandScript
        except Exception as e:
            self.appendErrorReport(44, exc_info=sys.exc_info())
            return self.reportStatus(CPluginScript.FAILED)
        if self.editComFile:
            self.displayEditor()
            return
        #return apply(self.startProcess, [self.command,[self.postProcess,{}]] , kw )
        try:
            #MN...apply has been deprecated in favour of the syntax below since python 2.3
            #rv = apply(self.startProcess, [self.command] , kw )
            rv = self.startProcess(*[self.command], **kw)
        except:
            self.appendErrorReport(48, exc_info=sys.exc_info())
            return self.reportStatus(CPluginScript.FAILED)
        else:
            if rv == CPluginScript.FAILED:
                return self.reportStatus(rv)
        if not self._ifAsync:
            return self.postProcess(processId=self._runningProcessId)
        else:
            return CPluginScript.SUCCEEDED

    def startProcess(self, command=None, handler=None, **kw):
        ''' Call PROCESSMANAGER().startProcess '''
        self._runningProcessId = None
        #print 'CPluginScript.startProcess',command
        if command is None or len(command) == 0:
            self.appendErrorReport(3)
            raise CException(self.__class__, 3)
        self.applyComFilePatches()
        if len(self.commandScript) > 0:
            inputFile = self.writeCommandFile(qualifier=kw.get('fileQualifier', None))
            #print 'CPluginScript.startProcess',command,inputFile
            if inputFile is None:
                return CPluginScript.FAILED
        else:
            inputFile = None
        logFile = self.makeFileName('LOG', qualifier=kw.get('fileQualifier', None))
        #print 'CPluginScript.startProcess',self.objectName(),command,self.commandLine,'inputFile:',inputFile
        cwd = kw.get('cwd', self.workDirectory)
        # First try to check the executable directly. Py3 only. Problem will get recorded twice (helpful for diag. in tests).
        cpth = shutil.which(command)
        if cpth == None:
            self.appendErrorReport(57, 'Process name: ' + self.command, stack=False)
        try:
            if self._ifAsync:
                handler = [self.postProcess, {}]
            else:
                handler = None
            #print 'CPluginScript.startProcess',self._ifAsync,handler; sys.stdout.flush()
            readyReadStandardOutputHandler = None
            if hasattr(self, '_readyReadStandardOutputHandler'):
                readyReadStandardOutputHandler = self._readyReadStandardOutputHandler
            self._runningProcessId = PROCESSMANAGER().startProcess(command=command, interpreter=self._interpreter, args=self.commandLine, inputFile=inputFile, logFile=logFile,
                                                                               ifAsync=self._ifAsync, timeout=self._timeout, jobId=self._dbJobId, jobNumber=self._dbJobNumber,
                                                                               projectId=self._dbProjectId, handler=handler, cwd=cwd, readyReadStandardOutputHandler=readyReadStandardOutputHandler)
            #print 'CPluginScript.startProcess',self._runningProcessId; sys.stdout.flush()
        except CException as e:
            self._errorReport.extend(e)
            self.appendErrorReport(9, 'Process name: ' + self.command, stack=False)
            self.handleProcessFailed(processId=-1)
            return CPluginScript.FAILED
        #except Exception as e:
        #  self.appendErrorReport(9,'Process name: '+self.command,exc_info=sys.exc_info(),stack=False)
        #  return CPluginScript.FAILED
        return CPluginScript.SUCCEEDED
        '''
        if not self._ifAsync:
          try:
            rv = self.postProcessCheck(self._runningProcessId)
          except Exception as e:
            #Failed to run the checking but lets let it continue
            self.appendErrorReport(47,exc_info=sys.exc_info())
          #print 'CPluginScript.startProcess after postProcessCheck',self.objectName(),rv

          try:
            rv=self.processOutputFiles()
          except CException as e:
            self.extendErrorReport(e)
            return CPluginScript.FAILED
          except Exception as e:
            self.appendErrorReport(45,exc_info=sys.exc_info())
            return CPluginScript.FAILED
          else:
            return rv
        '''

    def makeCommandAndScript(self, container=None):
        ''' Convert the parameters into program command line and script
        If COMTEMPLATE is set then try using the template to create command file '''
        from . import CCP4ComTemplate
        if container is None:
            container = self.container
        if self.COMTEMPLATEFILE != None:
            fileName = interpretFile(self.COMTEMPLATEFILE)   # Ok, this looks broken. COMTEMPxxFILE seems not to be used anywhere.
            if not os.path.exists(fileName):                  # & this looks like it will crash if it is.
                self.appendErrorReport(14, fileName)
            else:
                comTemplate = CCP4ComTemplate.CComTemplate(parent=self)
                err = comTemplate.loadTemplateFromFile(fileName)
                self.extendErrorReport(err)
                text, err = comTemplate.makeComScript(container)
                self.extendErrorReport(err)
                if len(text) > 0:
                    self.commandScript.append(text)
        elif self.COMTEMPLATE != None:
            comTemplate = CCP4ComTemplate.CComTemplate(parent=self, template=self.COMTEMPLATE)
            text, err = comTemplate.makeComScript(container)
            #print 'makeCommandAndScript from COMTEMPLATE',text
            self.extendErrorReport(err)
            if len(text) > 0:
                self.commandScript.append(text)
        else:
            #print 'CPluginScript.makeCommandAndScript should be reimplemented in sub-class'
            if self.RUNEXTERNALPROCESS:
                self.appendErrorReport(13)
        if self.COMLINETEMPLATE != None:
            comTemplate = CCP4ComTemplate.CComTemplate(parent=self, template=self.COMLINETEMPLATE)
            text,err = comTemplate.makeComScript(container)
            #print 'makeCommandAndScript from COMLINETEMPLATE',text
            self.extendErrorReport(err)
            if len(text) > 0:
                wordList = text.split()
                for word in wordList:
                    self.commandLine.append(word)
        else:
            if self.RUNEXTERNALPROCESS:
                self.appendErrorReport(15)

    @QtCore.Slot(dict)
    def postProcessWrapper(self, finishStatus):
        # Emit finished signal with whatever status has been emitted by the pdbset wrapper
        self.reportStatus(finishStatus)

    def postProcess(self, processId=-1, data={}):
        ''' Default callback method after running process '''
        #print 'CPluginScript.postProcess',processId
        if self.RUNEXTERNALPROCESS:
            try:
                rv, exitStatus, exitCode = self.postProcessCheck(processId)
            except Exception as e:
                self.appendErrorReport(47, exc_info=sys.exc_info())
                return self.reportStatus(CPluginScript.FAILED)
            #print 'CPluginScript.postProcess after postProcessCheck',self.objectName(),rv
            if rv != CPluginScript.SUCCEEDED:
                self.handleProcessFailed(processId, exitStatus, exitCode)
                self.appendErrorReport(47, 'exit status and  code: ' + str(exitStatus) + ' ' + str(exitCode))
                return self.reportStatus(rv)
        try:
            rv=self.processOutputFiles()
        except Exception as e:
            self.appendErrorReport(45, exc_info=sys.exc_info())
            return self.reportStatus(CPluginScript.FAILED)
        if rv != CPluginScript.SUCCEEDED:
            return self.reportStatus(rv)
        if self._finishHandler is not None and self._finishHandler[0] is not None:
            try:
                self._finishHandler[0](*[self._dbJobId,processId], **self._finishHandler[1])
            except  Exception as e:
                self.appendErrorReport(46, exc_info=sys.exc_info())
                return self.reportStatus(CPluginScript.FAILED)
        return self.reportStatus(CPluginScript.SUCCEEDED)

    #Method to handle process fail should be reimplemented in tasks
    def handleProcessFailed(self, processId=None, exitStatus=None, exitCode=None):
        print('CPluginScript.handleProcessFailed processId', processId, 'exitStatus', exitStatus, 'exitCode', exitCode, 'self._runningProcessId', self._runningProcessId)
        return

    def reportStatus(self, finishStatus=None, updateDb=True, traceback=False):
        #print 'CPluginScript.reportStatus',self.objectName(),finishStatus,updateDb,self._dbHandler
        if traceback:
            print('\n*Traceback from CPluginScript.reportStatus after error in a script*\n')
            traceback_module.print_exc()
            print('\n')
        #,updateDb,self._dbHandler,self._dbJobId,self.container.guiAdmin.dataOrder()
        #if self.container is not None: self.saveParams(finishStatus=finishStatus)
        #self.mergeLogFiles()
        self.setProgramVersion()
        if updateDb and (self._dbHandler is not None):
            self.setOutputFileContentFlags()
            self._dbHandler.updateJobStatus(jobId=self._dbJobId, finishStatus=finishStatus, container=self.container, dbOutputData=self.getDbOutputData())
        if self.container is not None:
            try:
                if isinstance(finishStatus, dict):
                    self.saveParams(finishStatus=finishStatus.get('finishStatus'))
                else:
                    self.container.guiAdmin.jobStatus = finishStatus
            except:
                self.appendErrorReport(26)
            if isinstance(finishStatus, dict):
                self.saveParams(finishStatus=finishStatus.get('finishStatus'))
            else:
                self.saveParams(finishStatus=finishStatus)
        if isinstance(self.parent(), CPluginScript):
            self.parent( ).errorReport.extend(self.errorReport)
        self.emitFinishSignal(finishStatus)
        return finishStatus

    def setOutputFileContentFlags(self):
        from . import CCP4File
        keyList = self.container.outputData.dataOrder()
        for key in keyList:
            obj0 = self.container.outputData.__getattr__(key)
            try:
                objList, xmlText, keyValues = obj0.saveToDb()
            except:
                objList, xmlText, keyValues = [], None, {}
            for obj in objList:
                if isinstance(obj, CCP4File.CDataFile) and obj.isSet():
                    obj.setContentFlag()

    def emitFinishSignal(self, finishStatus=None):
        #print 'CPluginScript.emitFinishSignal emiting FINISHED', self, self.objectName(),finishStatus
        if isinstance(finishStatus, dict):
            self.finished.emit(finishStatus)
        else:
            status = {'jobId' : self.jobId, 'pid' : self._runningProcessId, 'finishStatus' : finishStatus}
            self.finished.emit(status)

    def postProcessCheck(self, processId):
        if processId is None or processId < 0:
            self.appendErrorReport(10)
            return CPluginScript.FAILED, None, None
        exitStatus = PROCESSMANAGER().getJobData(processId, 'exitStatus')
        exitCode = PROCESSMANAGER().getJobData(processId, 'exitCode')
        #print 'postProcessCheck exitStatus',exitCode,exitStatus,type(exitCode),type(exitStatus)
        if exitStatus != 0:
            message = 'Process: ' + PROCESSMANAGER().getJobData(processId, 'command')
            if not PROCESSMANAGER().USEQPROCESS:
                try:
                    message = message + '\nError: ' + os.strerror(exitCode)
                except:
                    pass
            # PROCESSMANAGER.startQProcess has non-Qt exitCode=101 for program fails to start
            if exitCode == 101:
                self.appendErrorReport(11, message, stack=False)
            else:
                self.appendErrorReport(9, message, stack=False)
            return CPluginScript.FAILED, exitStatus, exitCode
        elif exitCode != 0:
            message = 'Process: ' + PROCESSMANAGER().getJobData(processId, 'command')
            self.appendErrorReport(56, message, stack=False)
            return CPluginScript.FAILED, exitStatus, exitCode
        else:
            return CPluginScript.SUCCEEDED, exitStatus, exitCode

    def mergeLogFiles(self):
        logFile = self.makeFileName('LOG')
        if os.path.exists(logFile):
            if len(self._makeHklinInput) > 1:
                try:
                    logText = CCP4Utils.readFile(logFile)
                    text = ''
                    for line in self._makeHklinInput:
                        text = text + line + '\n'
                    logText = text + logText
                    CCP4Utils.saveFile(fileName=logFile, text=logText)
                except Exception as e:
                    self.appendErrorReport(55, logFile + '/n' + str(e))
            return False
        else:
            logName = os.path.split(logFile)[1]
            text = ''
            #print 'CPluginScript.mergeLogFiles',self.getSubDirectories()
            for subD in self.getSubDirectories():
                try:
                    subLog = os.path.join(subD, logName)
                    if os.path.exists(subLog):
                        #print 'CPluginScript.mergeLogFiles',subLog
                        text = text + CCP4Utils.readFile(subLog)
                except:
                    self.appendErrorReport(49, subLog)
            if len(text) > 0:
                CCP4Utils.saveFile(logFile, text)
                return True
            else:
                return False

    def logFileText(self):
        logFile = self.makeFileName('LOG')
        try:
            logText = CCP4Utils.readFile(logFile)
        except:
            self.appendErrorReport(4, 'File name' + logFile)
            logText = ''
        return logText

    def makePluginObject(self, pluginName=None, reportToDatabase=True, dummy=False, mode=0, pluginTitle=None):
        plugin = None
        # Get jobNumber and workDirectory - have to go through database if mode==1 i.e. a top level job
        if mode == 1:
            try:
                jobId = self._dbHandler.createJob(pluginName, jobTitle=pluginTitle)
                jobNumber = self._dbHandler.db.getJobInfo(jobId, mode='jobnumber')
            except:
                jobId = None
                self.errorReport.append(self.__class__, 21, pluginName)
            workDir = ''
            try:
                workDir = os.path.join(self.projectDirectory(), 'CCP4_JOBS', 'job_' + str(jobNumber))
                if not os.path.exists(workDir):
                    os.mkdir(workDir)
            except:
                self.errorReport.append(self.__class__, 24, workDir)
                workDir = self.getWorkDirectory()
            name = str(self._dbProjectName) + '_' + str(jobNumber)
        elif mode == 2:
            # The plugin inherits the current job parameters
            jobId = self._dbJobId
            jobNumber = self._dbJobNumber
            workDir = self.getWorkDirectory()
            name =  str(self._dbProjectName) + '_' + str(jobNumber)
        else:
            #Conventional create child plugin
            self._dbChildJobNumber = self._dbChildJobNumber + 1
            #name = str(self.objectName())+'_'+str(self._dbChildJobNumber)+'_'+pluginName
            name = str(self.objectName()) + '_' + str(self._dbChildJobNumber)
            #print 'makePluginObject',pluginName,dummy,name
            if self._dbJobNumber is not None:
                jobNumber=self._dbJobNumber + '.' + str(self._dbChildJobNumber)
            else:
                jobNumber = None
            try:
                workDir = os.path.join(self.getWorkDirectory(), 'job_' + str(self._dbChildJobNumber))
                if not os.path.exists(workDir):
                    os.mkdir(workDir)
            #print 'makePluginObject workDir',workDir
            except:
                self.errorReport.append(self.__class__, 24, workDir)
                workDir = self.getWorkDirectory()
        # Create instance of CPluginScript (for dummy=True) or a subclass of CPluginScript
        if dummy:
            plugin = CPluginScript(parent=self, name=name, workDirectory=workDir, dummy=True)
            #print 'makePluginObject dummy',plugin
        else:
            from . import CCP4TaskManager
            cls = CCP4TaskManager.TASKMANAGER().getPluginScriptClass(pluginName)
            #print 'makePluginObject cls',pluginName,type(pluginName),cls,name,workDir
            if cls is None:
                self.errorReport.append(self.__class__,19,pluginName)
                return None
            else:
                try:
                    plugin = cls(parent=self,name=name, workDirectory=workDir)
                except:
                    self.errorReport.append(self.__class__, 19, pluginName)
                    return None
                if pluginTitle is not None and plugin.container is not None:
                    plugin.container.header.pluginTitle = pluginTitle
        if reportToDatabase and self._dbHandler is not None:
            if mode == 0:
                try:
                    jobId = self._dbHandler.createJob(pluginName, parentJobId=self._dbJobId, jobNumber=jobNumber, jobTitle=pluginTitle)
                except:
                    jobId = None
                    self.errorReport.append(self.__class__, 21, pluginName)
            elif mode == 2:
                # Update the task info in the db
                self._dbHandler.db.updateJob(jobId, 'taskName', pluginName)
                if pluginTitle is not None:
                    self._dbHandler.db.updateJob(jobId, 'jobTitle', pluginTitle)
                else:
                    self._dbHandler.db.updateJob(jobId, 'jobTitle', plugin.TASKTITLE)

            try:
                plugin.setDbData(handler=self._dbHandler, projectName=self._dbProjectName, projectId=self._dbProjectId, jobNumber=jobNumber, jobId=jobId)
            except:
                self.errorReport.append(self.__class__, 22, pluginName)
        #print 'CPluginScript.makePluginObject',pluginName,plugin,jobId,jobNumber
        return plugin

    def terminate(self):
        # Need to terminate children
        childPluginList = self.findChildren(CPluginScript)
        #print 'CPluginScript.terminate childPluginList',childPluginList
        for obj in childPluginList:
            obj.terminate()
        PROCESSMANAGER().terminateProcess(self._runningProcessId)
        self.reportStatus(CPluginScript.FAILED)

    def updateJobStatus(self, status=None, finishStatus=None):
        if self._dbHandler is not None and self._dbJobId is not None:
            self.setOutputFileContentFlags()
            self._dbHandler.updateJobStatus(jobId=self._dbJobId, status=status, finishStatus=finishStatus, container=self.container,
                                            dbOutputData=self.getDbOutputData)

    def saveParams(self, finishStatus=None):
        fileName = self.makeFileName('PARAMS')
        if os.path.exists(fileName):
            backup = CCP4Utils.backupFile(fileName, delete=True)
        if self.container.header is not None:
            self.container.header.projectName = self._dbProjectName
            self.container.header.projectId = self._dbProjectId
            self.container.header.jobNumber = self._dbJobNumber
            self.container.header.jobId = self._dbJobId
        try:
            self.container.saveDataToXml(fileName)
        except:
            self.errorReport.append(self.__class__, 23, fileName)
        if finishStatus == CPluginScript.INTERRUPTED:
            fileName = self.makeFileName('INTERRUPT')
            #print 'saveParams interruptStatus', fileName
            backup = CCP4Utils.backupFile(fileName, delete=True)
            try:
                self.container.saveDataToXml(fileName, subContainer='interruptStatus')
            except:
                self.errorReport.append(self.__class__, 23, fileName)

    def loadInterruptStatus(self):
        restartable = getattr(self, 'RESTARTABLE', False)
        #print 'CPluginScript.loadInterruptStatus restartable',restartable
        if not restartable:
            return
        fileName = self.makeFileName('INTERRUPT')
        if not os.path.exists(fileName):
            return
        rv = self.container.interruptStatus.loadDataFromXml(fileName, check=False)
        #print 'CPluginScript.loadInterruptStatus',fileName,self.container.interruptStatus,rv
        return rv


    def testForInterrupt(self):
        if os.path.exists(os.path.join(self.getWorkDirectory(), 'INTERRUPT')):
            return True
        else:
            return False

    def convertInputMiniMtzs(self):
        #Do any of the input miniMTZs need content conversion?
        from . import CCP4XtalData
        ret = {}
        for key in self.container.inputData.dataOrder():
            obj = self.container.inputData.get(key)
            if isinstance(obj, (CCP4XtalData.CObsDataFile, CCP4XtalData.CPhsDataFile)):
                requiredContentList = obj.qualifiers('requiredContentFlag')
                if requiredContentList is not None and requiredContentList is not NotImplemented:
                    obj.setContentFlag()
                    #print 'CPluginScript.convertInputMiniMtzs',obj.objectName(),requiredContentList,obj.contentFlag
                    if obj.contentFlag not in requiredContentList:
                        for requiredContent in requiredContentList:
                            fileName, err = obj.convert(requiredContent, os.path.join(self.workDirectory, obj.objectName().__str__() + '.mtz'), parentPlugin=self)
                            #print 'CPluginScript.convertInputMiniMtzs convert',obj.objectName(),fileName
                            if fileName is not None:
                                break
                        ret[obj.objectName().__str__()] = fileName
        return ret

    '''
    columnNames = obj.columnNames()
              colin = columnNames[0]
              colout = mtzName+'_'+columnNames[0]
              for col in columnNames[1:]:
                colin = colin +','+col
                colout = colout + ',' + mtzName+'_'+col
              infiles.append([str(obj.fullPath),colin,colout])
     '''

    def makeHklin0(self, miniMtzsIn=[], hklin='hklin', ignoreErrorCodes=[]):
        #  miniMtzsIn is list of either
        #     CDataFile object name in inputData
        #  or list of [ CDataFile object name in inputData, target content type]
        error = CErrorReport()
        outfile = os.path.join(self.workDirectory, hklin + '.mtz')
        allColout = ','
        infiles = []
        for miniMtz in miniMtzsIn:
            if isinstance(miniMtz, list):
                mtzName, targetContent = miniMtz
            else:
                mtzName = miniMtz
                targetContent = None
            obj = self.container.inputData.get(mtzName)
            #print 'makeHklin obj',mtzName,targetContent,obj,repr(obj)
            if obj is None:
                self.appendErrorReport(31, mtzName)
                error.append(self.__class__, 31, mtzName)
            elif not obj.isSet():
                self.appendErrorReport(30, mtzName)
                error.append(self.__class__, 30, mtzName)
            else:
                self.appendMakeHklinInput(obj)
                if targetContent is None:
                    try:
                        colin,colout = self.makeColinColout(mtzName, obj)
                        #print 'makeHklin makeColinColout',obj.objectPath(),obj,colin,colout
                        infiles.append([str(obj.fullPath), colin, colout])
                        allColout = allColout + colout + ','
                    except:
                        pass
                else:
                    filePath = None
                    obj.setContentFlag()
                    conversion, targetContent = obj.conversion(targetContent)
                    #print 'CPluginScript.makeHklin conversion,targetContent',obj.objectName(),conversion,targetContent
                    if conversion == 'no':
                        error.append(self.__class__, 35, mtzName)
                    elif conversion == 'convert':
                        #print 'Converting data from',obj.__str__(),'to',obj.columnNames(ifString=True,content=targetContent),' using ctruncate'
                        rv =  obj.convert(targetContent, parentPlugin=self)
                        filePath,convertError = rv
                        #print 'makeHklin from convert',filePath,convertError.report()
                        #print 'makeHklin0 using '+str(filePath)+' converted by ctruncate from '+str(obj)
                        error.extend(convertError)
                        if error.maxSeverity() <= Severity.WARNING:
                            colin,colout = self.makeColinColout(mtzName, obj, targetContent)
                            infiles.append([filePath, colin, colout])
                            allColout = allColout + colout + ','
                    elif conversion == 'mtzjoin':
                        #print 'Converting data from',obj.__str__(),'to',obj.columnNames(ifString=True,content=targetContent),' using cmtzjoin'
                        colin, colout = self.makeColinColout(mtzName, obj, targetContent)
                        infiles.append([obj.__str__(), colin, colout])
                        allColout = allColout + colout + ','
                    elif conversion == 'ok':
                        colin, colout = self.makeColinColout(mtzName, obj)
                        infiles.append([obj.__str__(), colin, colout])
                        allColout = allColout + colout + ','

        #print 'Merging MTZ file with cmtzjoin outfile,infiles',outfile,infiles
        status, ret = self.joinMtz(outfile, infiles)
        if status != CPluginScript.SUCCEEDED and ret not in ignoreErrorCodes:
            error.append(CPluginScript, ret, hklin)
            self.appendErrorReport(ret, hklin, cls=CPluginScript)
        return outfile, allColout.strip(','), error

    def makeColinColout(self, mtzName, obj, targetContent=None):
        columnNames = obj.columnNames(content=targetContent)
        colin = columnNames[0]
        colout = mtzName + '_' + columnNames[0]
        for col in columnNames[1:]:
            colin = colin + ',' + col
            colout = colout + ',' + mtzName + '_' + col
        return colin, colout

    def appendMakeHklinInput(self, obj):
        jList = re.findall(r'job_[0-9]*', str(obj))
        if len(jList) > 0:
            num = jList[0][4:]
            for j in jList[1:]:
                num = num + '.' + j[4:]
        else:
            num = ' '
        appd = '#   {:<15} {:<10} {:<10} {:<50}'.format(str(obj.qualifiers('guiLabel')), str(obj.objectName()), num, str(obj.annotation))
        self._makeHklinInput.append(appd)

    def makeHklin(self, miniMtzsIn=[], hklin='hklin', ignoreErrorCodes=[]):
        #  miniMtzsIn is list of either
        #     CDataFile object name in inputData
        #  or list of [ CDataFile object name in inputData, target content type]
        #print 'makeHklin',miniMtzsIn
        error = CErrorReport()
        outfile = os.path.join(self.workDirectory, hklin + '.mtz')
        infiles = []
        for miniMtz in miniMtzsIn:
            if isinstance(miniMtz, list):
                mtzName,targetContent = miniMtz
            else:
                mtzName = miniMtz
                targetContent = None
            obj = self.container.inputData.get(mtzName)
            if obj is None:
                self.appendErrorReport(31, mtzName)
                error.append(self.__class__, 31, mtzName)
            elif not obj.isSet():
                self.appendErrorReport(30, mtzName)
                error.append(self.__class__, 30, mtzName)
            else:
                self.appendMakeHklinInput(obj)
                if targetContent is None:
                    try:
                        infiles.append([str(obj.fullPath), obj.columnNames(True)])
                    except:
                        pass
                else:
                    filePath = None
                    obj.setContentFlag()
                    conversion,targetContent = obj.conversion(targetContent)
                    if conversion == 'no':
                        error.append(self.__class__, 35, mtzName)
                    elif conversion == 'convert':
                        #print 'Converting data from',obj.__str__(),'to',obj.columnNames(ifString=True,content=targetContent),' using ctruncate'
                        rv = obj.convert(targetContent, parentPlugin=self)
                        filePath, convertError = rv
                        #print 'makeHklin from convert',filePath,convertError.report()
                        print('makeHklin using ' + str(filePath) + ' converted by ctruncate from ' + str(obj))
                        error.extend(convertError)
                        if error.maxSeverity() <= Severity.WARNING:
                            infiles.append([filePath, obj.columnNames(ifString=True, content=targetContent)])
                    elif conversion == 'mtzjoin':
                        #print 'Converting data from',obj.__str__(),'to',obj.columnNames(ifString=True,content=targetContent),' using cmtzjoin'
                        infiles.append([obj.__str__(), obj.columnNames(ifString=True, content=targetContent)])
                    elif conversion == 'ok':
                        infiles.append([obj.__str__(), obj.columnNames(True)])
        #print 'CPluginScript.makeHklin joinMtz outfile,infiles',outfile,infiles
        status,ret = self.joinMtz(outfile, infiles)
        if status != CPluginScript.SUCCEEDED and ret not in ignoreErrorCodes:
            error.append(CPluginScript, ret, hklin)
            self.appendErrorReport(ret, hklin, cls=CPluginScript)
        return outfile, error

    def makeHklInput(self, miniMtzsIn=[], hklin='hklin', ignoreErrorCodes=[], extendOutputColnames=True, useInputColnames=False):
        # This function takes a list of mini-mtz files and runs mtzjoin to merge them into one file (the two flags are used to 
        # adjust the input to mtzjoin, the first adds the mtz file type to the output column names, the second uses fixed input 
        # column names). False, False is equiv. to the old makeHlin ftn and True, True is the same as makeHklin0. 
        # It's advised to use the defaults for new interfaces. The input miniMtzsIn is list of either CDataFile object names 
        # in inputData or list of [CDataFile object names in inputData, target content type]
        error = CErrorReport()
        infiles = []
        allColout = ','
        outfile = os.path.join(self.workDirectory, hklin + '.mtz')
        for miniMtz in miniMtzsIn:
            if isinstance(miniMtz, list):
                mtzName, targetContent = miniMtz
            else:
                mtzName = miniMtz
                targetContent = None
            obj = self.container.inputData.get(mtzName)
            if obj is None:
                self.appendErrorReport(31, mtzName)
                error.append(self.__class__, 31, mtzName)
            elif not obj.isSet():
                self.appendErrorReport(30, mtzName)
                error.append(self.__class__, 30, mtzName)
            else:
                self.appendMakeHklinInput(obj)
                self._buildInputVector(obj, mtzName, targetContent, error, extendOutputColnames, useInputColnames, infiles, allColout)
        status, ret = self.joinMtz(outfile, infiles)
        if status != CPluginScript.SUCCEEDED and ret not in ignoreErrorCodes:
            error.append(CPluginScript, ret, hklin)
            self.appendErrorReport(ret, hklin, cls=CPluginScript)
        return outfile, allColout.strip(','), error

    def _buildInputVector(self, obj, mtzName, targetContent, error, extendOutputColnames, useInputColnames, infiles, allColout):
        # Function designed to remove some of the busy-work from makeHklInput & keep it out the way. if you are writing an 
        # interface you should NOT be using this. The ftn essentially just constructs an input vector for joinMtz.
        if targetContent is None:
            try:
                infile = str(obj.fullPath)
                colin, ext_outputCol = self.makeColinColout(mtzName, obj)
                colout = obj.columnNames(True)
            except:
                pass # KJS : Keep this (unfortunately). Trying to avoid any changes cascading through i2.
        else:
            obj.setContentFlag()
            conversion, targetContent = obj.conversion(targetContent)
            infile = obj.__str__()
            if conversion == 'no':
                error.append(self.__class__, 35, mtzName)
            elif conversion == 'convert':
                rv = obj.convert(targetContent, parentPlugin=self)
                filePath = None
                filePath, convertError = rv
                error.extend(convertError)
                if error.maxSeverity() <= Severity.WARNING:
                    colin, ext_outputCol = self.makeColinColout(mtzName, obj, targetContent)
                    colout = obj.columnNames(ifString=True, content=targetContent)
                    infile = filePath
            elif conversion == 'mtzjoin':
                colin, ext_outputCol = self.makeColinColout(mtzName, obj, targetContent)
                colout = obj.columnNames(ifString=True, content=targetContent)
            elif conversion == 'ok':
                colin, ext_outputCol = self.makeColinColout(mtzName, obj)
                colout = obj.columnNames(True)
        # Save the correct mtz input to the list infiles (which later goes into mtzjoin).
        if not extendOutputColnames and not useInputColnames:   # should be equivalent to old makeHklin
            infiles.append([infile, colout])
        if extendOutputColnames and useInputColnames:           # equivalent to old makeHklin0
            infiles.append([infile, colin, ext_outputCol])
        if extendOutputColnames and not useInputColnames:
            infiles.append([infile, ext_outputCol])
        if not extendOutputColnames and useInputColnames:
            infiles.append([infile, colin, colout])
        allColout = allColout + infiles[-1][-1] + ','
        return

    def joinMtz(self, outfile, infiles):
        logFile = self.makeFileName('LOG', qualifier='mtzjoin')
        bin = shutil.which('cmtzjoin')
        arglist = ['-mtzout', outfile]
        try:
            if len(infiles[0]) == 2:
                for name, cols in infiles:
                    arglist.extend(('-mtzin', name))
                    if len(cols) > 0:
                        arglist.extend(('-colout', cols))
            else:
                for name, colin, colout in infiles:
                    arglist.extend(('-mtzin', name))
                    arglist.extend(('-colin', colin))
                    arglist.extend(('-colout', colout))
        except:
            self.appendErrorReport(28, str(infiles))
        #print 'joinMtz arglist',arglist
        pid = PROCESSMANAGER().startProcess(bin, arglist, logFile=logFile)
        status = PROCESSMANAGER().getJobData(pid)
        exitCode = PROCESSMANAGER().getJobData(pid, 'exitCode')
        #print 'CpluginScript.joinMtz',status,exitCode
        # MN Kludge because 101 errors are being generated unneccessarily by cmtzjoin when observations extend over
        # *all* reflections but FREER flags exclude systematic absences
        if status in [0, 101] and os.path.exists(outfile):
            return CPluginScript.SUCCEEDED, None
        elif status in [101]:
            return CPluginScript.UNSATISFACTORY, 36
        else:
            return CPluginScript.FAILED, 32

    def joinDicts(self, outfile, infiles):
        if(len(infiles)==0):
            return CPluginScript.SUCCEEDED, None
        if(len(infiles)==1):
            outfile.set(infiles[0])
            return CPluginScript.SUCCEEDED, None
        error = CErrorReport()
        
        outfile_string = outfile.fullPath.__str__()
        try:
            input_cif_list = []
            for dict in infiles:
               input_cif_list.append(dict.fullPath.__str__())

            from ..utils import dictionaryAccumulator
            dictionaryAccumulator.accumulate(input_cif_list,outfile_string)
            
            if not os.path.exists(outfile_string):
               raise

        except:
           print("Unable to use dictionaryAccumulator.py to merge dictionaries... will try libcheck instead...")
           print("Warning - the output dictionary from libcheck may not be self-consistent and thus may result in unexpected behaviour downstream.")

           bin = shutil.which('libcheck')
           outfile_string = outfile.fullPath.__str__()
           try:
               idx = 0
               for dictfile in infiles:
                  if(len(outfile_string)==0):
                     outfile_string = infiles[0].fullPath.__str__()
                  else:
                     outfile_old = outfile_string
                     idx = idx + 1
                     outfile_string = os.path.join(self.workDirectory, 'merged_dictionary_tmp'+str(idx))
                     print('Attempting to merge dictionary files:')
                     print('DictIn1: '+outfile_old)
                     print('DictIn2: '+dictfile.fullPath.__str__())
                     print('DictOut: '+outfile_string)
                     comFileText = '_N'+'\n'+'_FILE_L '+outfile_old+'\n'+'_FILE_L2 '+dictfile.fullPath.__str__()+'\n'+'_FILE_O '+outfile_string+'\n'+'_END'+'\n'
                     #print comFileText
                     logFile = self.makeFileName('LOG', qualifier='libcheck_'+str(idx))
                     pid = PROCESSMANAGER().startProcess(bin, inputText=comFileText, logFile=logFile)
                     status = PROCESSMANAGER().getJobData(pid)
                     exitCode = PROCESSMANAGER().getJobData(pid, 'exitCode')
                     outfile_string = os.path.join(self.workDirectory, 'merged_dictionary_tmp'+str(idx)+'.lib')
                     if status in [0, 101] and os.path.exists(outfile_string):
                        continue
                     elif status in [101]:
                        return CPluginScript.UNSATISFACTORY, 36
                     else:
                        return CPluginScript.FAILED, 32
           except:
               self.appendErrorReport(28, str(infiles))

        outfile_final = os.path.join(self.workDirectory, 'merged_dictionary.lib')
        if os.path.exists(outfile_string):
            shutil.copyfile(outfile_string, outfile_final)
        else:
            return CPluginScript.FAILED, 32
        if os.path.exists(outfile_final):
            outfile.set(outfile_final)
            print('Final merged dictionary: '+outfile.fullPath.__str__())
            return CPluginScript.SUCCEEDED, None
        else:
            return CPluginScript.FAILED, 32

    def splitHklout(self, miniMtzsOut=[], programColumnNames=[], infile=None, logFile=None):
        error = CErrorReport()
        if infile is None:
            infile = os.path.join(self.workDirectory, 'hklout.mtz')
        outfiles = []
        for ii in range(min(len(miniMtzsOut), len(programColumnNames))):
            mtzName = miniMtzsOut[ii]
            obj = self.container.outputData.get(mtzName)
            if obj is None:
                error.append(self.__class__, 33, mtzName)
                self.appendErrorReport(33, mtzName)
            elif not obj.isSet():
                pass
            else:
                outfiles.append([str(obj.fullPath), programColumnNames[ii], obj.columnNames(True)])
        status = self.splitMtz(infile, outfiles, logFile)
        if status != CPluginScript.SUCCEEDED:
            error.append(self.__class__, 34, str(outfiles))
            self.appendErrorReport(34, str(outfiles))
            return error
        print('splitHklout DONE')
        return error

    def splitHkloutList(self, miniMtzsOut=[], programColumnNames=[], outputBaseName=[], outputContentFlags=[], infileList=[], logFile=None):
        from . import CCP4File
        error = CErrorReport()
        fileTypeLabels = []
        guiLabels = []
        for  ii in range(len(miniMtzsOut)):
            objList = self.container.outputData.get(miniMtzsOut[ii])
            objList.unSet()
            #print 'splitHkloutList',miniMtzsOut[ii]#,objList.subItemQualifiers('label')
            label = objList.subItemQualifiers('label')
            if label is not NotImplemented and len(label) > 0:
                fileTypeLabels.append('_' + label)
            else:
                fileTypeLabels.append('')
            label = objList.subItemQualifiers('guiLabel')
            if label is not NotImplemented and len(label) > 0:
                guiLabels.append(label + ' ')
            else:
                guiLabels.append('Output ')
        #print 'splitHkloutList fileTypeLabels',fileTypeLabels
        indx = -1
        for infile in infileList:
            indx += 1
            if isinstance(infile, CCP4File.CDataFile):
                infilePath = infile.__str__()
            else:
                infilePath = infile
            outfiles = []

            for ii in range(min(len(miniMtzsOut), len(programColumnNames))):
                objList = self.container.outputData.get(miniMtzsOut[ii])
                while indx >= len(objList):
                    objList.append(objList.makeItem())
                    objList[-1].setFullPath(os.path.join(self.workDirectory, outputBaseName[ii] + '_' + str(indx + 1) + fileTypeLabels[ii] + '.mtz'))
                    objList[-1].annotation.set(guiLabels[ii] + str(indx + 1))
                    if len(outputContentFlags) > 0:
                        objList[-1].contentFlag.set(outputContentFlags[ii])
                outfiles.append([str(objList[-1].fullPath), programColumnNames[ii], objList[-1].columnNames(True)])
            status = self.splitMtz(infilePath, outfiles, logFile)
            if status != CPluginScript.SUCCEEDED:
                error.append(self.__class__, 34)    #, hklout)   KJS : hklout does not seem to exist here....
                self.appendErrorReport(34)          #, hklout)         Comment out.
                return error
        return error

    def splitMtz(self, infile, outfiles, logFile=None):
        print('CPluginScript.splitMtz', infile, outfiles)
        if logFile is None:
            logFile = self.makeFileName('LOG', qualifier='mtzsplit')
        bin = shutil.which('cmtzsplit')
        arglist = ['-mtzin', infile]
        try:
            for outfile in outfiles:
                if len(outfile) == 2:
                    name,colin = outfile
                    colout = ''
                else:
                    name, colin, colout = outfile
                arglist.append('-mtzout')
                arglist.append(name)
                arglist.append('-colin')
                arglist.append(colin)
                arglist.append('-colout')
                if len(colout) > 0:
                    arglist.append(colout)
                else:
                    arglist.append(colin)
        except:
            self.appendErrorReport(27, str(outfiles))
        #print 'CPluginScript.splitMtz',bin,arglist
        pid = PROCESSMANAGER().startProcess(bin, arglist, logFile=logFile)
        status = PROCESSMANAGER().getJobData(pid)
        exitCode = PROCESSMANAGER().getJobData(pid, 'exitCode')
        #print 'splitMtz',status,exitCode
        if status == CPluginScript.SUCCEEDED:
            for outfile in outfiles:
                if not os.path.exists(outfile[0]):
                    return CPluginScript.FAILED
            return CPluginScript.SUCCEEDED
        else:
            return CPluginScript.FAILED

    def timedCallback(self, timeout=60.0, handler=None, singleShot=True):
        t = QtCore.QTimer(self)
        t.setSingleShot(singleShot)
        t.setInterval(int(timeout*1000.0))
        t.timeout.connect(handler)
        if singleShot:
            t.timeout.connect(self.cleanupTimer)
        t.start()
        timerId = t.timerId()
        #print 'timedCallback',t,timerId
        return timerId

    def removeTimedCallback(self, timerId):
        #print 'CPLuginScript.removeTimedCallback'
        timerList = self.findChildren(QtCore.QTimer)
        for timer in timerList:
            if timer.timerId() == timerId:
                #print 'CPLuginScript.removeTimedCallback',timer
                timer.stop()
                timer.deleteLater()
            elif not timer.isActive():
                timer.deleteLater()

    @QtCore.Slot()
    def cleanupTimer(self):
        timerList = self.findChildren(QtCore.QTimer)
        for timer in timerList:
            if not timer.isActive():
                #print 'CPLuginScript.cleanupTimer deleting',timer
                timer.deleteLater()

    def fileWatcher(self):
        if (not hasattr(self, "fileSystemWatcher")) or self.fileSystemWatcher is None:
            self.fileSystemWatcher = QtCore.QFileSystemWatcher(parent=self)
            if PREFERENCES().FILESYSTEMWATCHERPOLLER:
                self.fileSystemWatcher.setObjectName("_qt_autotest_force_engine_poller")
            self.fileSystemWatcher.fileChanged.connect(self.filterAndDispatchFileUpdates)
            self.fileSystemWatcher.directoryChanged.connect(self.dispatchDirectoryUpdates)
        return self.fileSystemWatcher

    def watchedFiles(self):
        if not hasattr(self, "_watchedFiles"):
            self._watchedFiles = {}
        return self._watchedFiles

    def watchedDirectories(self):
        if not hasattr(self, "_watchedDirectories"):
            self._watchedDirectories = {}
        return self._watchedDirectories

    @QtCore.Slot(str)
    def filterAndDispatchFileUpdates(self, fileNameQString):
        fileName = fileNameQString.__str__()
        watchedFile = self.watchedFiles()[fileName]
        if watchedFile["minDeltaSize"] < 1:
            watchedFile["handler"](fileName)
            return
        newSize = 0
        try:
            newSize = os.stat(fileName).st_size
        except:
            return
        if (newSize-watchedFile["maxSizeYet"]) < watchedFile["minDeltaSize"]:
            return
        watchedFile["maxSizeYet"] = newSize
        if watchedFile["unwatchWhileHandling"]:
            self.fileWatcher().removePath(fileName)
        watchedFile["handler"](fileName)
        if watchedFile["unwatchWhileHandling"]:
            self.fileWatcher().addPath(fileName)
        return

    @QtCore.Slot(str)
    def dispatchDirectoryUpdates(self, directoryNameQString):
        directoryName = directoryNameQString.__str__()
        watchedDirectory = self.watchedDirectories()[directoryName]
        watchedDirectory["handler"](directoryName)
        return

    def watchFile(self, fileName, handler, minDeltaSize=0, unwatchWhileHandling=False):
        #print 'CpluginScript.watchFile',fileName,handler
        parentDirectoryPath, fileRoot = os.path.split(fileName)
        self.watchedFiles()[fileName] = {"parentDirectoryPath" : parentDirectoryPath, "handler" : handler, "maxSizeYet" : 0, "minDeltaSize" : minDeltaSize, "unwatchWhileHandling" : unwatchWhileHandling}
        if os.path.isfile(fileName):
            self.fileWatcher().addPath(fileName)
        self.watchDirectory(parentDirectoryPath, self.parentChanged)

    def parentChanged(self, directoryPath):
        for watchedFileName, watchedFile in list(self.watchedFiles().items()):
            if watchedFile["parentDirectoryPath"] == directoryPath:
                self.fileWatcher().addPath(watchedFileName)

    def watchDirectory(self, directoryName, handler):
        self.watchedDirectories()[directoryName] = {"handler" : handler}
        self.fileWatcher().addPath(directoryName)

    def displayEditor(self):
        #print 'displayEditor',self.commandScript
        if self.mainWindow is None:
            self.mainWindow = QtWidgets.QDialog()
            self.mainWindow.setModal(True)
            self.mainWindow.setLayout(QtWidgets.QVBoxLayout())
            #frame = QtWidgets.QFrame(self.mainWindow)
            #frame.setLayout(QtWidgets.QVBoxLayout())
            #self.mainWindow.setCentralWidget(frame)
            self.editor = QtWidgets.QTextEdit(self.mainWindow)
            self.mainWindow.layout().addWidget(self.editor)
            self.buttonBox = QtWidgets.QDialogButtonBox(self.mainWindow)
            b = self.buttonBox.addButton(QtWidgets.QDialogButtonBox.Save)
            b.clicked.connect(self.handleEditSave)
            b = self.buttonBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
            b.clicked.connect(self.handleEditContinue)
            self.mainWindow.layout().addWidget(self.buttonBox)
            #print 'Created edit window'
        text = ''
        for line in self.commandScript:
            if line[-1] == '\n':
                text = text + line
            else:
                text = text + line + '\n'
        self.editor.setPlainText(text)
        self.mainWindow.show()
        self.mainWindow.raise_()
        #print 'displayed edit window'

    @QtCore.Slot()
    def handleEditContinue(self):
        self.mainWindow.hide()
        self.startProcess(self.command)

    @QtCore.Slot()
    def handleEditSave(self):
        #print 'CDatabaseHandler.handleEditSave'
        text = self.editor.toPlainText()
        self.commandScript = []
        for line in text.split('\n'):
            self.commandScript.append(line)
        self.mainWindow.hide()
        self.startProcess([self.command])

    def recordInputFilesToDb(self):
        # Input files are recored by gui at run time but this handles possible
        # change to input files
        if self._dbHandler is None:
            return
        self._dbHandler.recordInputFilesToDb(jobId=self._dbJobId, container=self.container)

    def mergeDictToProjectLib(self, fileName=None, overwrite=False):
        #print 'CPluginScript.mergeDictToProjectLib',fileName
        from . import CCP4ModelData
        libFile = CCP4ModelData.CDictDataFile()
        #print 'CPluginScript.mergeDictToProjectLib',libFile,fileName
        libFile.setFullPath(libFile.defaultProjectDict(projectId=self.projectId()))
        #print 'CPluginScript.mergeDictToProjectLib',libFile,fileName
        err = libFile.fileContent.mergeFile(fileName=fileName, overwrite=overwrite)
        #print 'CPluginScript.mergeDictToProjectLib',libFile,fileName,err.report(),err.maxSeverity()
        self.extendErrorReport(err)


    def assertSame(self, otherContainer, diagnostic=False):
        from . import CCP4PerformanceData
        # Compare result of an instance of this task with another instance of the same task
        report = CErrorReport()
        # Loop over the objects that would be saved to database and create error report of all objects
        # not equivalent in self and other
        keyList = self.container.outputData.dataOrder()
        for key in keyList:
            obj0 = self.container.outputData.__getattr__(key)
            objList0, xmlText0, keyValues0 = obj0.testComparisonData()
            obj1 = otherContainer.outputData.__getattr__(key)
            objList1, xmlText1, keyValues1 = obj1.testComparisonData()
            if len(objList0) != len(objList1):
                apndtxt = str(obj0.objectPath()) + ' : ' + str(objList0)+' : ' + str(objList1)
                report.append(CPluginScript, 310, apndtxt, name=self.objectName())
            for i in range(min(len(objList0), len(objList1))):
                try:
                    result = objList0[i].assertSame(objList1[i])
                except:
                    print('Error attempting to test', objList0[i], type(objList0[i]), type(objList1[i]))
                    report.append(CPluginScript, 311, name=objList0[i].objectPath())
                else:
                    report.extend(result)
            if isinstance(obj0, CCP4PerformanceData.CPerformanceIndicator):
                try:
                    result = obj0.assertSame(obj1, diagnostic=diagnostic)
                except:
                    report.append(CPluginScript, 311, name=objList0[i].objectPath())
                else:
                    report.extend(result)
        return report

    def createWarningsXML(self, logfiles=[]):
        """
        This is an example which works well with refmac logfiles containing ATTENTION,WARNING,ERROR. 
        And at time of writing it is only called from derived prosmart_refmac plugin.
        We could probably make this more generic.
        It is up to wrapper authors to implement this appropriately and call when required.
        """
        warningsNode = etree.SubElement(self.xmlroot,"Warnings")
        for f in logfiles:
           fileNode = etree.SubElement(warningsNode,"logFile")
           fileNameNode = etree.SubElement(fileNode,"fileName")
           fileNameNode.text = f
           with open(f) as fh:
               lines = fh.readlines()
               il = 1
               for l in lines:
                   errors = re.findall("ERROR",l)
                   if len(errors)>0:
                       warningNode = etree.SubElement(fileNode,"warning")
                       typeNode = etree.SubElement(warningNode,"type")
                       warningTextNode = etree.SubElement(warningNode,"text")
                       typeNode.text = "ERROR"
                       lineNumberNode = etree.SubElement(warningNode,"lineNumber")
                       warningTextNode.text = l.rstrip("\n")
                       lineNumberNode.text = str(il)
                   warnings = re.findall("WARNING",l)
                   if len(warnings)>0:
                       warningNode = etree.SubElement(fileNode,"warning")
                       typeNode = etree.SubElement(warningNode,"type")
                       warningTextNode = etree.SubElement(warningNode,"text")
                       typeNode.text = "WARNING"
                       lineNumberNode = etree.SubElement(warningNode,"lineNumber")
                       warningTextNode.text = l.rstrip("\n")
                       lineNumberNode.text = str(il)
                   attentions = re.findall("ATTENTION",l)
                   if len(attentions)>0:
                       warningNode = etree.SubElement(fileNode,"warning")
                       typeNode = etree.SubElement(warningNode,"type")
                       warningTextNode = etree.SubElement(warningNode,"text")
                       typeNode.text = "ATTENTION"
                       lineNumberNode = etree.SubElement(warningNode,"lineNumber")
                       warningTextNode.text = l.rstrip("\n")
                       lineNumberNode.text = str(il)
                   il += 1
    
class CDatabaseHandler:

    def __init__(self, projectId=None, jobNumber=None, projectName=None):
        self.projectId = projectId
        self.projectName = projectName
        self.masterJobNumber = jobNumber
        self.db = None
        self.masterJobId = None
        # Keep list of sub-jobs for use by the file purge
        self.subJobList = []

    def openDb(self):
        from ..dbapi import CCP4DbApi
        #print 'CDatabaseHandler.openDb CCP4DbApi.CDbApi.insts',CCP4DbApi.CDbApi.insts
        try:
            if CCP4DbApi.CDbApi.insts is None:
                self.db = CCP4DbApi.CDbApi()
            else:
                self.db = CCP4DbApi.CDbApi.insts
        except Exception as e:
            print('Error opening database in CDatabaseHandler.openDb')
            print(e)
            return False
        else:
            #print 'CDatabaseHandler getting masterJobId',self.projectId,self.masterJobNumber
            try:
                self.masterJobId = self.db.getJobId(projectId=self.projectId,jobNumber=self.masterJobNumber)
            except:
                print('Error retrieving jobId for projectId,jobNumber',self.projectId,self.masterJobNumber)
                return False
            else:
                #print 'CDatabaseHandler.openDb success',self
                return True

    def createJob(self, taskName, parentJobId=None, jobNumber=None, status='Running', jobTitle=None):
        from ..dbapi import CCP4DbApi
        if status in CCP4DbApi.JOB_STATUS_TEXT:
            stat = CCP4DbApi.JOB_STATUS_TEXT.index(status)
        else:
            stat = CCP4DbApi.JOB_STATUS_PENDING
        jobId = self.db.createJob(self.projectId, taskName, jobTitle=jobTitle, parentJobId=parentJobId, jobNumber=jobNumber, status=stat)
        self.subJobList.append([jobNumber, taskName, jobId,CCP4DbApi.JOB_STATUS_FINISHED])
        #print 'CDatabaseHandler.createJob',taskName,jobId
        return jobId

    def updateJobStatus(self, jobId=None, status=None, finishStatus=None, container=None, dbOutputData=None):
        from ..dbapi import CCP4DbApi
        if status is None and finishStatus is not None:
            if isinstance(finishStatus, dict):
                finishStatus = finishStatus.get('finishStatus')
            if finishStatus == CPluginScript.SUCCEEDED:
                status = CCP4DbApi.JOB_STATUS_FINISHED
            elif finishStatus == CPluginScript.FAILED:
                status = CCP4DbApi.JOB_STATUS_FAILED
            elif finishStatus == CPluginScript.INTERRUPTED:
                status = CCP4DbApi.JOB_STATUS_INTERRUPTED
            elif finishStatus == CPluginScript.MARK_TO_DELETE:
                status = CCP4DbApi.JOB_STATUS_TO_DELETE
            elif finishStatus == CPluginScript.UNSATISFACTORY:
                status = CCP4DbApi.JOB_STATUS_UNSATISFACTORY
        if container is not None:
            #print 'CDatabaseHandler.updateJobStatus calling gleanJobFiles jobId',jobId
            try:
                e = self.db.gleanJobFiles(jobId=jobId, container=container, dbOutputData=dbOutputData, roleList=[CCP4DbApi.FILE_ROLE_OUT], unSetMissingFiles=True)
                if len(e) > 0:
                    print('Error report from extracting output data to database')
                    print(e.report())
            except:
                print('Error in CDatabaseHandler.updateJobStatus calling CDbApi.gleanJobFiles')
        try:
            self.db.updateJobStatus(jobId=jobId, status=status)
        except CException as e:
            print(e.report())
        except Exception as e:
            print(e)

    def recordInputFilesToDb(self, jobId=None, container=None):
        from ..dbapi import CCP4DbApi
        # Delete any previous recorded fileUses
        self.db.deleteFileUses(jobId=jobId)
        e = self.db.gleanJobFiles(jobId=jobId, container=container, roleList=[CCP4DbApi.FILE_ROLE_IN])
        if len(e) > 0:
            print('Error report from extracting input data to database')
            print(e.report())


class CInternalPlugin(CPluginScript):

    ERROR_CODES = {101 : { 'description' : 'No jobId or projectId provided to internal plugin'}}

    def __init__(self, parent=None, jobId=None, projectId=None, jobTitle=None, **kw):
        from ..dbapi import CCP4DbApi
        db = PROJECTSMANAGER().db()
        if jobId is None:
            if projectId is None:
                raise CException(self.__class__, 101)
            jobId = db.createJob(projectId, self.TASKNAME, parentJobId=None, jobTitle=jobTitle, status=CCP4DbApi.JOB_STATUS_PENDING)
        jobInfo = db.getJobInfo(jobId=jobId, mode=['projectid', 'projectname', 'jobnumber'])
        workDirectory = os.path.join(db.getProjectDirectory(projectId=jobInfo['projectid']), 'CCP4_JOBS', 'job_' + str(jobInfo['jobnumber']))
        if not os.path.exists(workDirectory):
            os.mkdir(workDirectory)
        name = str(jobInfo['projectname']) + '_' + str(jobInfo['jobnumber'])
        #print 'CInternalPlugin',jobInfo,workDirectory,name
        CPluginScript.__init__(self, parent=parent, name=name, workDirectory=workDirectory)
        dbHandler = CDatabaseHandler(projectId=jobInfo['projectid'], jobNumber=jobInfo['jobnumber'], projectName=jobInfo['projectname'])
        dbHandler.openDb()
        self.setDbData(handler=dbHandler, projectId=jobInfo['projectid'], jobNumber=jobInfo['jobnumber'], projectName=jobInfo['projectname'], jobId=jobId)


class CRunPlugin(CObject):

    finished = QtCore.Signal()

    ERROR_CODES = {1 : { 'description' : 'Failed importing plugin module'},
                   2 : { 'description' : 'Failed instantiating plugin object'},
                   3 : { 'description' : 'Failed running plugin process'},
                   4 : { 'description' : 'Failed assigning database handler to plugin'},
                   5 : { 'description' : 'Failed to run plugin reportStatus - database probably has wrong status'},
                   6 : { 'description' : 'Unknown error opening database'},
                   7 : { 'description' : 'Failed to find or failed to load plugin module'},
                   8 : { 'description' : 'Failed to load data from params file'},
                   9 : { 'description' : 'Failed reading compressed job data file'},
                   10 : { 'description' : 'Failed creating temporary database'},
                   11 : { 'description' : 'Failed reading com file'}}

    def __init__(self, parent=None, ccp4i2Path=None, comFilePath=None, compressedFile=None, masterWorkDir=None, dbXmlFile=None):
        CObject. __init__(self, parent)
        self.ccp4i2Path = ccp4i2Path
        self.comFilePath = comFilePath
        if compressedFile is None:
            self.compressedFile = None
        else:
            self.compressedFile = os.path.normpath(compressedFile)
        #print 'CRunPlugin.__init__ comFilePath',self.comFilePath,'dbXmlFile',dbXmlFile,'compressedFile',self.compressedFile
        if masterWorkDir is not None:
            self.masterWorkDir = os.path.normpath(masterWorkDir)
        elif comFilePath is not None:
            self.masterWorkDir = os.path.split(comFilePath)[0]
        else:
            #self.masterWorkDir = tempfile.mkdtemp(suffix='.ccp4i2_temp')
            # Need a definitely path name for remote running to find the program.xml
            # remove '_setup.ccp4db.zip' from compressed file name
            s = os.path.splitext(os.path.splitext(self.compressedFile)[0])[0]
            self.masterWorkDir = re.sub(r'_setup$',"",s) + '_work'    # self.db is set if we are running from a compressed file (i.e assumed to be 'remote' job)
            if os.path.exists(self.masterWorkDir):
                shutil.rmtree(self.masterWorkDir)
            #print 'CRunPlugin.__init__ making masterWorkDir',self.masterWorkDir
            os.mkdir(self.masterWorkDir)
        self.dbXmlFile = dbXmlFile
        #print 'CRunPlugin dbXmlFile',self.dbXmlFile
        self.db = None

    def unCompressFile(self):
        from ..qtcore import CCP4Export
        self.importThread = CCP4Export.ImportProjectThread(self, projectDir=os.path.join(self.masterWorkDir, 'project'), compressedFile=self.compressedFile)
        self.importThread.run()
        if len(self.importThread.errReport) > 0:
            print('ERROR unpacking compressed file')
            print(self.importThread.errReport.report())
        comFileList = glob.glob(os.path.join(self.masterWorkDir, 'project', 'CCP4_JOBS', 'job_*', 'input_params.xml'))
        print('comFileList', comFileList)
        if len(comFileList) >= 1:
            self.comFilePath = comFileList[0]
        else:
            print('ERROR wrong number of com files in compressed file')

    def createDb(self):
        #print 'CRunPlugin.createDb',self.masterWorkDir
        # This code is currently unused and would only be used if running a job remotely
        # in which case we need a temporary db
        # BEWARE should only ever have one CDbXml at a time as it creates temporary tables
        from ..dbapi import CCP4DbApi
        dotDir = os.path.join(self.masterWorkDir, 'dotCCP4I2')
        if os.path.exists(dotDir):
            shutil.rmtree(dotDir)
        os.mkdir(dotDir)
        os.mkdir(os.path.join(dotDir, 'db'))
        dbFile = os.path.join(dotDir, 'db', 'database.sqlite')
        pm = startup.startProjectsManager(dbFileName=dbFile)
        self.db = pm.db()
        if self.dbXmlFile is None:
            # Running from a complete compressed jobball
            self.dbXmlFile = pm.extractDatabaseXml(self.compressedFile, tempDir=self.masterWorkDir)
            projectDirectory = os.path.join(self.masterWorkDir, 'project')
        else:
            # Running remote with shared file system - projectdir is two dirs up from the job masterWorkDir
            projectDirectory = os.path.split(os.path.split(self.masterWorkDir)[0])[0]
        self.dbXml = CCP4DbApi.CDbXml(db=self.db, xmlFile=self.dbXmlFile)
        projectInfo = self.dbXml.loadProjectInfo()
        self.dbXml.projectDirectory = projectDirectory
        self.dbXml.createProject()
        commited = self.dbXml.loadTable()
        if len(self.dbXml.errReport) > 0:
            print('Import Error Report')
            self.dbXml.errReport.report()

    def setupDatabase(self, ccp4i2Path):
        # Done here rather than in startup script so the database is in the CRunPlugin thread
        #print 'CRunPlugin.setupDatabase',self.db
        if self.db is not None:
            return
        pm = startup.startProjectsManager()
        #print 'CRunPlugin.setupDatabase pm',pm,pm.db()
        self.db = pm.db()
        self.db.setDiagnostic(False)

    def compressJobData(self, jobNumber=None):
        from ..qtcore import CCP4Export
        path,base= os.path.split(self.compressedFile)
        splitBase = base.split('.')
        if splitBase[0][-6:] == '_setup':
            outFile = os.path.join(path, splitBase[0][0:-6]) + '_finished'
        else:
            outFile = os.path.join(path, splitBase[0]) + '_finished'
        for item in splitBase[1:]:
            outFile = outFile + '.' + item
        print('Saving job to:',outFile)
        try:
            PROJECTSMANAGER().cleanupJob(jobDirectory=os.path.join(self.dbXml.projectDirectory, 'CCP4_JOBS','job_' + str(jobNumber)))
        except Exception as e:
            print('ERROR cleaning up job directory before compressing\n' + str(e))
        try:
            finalDbXml = os.path.join(self.masterWorkDir, 'DATABASE_final.db.xml')
            jobNumberList, errReport = self.db.exportProjectXml(self.dbXml.projectId, fileName=finalDbXml)
            if errReport.maxSeverity() > Severity.WARNING:
                print('ERROR in exporting project database to xml\n', errReport.report())
            self.exportThread = CCP4Export.ExportProjectThread(self, projectDir=self.dbXml.projectDirectory, dbxml=finalDbXml, target=outFile, jobList=[jobNumber], directoriesList=[])
            self.exportThread.finished.connect(self.compressJobData1)
            self.exportThread.run()
            print('exportThread errorReport', self.exportThread.errorReport.report())
        except Exception as e:
            print('ERROR in exporting job files\n' + str(e))
        self.compressJobData1()
        print('DONE CRunPlugin.compressJobData')

    @QtCore.Slot()
    def compressJobData1(self):
        #print 'compressJobData1'
        self.writeFinishedFlagFile(self.compressedFile)
        self.emitFinishedSignal(0)

    def setupCom(self):
        from . import CCP4File
        # Slight kludge -- assume control file is in the work directory
        #print 'CRunPlugin.setupComAndLog',self.comFilePath
        self.workDirectory = os.path.split(self.comFilePath)[0]
        self.comFile = CCP4File.CI2XmlDataFile(self.comFilePath)
        self.comFile.loadFile()
        #print 'CRunPlugin.setupComAndLog self.comFile.header',self.comFile.header

    def setupLog(self, fileName=None):
        from . import CCP4File
        if fileName is None:
            fileName = os.path.join(self.workDirectory, 'diagnostic.xml')
        CCP4Utils.backupFile(fileName, delete=True)
        self.logFile = CCP4File.CI2XmlDataFile(fileName)
        self.logFile.header.setCurrent()
        try:
            self.logFile.header.set(self.comFile.header)
        except:
            pass
        self.logFile.header.function.set('LOG')

    def run(self, ifTrapErrors=True):
        from . import CCP4TaskManager
        self.plugin = None
        self._dbHandler = None
        self.jobNumber = None
        self.errorReport = CErrorReport()
        if self.compressedFile is not None:
            try:
                self.createDb()
            except Exception as e:
                self.errorReport.append(self.__class__, 10, details=str(e), stack=False)
            else:
                try:
                    self.unCompressFile()
                except Exception as e:
                    self.errorReport.append(self.__class__, 9, details=str(e), stack=False)
        elif self.dbXmlFile is not None:
            try:
                self.createDb()
            except Exception as e:
                self.errorReport.append(self.__class__, 10, details=str(e), stack=False)
        else:
            try:
                self.setupDatabase(self.ccp4i2Path)
            except Exception as e:
                self.errorReport.append(self.__class__, 10, details=str(e), stack=False)
        if self.errorReport.maxSeverity() > Severity.WARNING:
            self.setupLog(fileName=os.path.splitext(self.compressedFile)[0] + '.diagnostic.xml')
            self.reportFailedInitialisation()
            return
        print('CRunPlugin.run before setupComAndLog')
        try:
            self.setupCom()
            self.setupLog()
        except Exception as e:
            self.errorReport.append(self.__class__, 11, details=str(e), stack=False)
            self.reportFailedInitialisation()
            return
        self.pluginName = str(self.comFile.header.pluginName)
        projectName= str(self.comFile.header.projectName)
        projectId= str(self.comFile.header.projectId)
        self.jobNumber = str(self.comFile.header.jobNumber)
        # Initialise _dbHandler
        try:
            self._dbHandler = CDatabaseHandler(projectId= projectId, jobNumber=self.jobNumber, projectName=projectName)
            # Db is already running so _dbHandler should pick up that one
            dbOk = self._dbHandler.openDb()
        except CException as e:
            self.errorReport.extend(e)
            self.reportFailedInitialisation()
            return
        except Exception as e:
            self._dbHandler = None
            self.errorReport.append(self.__class__, 6, details=str(e), stack=False)
            self.reportFailedInitialisation()
            return
        name = str(projectName) + '_' + str(self.jobNumber)
        #print 'CRunPlugin.run',name
        cls = CCP4TaskManager.TASKMANAGER().getPluginScriptClass(self.pluginName)
        #print 'CPluginScript.run cls from TASKMANAGER',self.pluginName,cls
        if cls is None:
            self.errorReport.append(self.__class__, 7, self.pluginName)
            self.reportFailedInitialisation()
            return
        try:
            self.plugin = cls(parent=self, name=name, workDirectory=self.workDirectory, taskName=self.pluginName)
            #print 'CPluginScript.run plugin from TASKMANAGER',self.plugin
        except CException as e:
            self.errorReport.extend(e, stack=True)
            self.reportFailedInitialisation()
            return
        except Exception as e:
            self.errorReport.append(self.__class__, 2, self.pluginName, exc_info=sys.exc_info())
            self.reportFailedInitialisation()
            return
        self.projectDirectory = self.plugin.projectDirectory()
        self.plugin.finished.connect(self.postRun)
        if self.plugin.container is not None:
            try:
                rv = self.plugin.container.loadDataFromXml(str(self.comFile), guiAdmin=True, check=False)
            except:
                self.errorReport.append(self.__class__, 8, 'Params file:' + str(self.comFile), exc_info=sys.exc_info())
                self.reportFailedInitialisation()
                return
            else:
                self.errorReport.extend(rv)
        self.errorReport.extend(self.plugin.loadInterruptStatus())
        try:
            self.plugin.setDbData(handler=self._dbHandler, projectId=projectId, projectName=projectName, jobId=self._dbHandler.masterJobId, jobNumber=self.jobNumber)
        except:
            self.errorReport.append(self.__class__, 4, self.pluginName, stack=True)
        try:
            self.plugin.process()
        except CException as e:
            self.errorReport.extend(e)
            self.plugin.reportStatus(CPluginScript.FAILED)
        except Exception as e:
            self.errorReport.append(self.__class__, 3, self.pluginName, exc_info=sys.exc_info())
            self.plugin.reportStatus(CPluginScript.FAILED)
        else:
            pass

    def reportFailedInitialisation(self):
        # Early error means we never created a plugin so need to report error
        maxSeverity = self.errorReport.maxSeverity()
        #print 'CRunPlugin.run maxSeverity',maxSeverity
        from ..dbapi import CCP4DbApi
        status = CCP4DbApi.JOB_STATUS_FAILED
        #print self.errorReport.report()
        #print 'CRunPlugin.run maxSeverity',maxSeverity,'status',status
        try:
            if self.plugin is not None and hasattr(self.plugin,'container') and self._dbHandler is not None:
                self._dbHandler.updateJobStatus(jobId=self._dbHandler.masterJobId, status=status, container=self.plugin.container)
            else:
                self._dbHandler.updateJobStatus(jobId=self._dbHandler.masterJobId, status=status)
        except:
            print('Error in calling dbHandler.updateJobStatus from reportFailedInitialisation')
        self.postRun(CPluginScript.FAILED)

    @QtCore.Slot(int)
    def postRun(self, status):
        # We have called _dbHandler.updateJobStatus() in plugin.reportStatus() so we are done
        #print 'CRunPlugin.postRun',status,'dbXmlFile',self.dbXmlFile,PREFERENCES().RETAIN_DIAGNOSTIC_FILES
        if self.plugin:
            self.errorReport.extend(self.plugin.errorReport, stack=True)
        if not PREFERENCES().RETAIN_DIAGNOSTIC_FILES:
            if self.compressedFile is not None:
                self.cleanup(status=status, context='script_finish_remote')
            elif status == CPluginScript.FAILED:
                self.cleanup(status=status, context='script_finish_fail')
            else:
                self.cleanup(status=status, context='script_finish')
        self.makeLog()
        if self.compressedFile is not None and self.jobNumber is not None:
            # Construct a name for the output compressed file swapping 'setup' to 'finished'
            self.compressJobData(jobNumber=self.jobNumber)
            #print 'CRunPlugin.run DONE'; sys.stdout.flush()
            #self.writeFinishedFlagFile(self.compressedFile)
        elif self.dbXmlFile is not None:
            finalDbXml = os.path.join(os.path.split(self.dbXmlFile)[0], 'DATABASE_final.db.xml')
            jobNumberList, errReport = self.db.exportProjectXml(self.dbXml.projectId, fileName=finalDbXml)
            if errReport.maxSeverity() > Severity.WARNING:
                print(errReport.report())
            self.writeFinishedFlagFile(self.dbXmlFile)
            self.emitFinishSignal(0)
        else:
            self.emitFinishSignal(0)
        return

    def writeFinishedFlagFile(self, filename):
        # Save file as flag to remote polling mechanism that job completed
        # only write this after all other big files written
        #print 'writeFinishedFlagFile',filename
        # beware ccp4db.zip files are all in a tmp directory so need to distinguish project/job
        if filename.count('ccp4db'):
            filename = os.path.join(os.path.split(filename)[0], os.path.split(filename)[1][0:-17] + '.FINISHED')
        else:
            filename = os.path.join(os.path.split(filename)[0], 'FINISHED')
        #print 'writeFinishedFlagFile',filename
        CCP4Utils.saveFile(filename,'Temporary file to indicate to remote client that job is finished')
        self.emitFinishSignal(0)

    def makeLog(self):
        progTree = etree.Element('programVersions')
        try:
            progVersions = self.plugin.getProgramVersions()
        except:
            print('ERROR trying to get program versions')
        else:
            for pV in list(progVersions.items()):
                ele = etree.Element('programVersion')
                progTree.append(ele)
                e = etree.Element('program')
                e.text = pV[0]
                ele.append(e)
                e = etree.Element('version')
                e.text = pV[1]
                ele.append(e)
        body = self.errorReport.getEtree()
        body.append(progTree)
        try:
            self.logFile.saveFile(bodyEtree=body)
        except:
            print('Error saving diagnostic log file')
            print('The contents are:')
            try:
                print(etree.tostring(self.errorReport.getEtree(), pretty_print=True))
            except:
                print('Error writing contents')

    def cleanup(self, status=None, context='script_finish'):
        try:
            os.remove(os.path.join(self.workDirectory, 'INTERRUPT'))
        except:
            pass
        if self._dbHandler is None:
            return
        from . import CCP4ProjectsManager
        purger = CCP4ProjectsManager.CPurgeProject(projectId=self._dbHandler.projectId, db=self._dbHandler.db)
        if hasattr(self,'projectDirectory'):
            purger.projectInfo['projectdirectory'] = self.projectDirectory
        purger.taskLookup = self._dbHandler.subJobList
        #print 'CRunPlugin.cleanup',self._dbHandler.masterJobNumber,status,self.pluginName
        #print 'CRunPlugin.cleanup subJobList',self._dbHandler.subJobList
        purger.purgeJob(jobId=self._dbHandler.masterJobId, jobNumber=self._dbHandler.masterJobNumber, jobStatus=status,
                        taskName=self.pluginName, purgeSubJobs=True, context=context, reportMode='skip')
        #print 'CRunPlugin.cleanup done purge'

    def emitFinishSignal(self, rv):
        #print 'CRunPlugin.emitFinishSignal',rv
        self.finished.emit()
