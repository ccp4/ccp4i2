"""
Copyright (C) 2011 University of York
Sept 2011 Liz Potterton - rewrite using Python subprocess or Qt QProcess
"""

import os
import re
import shutil
import subprocess
import sys
import time

from PySide2 import QtCore

from . import CCP4Utils
from ..utils.QApp import QTAPPLICATION
from .CCP4Config import CONFIG
from .CCP4ErrorHandling import CErrorReport, CException, Severity
from .CCP4Preferences import PREFERENCES
from .CCP4ProjectsManager import PROJECTSMANAGER
from .CCP4TaskManager import TASKMANAGER


def PROCESSMANAGER():
    if CProcessManager.insts is None:
        CProcessManager()
    return CProcessManager.insts


#NOQ class CProcessManager():
class CProcessManager(QtCore.QObject):
    insts = None
    USEQPROCESS = True
    ERROR_CODES = {101 : {'description' : 'Error creating temporary command file'},
                   102 : {'description' : 'Process input file does not exist'},
                   103 : {'severity' : Severity.WARNING, 'description' : 'Exisiting log file has been deleted'},
                   104 : {'description' : 'Error opening input file'},
                   105 : {'description' : 'Error opening log file'},
                   106 : {'description' : 'Error starting sub-process'},
                   107 : {'description' : 'Can not run process - no executable with name'},
                   108 : {'severity' : Severity.WARNING, 'description' : 'Creating temporary log file for sub-process'},
                   109 : {'description' : 'Error opening stderr file'},
                   110 : {'description' : 'Error handling finished sub-process'},
                   111 : {'description' : 'Error calling handler after finished sub-process'}}

    def __init__(self, parent=None):
        #NOQ -- remove following lines
        if parent is None:
            parent = QTAPPLICATION()
        QtCore.QObject.__init__(self, parent)
        #NOQ -- remove above lines
        if not CProcessManager.insts:
            CProcessManager.insts = self
        self.lastProcessId = 0
        self.processInfo = {}
        self.ifAsync = False
        self.timeout = 999999
        self._processEnvironment = None
        self.runningProcesses = []
        self.pendingProcesses = []

    ##If time is positive all jobs subsequently run with startProcess() will wait for finish
    #This is useful in testing scripts with unittest. See QProcess.waitForFinish()
    #@brief All processes subsequently started with startProcess() will wait for finish
    #@param time Maximum time to wait (msec)

    def processEnvironment(self):
        if self._processEnvironment is None:
            if 'darwin' in sys.platform:
                pathVar = 'DYLD_FALLBACK_LIBRARY_PATH'
            elif 'linux' in sys.platform:
                pathVar = 'LD_LIBRARY_PATH'
            self._processEnvironment = QtCore.QProcessEnvironment.systemEnvironment()
            path = self._processEnvironment.value(pathVar).__str__()
            pathList = path.split(':')
            #print 'CProcessManager.processEnvironment path',path,pathList
            libDir = os.path.join(os.environ['CCP4MG'], 'lib')
            newPath = ''
            for item in pathList:
                try:
                    if not os.path.samefile(item, libDir):
                        newPath = newPath + item + ':'
                    else:
                        print('processEnvironment removing', item, 'from', pathVar)
                except:
                    pass
            if len(newPath) > 0:
                newPath = newPath[0:-1]
            #print 'CProcessManager.processEnvironment newPath', newPath
            self._processEnvironment.insert(pathVar, newPath)
        return self._processEnvironment

#------------------------------------------------------------------------------
    def setWaitForFinished(self, timeout=-1):
#------------------------------------------------------------------------------
        if timeout < 0:
            self.ifAsync = True
        else:
            self.ifAsync = False
            self.timeout = timeout


## Start an external process
# @param interpreter string optional name of interpreter (only 'python' supported)
# @param command string command to run
# @param args list of words Arguments for command
# @handler [callable object, dictionary of argments] The object tot call when the process completes
#   and a set of arguments (as a dictionary) to pass as the second argment to the handler. The
#   first argument to the handler is the jobId.
# @param inputFile string optional name of file with input parameters for the process
# @param inputText string optional text to input to process on stdin
# @param logFile   string optional name for log file to which stdout is directed
  
#------------------------------------------------------------------------
    def startProcess(self, command=None, args=[], inputFile=None, logFile=None, interpreter=None,
                     inputText=None, handler=[], resetEnv=True, readyReadStandardOutputHandler=None, **kw):
#------------------------------------------------------------------------
        #Use Python subprocess module or QProcess
        #print 'PROCESSMANAGER.startProcess',command,args 
        ifAsync = kw.get('ifAsync', self.ifAsync)
        timeout = kw.get('timeout', self.timeout)
        if ifAsync:
            useQProcess = True
            print(f'Running process {command} asyncronously using QProcess')
        else:
            useQProcess = False
        argList = []
        cmd = None
        if interpreter is not None:
            if interpreter == 'python':
                cmd = CCP4Utils.pythonExecutable()
            argList = [cmd, command]
        else:
            argList = [command]
            if useQProcess and 'win32' in sys.platform:
                cmd = shutil.which(command)
            else:
                cmd = command
        if isinstance(args,list):
            argList.extend(args)
        else:
            argList.extend(args.split())
        self.lastProcessId = self.lastProcessId + 1
        pid = self.lastProcessId
        self.processInfo[pid] = {'command' : cmd, 'argList' : argList, 'handler' : handler, 'readyReadStandardOutputHandler': None,
                                 'errorReport' : CErrorReport(), 'ifAsync' :ifAsync, 'timeout' : timeout, 'resetEnv' : resetEnv,
                                 'startTime' : None, 'inputFile' : None, 'logFile' : None, 'jobId' : kw.get('jobId',None),
                                 'jobNumber' : kw.get('jobNumber',None), 'projectId' : kw.get('projectId',None),
                                 'editEnv' : kw.get('editEnv',[]), 'cwd' : kw.get('cwd',None), 'finishTime' : None,
                                 'exitStatus': None, 'exitCode' : None}
        if inputFile is None and inputText is not None:
            tmpFile = None
            try:
                fd,tmpFile = CCP4Utils.makeTmpFile(mode='w')
                fd.write(inputText)
                fd.close()
                inputFile = tmpFile
            except:
                self.processInfo[pid]['errorReport'].append(self.__class__, 101, str(tmpFile))
        if inputFile is not None:
            if not os.path.exists(inputFile):
                self.processInfo[pid]['errorReport'].append(self.__class__, 102, str(inputFile))
            else:
                self.processInfo[pid]['inputFile'] = inputFile
        if logFile is None:
            logFile = CCP4Utils.makeTmpFile()
            self.processInfo[pid]['errorReport'].append(self.__class__, 108, str(logFile))
        else:
            if os.path.exists(logFile):
                os.remove(logFile)
                self.processInfo[pid]['errorReport'].append(self.__class__, 103, str(logFile))
            self.processInfo[pid]['logFile'] = logFile
        # If the invoking object has provided a callback for the "readyReadStandardOutput" signal, then
        # Provide this to the qprocess and remove the logFile provision
        if readyReadStandardOutputHandler is not None and useQProcess:
            self.processInfo[pid]['readyReadStandardOutputHandler'] = readyReadStandardOutputHandler
            self.processInfo[pid]['logFile'] = None
            useQProcess = True
        try:
            textOut = cmd
            for item in argList[1:]:
                textOut = textOut + ' ' + item
            if inputFile is not None:
                textOut = textOut + ' < ' + inputFile
            if logFile is not None:
                textOut = textOut + ' > ' + logFile
            print('\nPROCESSMANAGER running command:\n' + textOut + '\n')
            sys.stdout.flush()
        except:
            print('Error printing program run command for ' + str(cmd))
        #print 'PROCESSMANAGER.startProcess processInfo',pid,self.processInfo[pid]
        if not useQProcess:
            # method = subprocess.Popen
            # See http://bugs.python.org/issue1068268  for problem that seems to affect on OSX
            self.processInfo[pid]['startTime'] = time.time()
            try:
                callDict= {}
                if self.processInfo[pid]['inputFile'] is not None:
                    callDict['stdin'] = open(self.processInfo[pid]['inputFile'])
                if self.processInfo[pid]['logFile'] is not None:
                    callDict['stdout'] = open(self.processInfo[pid]['logFile'],'w')
                callDict['env'] = self.ccp4Env(self.processInfo[pid]['resetEnv'])
                if self.processInfo[pid]['cwd'] is not None:
                    callDict['cwd'] = self.processInfo[pid]['cwd']
                if 'win32' in sys.platform:
                    callDict['shell'] = 'True'
                #print 'calling subprocess',argList,callDict
                rv = subprocess.call(*[argList], **callDict)
                self.handleFinish(pid, rv, 0)
            except subprocess.CalledProcessError as e:
                self.processInfo[pid]['exitStatus'] = -2
                self.processInfo[pid]['exitCode'] = e.errno
            except KeyError as e:
                print('subprocess KeyError')
                self.processInfo[pid]['finishTime'] = time.time()
                self.processInfo[pid]['exitStatus'] = rv
            except Exception as e:
                print('SUBPROCESS EXCEPTION', e, type(e))
                try:
                    print(e.errno, e.strerror, os.strerror(e.errno))
                except:
                    pass
                self.processInfo[pid]['finishTime'] = time.time()
                self.processInfo[pid]['exitStatus'] = -1
                try:
                    self.processInfo[pid]['exitCode'] = e.errno
                except:
                    pass
            else:
                self.processInfo[pid]['finishTime'] = time.time()
                self.processInfo[pid]['exitStatus'] = rv
        else:
            # Use QProcess
            #print('before startQProcess', len(self.runningProcesses))
            if len(self.runningProcesses) < CONFIG().maxRunningProcesses:
                self.startQProcess(pid)
            else:
                self.pendingProcesses.append(pid)
            #print 'processManager running QProcess exiting',time.time()
        return pid

    def editProcessEnvironment(self, pid, p=None):
        #print 'editProcessEnvironment',pid,p
        if p is None:
            return
        editEnv = self.processInfo[pid]['editEnv']
        pwdDir = None
        ii=0
        while ii < len(editEnv):
            if editEnv[ii][0] in ('PWD', 'pwd'):
                pwdDir = editEnv[ii][1]
                del editEnv[ii]
            else:
                ii +=1
        # Beware handling cwd argument
        if self.processInfo[pid]['cwd']:
            pwdDir = self.processInfo[pid]['cwd']
        if self.processInfo[pid]['command'].count('coot'):
            if not pwdDir is not None and self.processInfo[pid]['projectId'] is not None:
                pwdDir = os.path.join(PROJECTSMANAGER().getProjectDirectory(projectId=self.processInfo[pid]['projectId']), 'CCP4_COOT')
                if not os.path.exists(pwdDir):
                    try:
                        os.mkdir(pwdDir)
                    except:
                        pwdDir = None
            editEnv.extend([['PYTHONHOME'], ['PYTHONPATH'], ['PYTHONSTARTUP']])
        #print 'editProcessEnvironment editEnv',editEnv
        if len(editEnv) > 0:
            processEnvironment = QtCore.QProcessEnvironment.systemEnvironment()
            for editItem in editEnv:
                if len(editItem) == 1:
                    processEnvironment.remove(editItem[0])
                else:
                    processEnvironment.insert(editItem[0], editItem[1])
            #if self.processInfo[pid]['command'].count('coot') and sys.platform == 'win32':
            #  processEnvironment = self.setCootWindowsEnvironment(processEnvironment)
            p.setProcessEnvironment(processEnvironment)
        if pwdDir is not None:
            p.setWorkingDirectory(pwdDir)

    def setCootWindowsEnvironment(self, p):
        cootDir = str(PREFERENCES().COOT_EXECUTABLE)
        COOT_GUILE_PREFIX = re.sub(r"\\\\",r"/", cootDir)
        #print 'setCootWindowsEnvironment', cootDir, COOT_GUILE_PREFIX
        coot_locations = {'COOT_PREFIX' : cootDir, 'COOT_GUILE_PREFIX' : COOT_GUILE_PREFIX, 'COOT_HOME': cootDir,
                          'COOT_BACKUP_DIR' : os.path.join(cootDir, 'coot-backup'), 'COOT_SHARE' : os.path.join(cootDir, 'share'),
                          'COOT_SCHEME_DIR' : os.path.join(cootDir, 'share', 'coot', 'scheme'),
                          'COOT_STANDARD_RESIDUES' : os.path.join(cootDir, 'share', 'coot', 'standard-residues.pdb'),
                          'COOT_PIXMAPS_DIR' : os.path.join(cootDir, 'share', 'coot', 'pixmaps'),
                          'COOT_RESOURCES_FILE' : os.path.join(cootDir, 'share', 'coot', 'cootrc'),
                          'COOT_DATA_DIR' : os.path.join(cootDir, 'share', 'coot'),
                          'COOT_REF_STRUCTS' :  os.path.join(cootDir, 'share', 'coot', 'reference-structures'),
                          'COOT_PYTHON_DIR' : os.path.join(cootDir, 'share', 'coot', 'python'),
                          'COOT_REF_SEC_STRUCTS' : os.path.join(cootDir, 'share', 'coot', 'ss-reference-structures'),
                          'PYTHONPATH' : os.path.join(cootDir, 'share', 'coot', 'python'),
                          'PYTHONHOME' : os.path.join(cootDir, 'bin'),
                          'SYMINFO': os.path.join(cootDir, 'share', 'coot', 'syminfo.lib'),
                          'GUILE_LOAD_PATH' : os.path.join(COOT_GUILE_PREFIX, 'share', 'guile', '1.8') + ';' + \
                          os.path.join(COOT_GUILE_PREFIX,'share','guile') + ';' + \
                          os.path.join(COOT_GUILE_PREFIX,'share','guile','gtk-2.0') + ';' + \
                          os.path.join(COOT_GUILE_PREFIX,'share','guile','gui')+';' + \
                          os.path.join(COOT_GUILE_PREFIX,'share','guile','www')+';' + \
                          os.path.join(COOT_GUILE_PREFIX,'share','guile','site')}
        for key,value in list(coot_locations.items()):
            p.insert(key, value)
        if p.contains('PATH'):
            p.insert('PATH', os.path.join(cootDir,'bin') + ';' + os.path.join(cootDir, 'lib') + ';' + p.value('PATH', ''))
        else:
            p.insert('PATH', os.path.join(cootDir,'bin') + ';' + os.path.join(cootDir, 'lib'))
        return p

    @QtCore.Slot(str,str)
    def printFinished(self, code, stat):
        print('startQProcess process says finished', code, stat)

    def startQProcess(self, pid):
        #print 'startQProcess',pid,self.processInfo[pid].get('logFile',None)
        self.processInfo[pid]['startTime'] = time.time()
        qArgList = []
        for item in self.processInfo[pid]['argList'][1:]:
            qArgList.append(item)
        p = QtCore.QProcess(self)
        self.editProcessEnvironment(pid,p)
        if self.processInfo[pid]['inputFile'] is not None:
            p.setStandardInputFile(self.processInfo[pid]['inputFile'])
        else:
            p.setStandardInputFile(QtCore.QProcess.nullDevice())
        if self.processInfo[pid]['logFile'] is not None:
            p.setStandardOutputFile(self.processInfo[pid]['logFile'])
        if self.processInfo[pid]['readyReadStandardOutputHandler'] is not None:
            p.readyReadStandardOutput.connect(self.processInfo[pid]['readyReadStandardOutputHandler'])
        p.start(self.processInfo[pid]['command'], qArgList)
        p.finished.connect(self.printFinished)
        p.finished.connect(lambda exitCode,exitStatus: self.handleFinish(pid,exitCode,exitStatus))
        ok = p.waitForStarted(1000)
        #print 'startQProcess waitForStarted',ok
        if not ok:
            self.handleFinish(pid, 1, 101)
            return
        if not self.processInfo[pid]['ifAsync']:
            p.waitForFinished(self.processInfo[pid]['timeout'])
        self.processInfo[pid]['qprocess'] = p
        self.runningProcesses.append(pid)
        #print 'PROCESSMANAGER.startQProcess state error',p.state(), p.error(); sys.stdout.flush()

    def runHandler(self, pid):
        handler = self.processInfo[pid]['handler']
        #print 'runHandler',pid,handler
        if handler is None or len(handler) == 0:
            return
        try:
            if len(handler) > 1:
                handler[0](*[pid], **handler[1])
            else:
                handler[0](*[pid], **{})
        except CException as e:
            self.processInfo[pid]['errorReport'].extend(e)
        except Exception as e:
            print('runHandler Error', e)
            self.processInfo[pid]['errorReport'].appendPythonException(self.__class__, str(e))

    @QtCore.Slot(str,int,int)
    def handleFinish(self, pid, exitCode=0, exitStatus=0):
        print('Process finished:', pid, 'exit code:', exitCode, 'exit status:', exitStatus,'time:', time.strftime('%H:%M:%S %d/%b/%Y', time.localtime(time.time())))
        self.processInfo[pid]['finishTime'] = time.time()
        self.processInfo[pid]['exitStatus'] = int(exitStatus)
        self.processInfo[pid]['exitCode'] = exitCode
        if "logFile" in self.processInfo[pid] and self.processInfo[pid]["logFile"]:
            if "jobId" in self.processInfo[pid] and self.processInfo[pid]["jobId"]:
                jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=self.processInfo[pid]["jobId"])
                try:
                    logFileHandle = open(self.processInfo[pid]["logFile"],'a')
                    logFileHandle.write("JOB TITLE SECTION (PROCESSMANAGER)\n")
                    if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                        logFileHandle.write(str(jobInfo["jobtitle"])+"\n")
                    else:
                        logFileHandle.write(str(TASKMANAGER().getShortTitle(jobInfo['taskname']))+"\n")
                    while "parentjobid" in jobInfo and jobInfo["parentjobid"]:
                        jobInfo = PROJECTSMANAGER().db().getJobInfo(jobId=jobInfo["parentjobid"])
                        if "jobtitle" in jobInfo and jobInfo["jobtitle"]:
                            logFileHandle.write(str(jobInfo["jobtitle"])+"\n")
                        else:
                            logFileHandle.write(str(TASKMANAGER().getShortTitle(jobInfo['taskname']))+"\n")
                    logFileHandle.close()
                except:
                    print("Could not append job title info to log file."); sys.stdout.flush()
        if self.processInfo[pid].get('qprocess') is not None:
            if pid in self.runningProcesses:
                self.runningProcesses.remove(pid)
            self.processInfo[pid]['qprocess'].deleteLater()
            del self.processInfo[pid]['qprocess']
            if len(self.pendingProcesses) > 0:
                nextPid = self.pendingProcesses.pop(0)
                self.startQProcess(nextPid)
        self.runHandler(pid)

    def getJobData(self, pid, attribute='exitStatus'):
        return self.processInfo.get(pid, {}).get(attribute, None)

    def deleteJob(self, cid):
        if cid in self.processInfo:
            del self.processInfo[cid]

    def ccp4Env(self, resetEnv=True):
        # Remove any libs in ccp4mg from path as these are built with
        # different system
        env = {}
        env.update(os.environ)
        if resetEnv:
            for envVar in ['LD_LIBRARY_PATH', 'DYLD_FALLBACK_LIBRARY_PATH']:
                if envVar in env:
                    envList = env[envVar].split(':')
                    path = ''
                    for item in envList:
                        if item.count('ccp4mg') <= 0 or item.count('pythondist') > 0:
                            path = path + ':' + item
                        else:
                            pass
                            #print 'Modified',envVar,'removing',item
                    env[envVar] = path.strip(':')
                    del env[envVar]
        return env

    def terminateProcess(self, pid):
        #print 'CProcessManager.terminateProcess',pid,self.processInfo[pid].get('qprocess',None);sys.stdout.flush()
        if self.processInfo.get(pid, None) is None:
            return 1
        if self.processInfo[pid].get('qprocess', None) is not None:
            self.processInfo[pid]['qprocess'].kill()
            #self.processInfo[pid]['qprocess'].deleteLater()
            #del self.processInfo[pid]['qprocess']
            #print 'CProcessManager.terminateProcess DONE'
            return 0
        else:
            return 2
