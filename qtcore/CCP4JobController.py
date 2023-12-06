from __future__ import print_function

"""
     CCP4JobController.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

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
   Liz Potterton May 2010 - Non-graphical controller for running scripts
   Liz Potterton May 2011 - Rewrite to use CCP4DbApi to RDMS database and separate jobs as processes
"""

import os
import sys
import re
import shutil
import copy
import functools
import subprocess
import traceback

from PySide2 import QtCore
from core import CCP4Utils
from core import CCP4Container
from core import CCP4JobServer
from core.CCP4ErrorHandling import *

class CServerSetup(CCP4Container.CContainer):

    serverSetupSaved = QtCore.Signal()

    insts = None
    
    def __init__(self,source=None):
        CCP4Container.CContainer.__init__(self,name='SERVER_SETUP')
        CServerSetup.insts = self
        self.__dict__['source'] = None
        from core import CCP4Modules
        defFile = CCP4Modules.TASKMANAGER().searchDefFile('serverSetup')
        self.loadContentsFromXml(defFile)
        self.load(source=source)

    def load(self,source=None):
        prefFile,source = self.preferencesFile(source=source)
        if os.path.exists(prefFile):
            # Beware if params file has >1 serverGroup we need to add the 2+ groups to the
            # the container as they are not in the def file 
            from core import CCP4File, CCP4Annotation
            fObj = CCP4File.CI2XmlDataFile(prefFile)
            for sGEle in fObj.getBodyEtree():
                if self.get(str(sGEle.tag)) is None:
                    self.setContents( { str(sGEle.tag) : { 'class' :CCP4Annotation.CServerGroup }} )
            self.loadDataFromXml(prefFile)
            self.__dict__['source'] = source

    def writeAccess(self,source):
        dir = os.path.split(self.preferencesFile(source)[0])[0]
        return os.access(dir , os.W_OK | os.X_OK)

    def preferencesFile(self,source=None):
        if source is None or source == 'user':
            filename = str(os.path.join(CCP4Utils.getDotDirectory(),'configs','serverSetup.params.xml'))
            if os.path.exists(filename) or source == 'user':
                return filename,'user'
        filename = str(os.path.join(CCP4Utils.getCCP4I2Dir(),'local_setup','serverSetup.params.xml'))
        return filename,'installation'

    def save(self,source=None):
        prefFile,source = self.preferencesFile(source=source)
        if os.path.exists(prefFile):
            shutil.copyfile(prefFile,prefFile+'.bak')
        #print 'CServerSetup.save',prefFile
        self.saveDataToXml(fileName=prefFile)
        self.__dict__['source'] = source
        try:
            self.serverSetupSaved.emit()
        except:
            pass


class CJobController(CCP4JobServer.CJobServer):

    serverJobFailed = QtCore.Signal(tuple)
    serversEnabledReset = QtCore.Signal(bool)
    logFileUpdated = QtCore.Signal(str)
    remoteJobUpdated = QtCore.Signal(tuple)

    insts = None
    MAXTHREADS = 3
    CHECKINTERVAL = 1000
    SERVERCHECK = 10
    REPORTCHECK = 20
    # Beware using QProcess does handle the environment customisation
    USE_QPROCESS = True
    SERVERSENABLED = None

    ERROR_CODES = {101 : {'description' : 'Control file not found'},
                   102 : {'description' : 'Error attempting to start job'},
                   103 : {'description' : 'No spare job capacity'},
                   104 : {'description' : 'Job queue does not provide next job id'},
                   105 : {'description' : 'No plugin name in control file'},
                   106 : {'description' : 'Could not find plugin'},
                   111 : {'description' : 'Request to watch job that is not running'},
                   112 : {'description' : 'Failed killing job'},
                   113 : {'description' : 'No job process found - job probably finished'},
                   114 : {'description' : 'Failed trying to kill sub-processes of the task process'},
                   115 : {'description' : 'Error in handling remote job'}}

    def __init__(self,parent=None,db=None):
        from core import CCP4Modules
        CCP4JobServer.CJobServer.__init__(self)
        if parent is None:
            parent = CCP4Modules.QTAPPLICATION()
        if db is None:
            db = CCP4Modules.PROJECTSMANAGER().db()
        QtCore.QObject.__init__(self,parent)
        if CJobController.insts is None:
            CJobController.insts = self
        self.db = db
        # Update the serverParams list if job deleted
        # When serverParams are moved to db then this should be unnecessary
        self.db.jobDeleted.connect(self.handleJobDeleted)
        self.failedOpenConnection.connect(self.handleFailTuple)
        self.failedRemoteCommand.connect(self.handleFailTuple)
        self.configFile = None
        self._watchedJobs = {}
        self._errorReport = CErrorReport()
        self._diagnostic = False
        self._dbFile = None
        self.blockExit = False

    def setDbFile(self, dbFile):
        self._dbFile = dbFile

    @QtCore.Slot(dict)
    def handleJobDeleted(self, args):
        self.deleteServerParams(args['jobId'])
        #print 'CJobController.handleJobDeleted done',self.serverParams.keys()

    def serversEnabled(self):
#Quick check is to count serverList in file
        filename = str(os.path.join(CCP4Utils.getDotDirectory(), 'configs', 'serverSetup.params.xml'))
        if os.path.exists(filename):
            with open(filename) as fconfig:
                dconfig = fconfig.read()
                if "serverList" in dconfig:
                    CJobController.SERVERSENABLED = True
                    return True

        filename = str(os.path.join(CCP4Utils.getCCP4I2Dir(), 'local_setup', 'serverSetup.params.xml'))
        if os.path.exists(filename):
            with open(filename) as fconfig:
                dconfig = fconfig.read()
                if "serverList" in dconfig:
                    CJobController.SERVERSENABLED = True
                    return True
        return False

    def resetServersEnabled(self):
        if CJobController.SERVERSENABLED is None:
            self.serversEnabled()
        else:
            previous = copy.deepcopy(CJobController.SERVERSENABLED)
            current = self.serversEnabled()
            if current != previous:
                self.serversEnabledReset.emit(current)

    def setDiagnostic(self, mode):
        self._diagnostic = mode

    def setConfigFile(self, fileName):
        # Set the config file to be passed to the sub-process
        # This is most mechanism to specify alternative config file that specifies database used for testing purposes
        self.configFile = fileName

    def startTimer(self):
        self.timer = QtCore.QTimer(self)
        self.timer.timeout.connect(self.doChecks)
        self.timer.start(CJobController.CHECKINTERVAL)

    def pollRemoteReportUpdate(self):
        for jobId, sP in list(self._serverParams.items()):
            if sP.pollReport > 0:
                return True
        return False

    def pollRemoteFinish(self):
        for jobId ,sP in list(self._serverParams.items()):
            if sP.pollFinish > 0:
                return True
        return False

    def pollRemoteFinishCount(self):
        modes = [0, 0, 0, 0, 0]
        for jobId,sP in list(self._serverParams.items()):
            if sP.pollFinish >= 0:
                modes[sP.pollFinish] += 1
        if self._diagnostic:
            print('pollRemoteFinishCount', modes)
        return modes

    @QtCore.Slot()
    def doChecks(self):
        self.blockExit = True
        self.checkJobQueue()
        if self.pollRemoteFinish():
            # Perhaps only poll the server every 2 (or more) calls to doCheck()
            self.finishPollCount = getattr(self, 'finishPollCount', 0) + 1
            if self.finishPollCount >= CJobController.SERVERCHECK:
                self.finishPollCount = 0
                modes = self.pollRemoteFinishCount()
                #print 'poll remote jobs modes',modes
                if modes[2] > 0:
                    self.pollQsubStat()
                if modes[1] > 0:
                    self.pollForFinishedFlagFile()
                if modes[3] + modes[4] > 0:
                    for jobId,sP in list(self._serverParams.items()):
                        if sP.pollFinish in [3, 4]:
                            try:
                                status = self.customHandler(jobId).pollForFinishedJob(jobId)
                                if status > 0:
                                    self.handleFinishedServerJob(jobId, status=status)
                            except:
                                pass
        if self.pollRemoteReportUpdate():
            self.reportPollCount = getattr(self, 'reportPollCount', 0) + 1
            #print 'reportPollCount',self.reportPollCount
            if self.reportPollCount >= CJobController.REPORTCHECK:
                self.reportPollCount = 0
                self.updateReports()
        """
        This is potential cause of database access error in the running script???
        if len(CCP4ProjectViewer.FILEWATCHER().jobsByUpdateInterval)>0:
          self.triggerReportCount = getattr(self,'triggerReportCount',0) + 1
          if self.triggerReportCount > 5:
            self.triggerReportCount = 0
            CCP4ProjectViewer.FILEWATCHER().triggerJobsByUpdateInterval()
        """
        self.blockExit = False
        #self.watchLogFiles()

    def checkJobQueue(self):
        # getJobsByStatus defaults to getting queued jobs
        queuedJobs = self.db.getJobsByStatus()
        #print 'JOBCONTROLLER.checkJobQueue',queuedJobs
        if len(queuedJobs) > 0:
            jobId = queuedJobs[0]
            if self.serverParams(jobId) is None:
                sp = self.runTask(jobId=jobId)
            else:
                self.runOnServer(jobId=jobId)

    def handleFinishedServerJob(self, jobId, status=None):
        #print 'handleFinishedServerJob',jobId
        sP = self.serverParams(jobId)
        if sP is None:
            return
        self.setServerParam(jobId, 'pollFinish', 0)
        self.setServerParam(jobId, 'pollReport', 0)
        if  sP.mechanism in ['qsub_shared', 'qsub_local', 'ssh_shared'] :
            status = self.loadRemoteRun(jobId, xmlDbFile=sP.dbXml)
        elif sP.mechanism != 'custom':
            if sP.mechanism != 'test':
                jobDir = os.path.join(sP.projectDirectory, 'CCP4_JOBS', 'job_' + sP.jobNumber)
                self.transportFiles(jobId, [[sP.local_finished_tarball, sP.remote_finished_tarball],
                                            [os.path.join(jobDir,'stdout.txt'), sP.remotePath+'stdout.txt'],
                                            [os.path.join(jobDir,'stderr.txt'), sP.remotePath+'stderr.txt']
                                            ], 'get', failSignal=False, finishHandler=self.handleFinishedServerJob1)
                return
            else:
                try:
                    shutil.copyfile(sP.remote_finished_tarball, sP.local_finished_tarball)
                except Exception as e:
                    print('ERROR: Copying tarball for local test:', sP.remote_finished_tarball,'to', sP.local_finished_tarball)
                    print(e)
                else:
                    if self._diagnostic:
                        print('Copied tarball for local test:', sP.remote_finished_tarball, 'to', sP.local_finished_tarball)
            status = self.loadRemoteRun(jobId, sP.local_finished_tarball)
        else:
            try:
                self.customHandler(jobId).handleFinishedJob(jobId, status=status)
            except Exception as e:
                self.handleServerFail(jobId=jobId, exception=e)
            else:
                status = self.loadRemoteRun(jobId, sP.local_finished_tarball)
        if status:
            self.deleteJobBalls(jobId)
        self.deleteServerParams(jobId)

    def handleFinishedServerJob1(self,jobId):
        #print 'handleFinishedServerJob1',jobId
        sP = self.serverParams(jobId)
        status = self.loadRemoteRun(jobId, sP.local_finished_tarball)
        if status:
            self.deleteJobBalls(jobId)
        self.deleteServerParams(jobId)

    def deleteJobBalls(self,jobId,dbXml=False):
        sP = self.serverParams(jobId)
        if dbXml:
            try:
                os.remove(str(sP.dbXml))
            except:
                print('ERROR deleting ',sP.dbXml)
        else:
            try:
                os.remove(str(sP.local_finished_tarball))
            except:
                print('ERROR deleting ',sP.local_finished_tarball)
            try:
                os.remove(str(sP.local_tarball))
            except:
                print('ERROR deleting ',sP.local_tarball)

    def watchLogFiles(self):
        # Initialise the runningJobs list is necessary
        '''
        This is for watching all jobs - probably wrong approach
        if self.runningJobs is None:
          runningJobs = self.db.getJobsByStatus(status=CCP4DbApi.JOB_STATUS_RUNNING)
          self.runningJobs = {}
          for jobId in runningJobs:
            logFile = self.db._makeJobFileName(jobId=jobId,mode='LOG')
            size = os.path.getsize(logFile)
            self.runningJobs[jobId] = { 'logFile' : logFile, 'size' : size }
        '''
        for jobId, details in list(self._watchedJobs.items()):
            size = self.getFileSize(details['logFile'])
            #print 'JOBCONTROLLER.watchLogFiles',jobId,details['logFile'],size
            if size is not None and (details['size'] is None or size > details['size']):
                #print 'watchLogFiles logUpdated',jobId
                self._watchedJobs[jobId]['size'] = size
                self.logFileUpdated.emit(jobId)
    
    def getFileSize(self, fileName):
        try:
            sz = os.path.getsize(fileName)
            size = int(sz.strip('L'))
        except:
            size = None
        return size

    def watchJob(self, jobId=None, mode=True):
        if mode:
            if jobId not in self._watchedJobs:
                from core import CCP4Modules
                status = self.db.getJobInfo(jobId=jobId, mode='status')
                if not status in ['Queued', 'Running']:
                    raise CException(self.__class__, 111, 'Job id' + str(jobId))
                logFile = CCP4Modules.PROJECTSMANAGER().makeFileName(jobId=jobId, mode='LOG')
                size = self.getFileSize(logFile)
                self._watchedJobs[jobId] = {'logFile' : logFile, 'size' : size}
        else:
            if jobId in self._watchedJobs:
                del self._watchedJobs[jobId]
        #print 'JOBCONTROLLER.watchJob',self._watchedJobs

    def getArgList(self, fileName, ccp4Dir=None, ccp4i2Dir=None):
        if ccp4Dir is not None:
            # This is only used for 'remote' jobs which will only be run on Linux so can assume '/' separator
            runTask = "${CCP4I2}/bin/runTask.py"
            exe = ccp4Dir + '/bin/ccp4-python'
        else:
            # running locally
            ccp4i2Dir = CCP4Utils.getCCP4I2Dir()
            runTask = os.path.join(ccp4i2Dir, 'bin', 'runTask.py')
            # For some reason, ccp4-python fails to start qprocesses on windows,so fall back to sys.executable...something fishy here: MN
            if sys.platform.startswith('win'):
                exe = sys.executable
            else:
                exe = CCP4Utils.pythonExecutable()
        argList = [exe, runTask, fileName]
        if self.configFile is not None:
            argList.extend(['-config', self.configFile])
        if self._dbFile is not None:
            argList.extend(['-db', self._dbFile])
        return argList

    def runTask(self, jobId=None, wait=None):
        from core import CCP4Modules
        from qtcore import CCP4HTTPServerThread
        from dbapi import CCP4DbApi
        #print 'CJobController.runTask pythonExecutable',CCP4Utils.pythonExecutable()
        controlFile = self.db._makeJobFileName(jobId=jobId, mode='JOB_INPUT')
        argList = self.getArgList(controlFile)
        db = CCP4Modules.PROJECTSMANAGER().db()
        taskname = db.getJobInfo(jobId=jobId)['taskname']
        if taskname == "moorhen_rebuild":
            argList = argList + ["--graphical"]
        if self._diagnostic:
            print('JOBCONTROLLER starting job:', argList , taskname)
        callDict = {}
        path, name = os.path.split(controlFile)
        if not self.USE_QPROCESS or taskname == "dui" or taskname == "dials_image" \
                                                       or taskname == "dials_rlattice":
            callDict['stdout'] = open(os.path.join(path, 'stdout.txt'), 'w')
            callDict['stderr'] = open(os.path.join(path, 'stderr.txt'), 'w')
        #MN run tasks with modified environment carrying :
        #-port of HTTP Server
        #-projectid and
        #-projectname
        try:
            httpServerPort = CCP4Modules.HTTPSERVER(fileName=self._dbFile).port
        except:
            httpServerPort = CCP4HTTPServerThread.DEFAULT_PORT
        my_env = os.environ.copy()
        my_env['CCP4I2_HTTP_PORT'] = str(httpServerPort)
        if jobId is not None:
            jobInfo = db.getJobInfo(jobId=jobId, mode=['projectid', 'projectname', 'jobnumber'])
            #print '\*\* jobInfo'
            #print jobInfo
            my_env['CCP4I2_PROJECT_ID'] = jobInfo['projectid']
            my_env['CCP4I2_PROJECT_NAME'] = jobInfo['projectname']
        callDict['env'] = my_env
        if not self.USE_QPROCESS or taskname == "dui" or taskname == "dials_image" \
                                                       or taskname == "dials_rlattice":
            print("RUNNING DUI (or DIALS viewers) IN SUBPROCESS PLUGIN MODE")   # KJS : tmp diag print.
            try:
                sp = subprocess.Popen(argList, **callDict)
            except Exception as e:
                # Failed to start job --
                sp = None
                self._errorReport.append(self.__class__, 102, 'Control file: ' + str(controlFile))
            self.db.updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_RUNNING)
            self.db.updateJob(jobId, 'processId', int(sp.pid))
            return sp
        else:
            qArgList = []
            for item in argList[1:]:
                qArgList.append(item)
            p = QtCore.QProcess(self)
            p.setObjectName(str(jobId))
            stdoutFile = os.path.join(path, 'stdout.txt')
            stderrFile = os.path.join(path, 'stderr.txt')
            if stdoutFile is not None:
                if os.path.exists(stdoutFile):
                    CCP4Utils.backupFile(stdoutFile, delete=True)
                p.setStandardOutputFile(stdoutFile)
            if stderrFile is not None:
                if os.path.exists(stderrFile):
                    CCP4Utils.backupFile(stderrFile, delete=True)
                p.setStandardErrorFile(stderrFile)
            #MN Changed here to apply environment edits to inform Coot or other tasks of how to talk to the http server
            processEnvironment = QtCore.QProcessEnvironment.systemEnvironment()
            for editItem in [('CCP4I2_HTTP_PORT', str(httpServerPort)), ('CCP4I2_PROJECT_ID', jobInfo['projectid']), ('CCP4I2_PROJECT_NAME', jobInfo['projectname'])]:
                processEnvironment.insert(editItem[0],editItem[1])
            # Fudge for MRBUMP task on OS X, because of rosetta requiring cctbx to set DYLD_LIBRARY_PATH.
            dbtn = CCP4Modules.PROJECTSMANAGER().db()
            jobInfoTN = dbtn.getJobInfo(jobId=jobId,mode=['taskname'])
            if jobInfoTN == "mrbump_basic" and sys.platform == "darwin":
                if 'CCP4' in os.environ:
                    print("Copying $CCP4/lib to DYLD_LIBRARY_PATH")
                    if 'DYLD_LIBRARY_PATH' in os.environ:
                        processEnvironment.insert('DYLD_LIBRARY_PATH', os.path.join(os.environ['CCP4'], 'lib') + os.pathsep + os.environ['DYLD_LIBRARY_PATH'])
                    else:
                        processEnvironment.insert('DYLD_LIBRARY_PATH', os.path.join(os.environ['CCP4'], 'lib'))
            #MN end edit
            p.finished.connect(lambda exitCode,exitStatus: self.handleFinish(jobId,exitCode,exitStatus))
            p.setProgram(argList[0])
            p.setArguments(qArgList)
            p.setProcessEnvironment(processEnvironment)
            p.startDetached()
            if wait is not None:
                p.waitForFinished(wait)
            self.db.updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_RUNNING)
            self.db.updateJob(jobId, 'processId', int(p.pid()))
            #print 'runTask processId',p.pid(),type(p.pid())
            return p

    def saveSh(self,jobId,argList,local_sh,pidfile=None):
        #pyi2 = '$CCP4/share/ccp4i2/bin/pyi2'
        pyi2 = '$CCP4/bin/ccp4-python'
        command = pyi2+" "+' '.join(argList[1:])
        #print 'saveSh',command
        tx = '''#!/bin/sh
CCP4='''+self.getServerParam(jobId,'ccp4Dir')+'''
export CCP4
.  ${CCP4}/bin/ccp4.setup-sh
if test -d $CCP4/share/ccp4i2; then
  CCP4I2=$CCP4/share/ccp4i2
elif test -d $CCP4/Frameworks/Python.framework/Versions/Current/lib/python3.7/site-packages/ccp4i2; then
  CCP4I2=$CCP4/Frameworks/Python.framework/Versions/Current/lib/python3.7/site-packages/ccp4i2
elif test -d $CCP4/Frameworks/Python.framework/Versions/Current/lib/python3.9/site-packages/ccp4i2; then
  CCP4I2=$CCP4/Frameworks/Python.framework/Versions/Current/lib/python3.9/site-packages/ccp4i2
elif test -d $CCP4/lib/python3.7/site-packages/ccp4i2; then
  CCP4I2=$CCP4/lib/python3.7/site-packages/ccp4i2
elif test -d $CCP4/lib/python3.9/site-packages/ccp4i2; then
  CCP4I2=$CCP4/lib/python3.9/site-packages/ccp4i2
elif test -d $CCP4/lib/site-packages/ccp4i2; then
  CCP4I2=$CCP4/lib/site-packages/ccp4i2
else
  CCP4I2=`python -c "import sysconfig; print(sysconfig.get_path('purelib'))"`/ccp4i2
fi
echo "Done CCP4 setup"
''' + command + '\n'
        #print 'saveSh',pidfile
        if pidfile is not None:
            tx = tx[0:-1] + \
''' &
pid=$!
/bin/cat <<EOM > ''' + pidfile + '''
$pid
EOM
echo "PID=$pid"
    '''
        CCP4Utils.saveFile(local_sh, tx )
        #print 'saved local sh file',local_sh

    def qsubOptionsFile(self,jobId=None):
        # Look for qsub option file in installation, user .CCP4I2/configs
        # job directory could also be searched
        for filename in [ os.path.join(CCP4Utils.getDotDirectory(),'configs','qsub_options'),
                          os.path.join(CCP4Utils.getCCP4I2Dir(),'local_setup','qsub_options') ]:
            if os.path.exists(filename):
                return filename
        return None

    def slurmOptionsFile(self,jobId=None):
        # Look for Slurm option file in installation, user .CCP4I2/configs
        # job directory could also be searched
        for filename in [ os.path.join(CCP4Utils.getDotDirectory(),'configs','slurm_options'),
                          os.path.join(CCP4Utils.getCCP4I2Dir(),'local_setup','slurm_options') ]:
            if os.path.exists(filename):
                return filename
        return None

    def runOnServer(self,jobId=None):
        '''Run ssh or qsub (SGE) or sbatch (Slurm) with or without shared file system
          always use a temporary db on remote machine to avoid sqlite/NFS issues'''
        from dbapi import CCP4DbApi
        from core import CCP4Modules
        self.db.updateJobStatus(jobId=jobId,status=CCP4DbApi.JOB_STATUS_REMOTE)
        sP = self.serverParams(jobId)
        #print 'runOnServer',sP.machine,sP.mechanism,sP.username
        controlFile = self.db._makeJobFileName(jobId=jobId,mode='JOB_INPUT')
        local_sh =  os.path.join(os.path.split(controlFile)[0],'remote.sh')

        if sP.mechanism in [ 'ssh_shared', 'qsub_local', 'qsub_shared' ]:
            ccp4i2Dir = os.environ.get('CCP4I2_REMOTE','$CCP4/share/ccp4i2')
            #print 'runOnServer ccp4i2Dir',ccp4i2Dir
            #argList = self.getArgList(controlFile,ccp4Dir='$CCP4',ccp4i2Dir='$CCP4/share/ccp4i2-devel')
            argList = self.getArgList(controlFile,ccp4Dir='$CCP4',ccp4i2Dir=ccp4i2Dir)
            # argList add the -dbxml arg for running on shared file system
            argList.extend ( [ '-dbxml' , sP.dbXml ] )
            # Add redirect stderr/stdout  
            argList.extend ( [ '>',os.path.join(os.path.split(controlFile)[0],'stdout.txt'),
                           '2>', os.path.join(os.path.split(controlFile)[0],'stderr.txt')])
            if self._diagnostic: print('JOBCONTROLLER starting remote job:',argList)
            self.saveSh(jobId,argList=argList,local_sh=local_sh,pidfile=self.pidFile(jobId))
#This is surely not necessary and indeed bad on Python3.
            #sP.local_report_file = os.path.join( sP.projectDirectory,'CCP4_JOBS','job_'+sP.jobNumber,'program.xml')
            sP.remoteSh = local_sh
        elif sP.mechanism != 'custom':
            # Not share disks - need to create and transport jobball
            # Use same remote dir and base of filename for all files
            # Beware remotePath setup assuming server is linux
            self.setServerParam(jobId,'remotePath',re.sub(r'\$USER',sP.username,sP.tempDir)+'/'+sP.projectName +'_'+sP.jobNumber+'_'+str(jobId)[0:4]+'_')
            #sP.local_tarball = os.path.join( sP.projectDirectory,'CCP4_TMP','job_'+sP.jobNumber+'_setup.ccp4db.zip')
            remote_tarball = sP.remotePath+'setup.ccp4db.zip'
            #sP.local_finished_tarball=os.path.join( sP.projectDirectory,'CCP4_TMP','job_'+sP.jobNumber+'_finished.ccp4db.zip')
            #sP.remote_finished_tarball=sP.remotePath+'finished.ccp4db.zip'
            #sP.remote_report_file = sP.remotePath+'work/project/CCP4_JOBS/job_'+sP.jobNumber+'/program.xml'
            #print 'remote_report_file',sP.remote_report_file
            #sP.local_report_file = os.path.join( sP.projectDirectory,'CCP4_JOBS','job_'+sP.jobNumber,'program.xml')
            ccp4i2Dir = os.environ.get('CCP4I2_REMOTE','$CCP4/share/ccp4i2')
            #print 'runOnServer ccp4i2Dir',ccp4i2Dir
            argList = self.getArgList(remote_tarball,ccp4Dir='$CCP4',ccp4i2Dir=ccp4i2Dir)   
            argList.extend ( [ '>',sP.remotePath+'stdout.txt',
                           '2>', sP.remotePath+'stderr.txt' ])
            if self._diagnostic: print('JOBCONTROLLER remote command:',argList)
            self.saveSh(jobId,argList=argList,local_sh=local_sh,pidfile=self.pidFile(jobId))
            sP.remoteSh = sP.remotePath+'remote.sh'
            if sP.mechanism not in ['test']:
                try:
                    self.transportFiles(jobId,[[local_sh,sP.remoteSh],[str(sP.local_tarball),remote_tarball]],'put',self.runOnServer2)
                except CException as e:
                    self.handleServerFail(jobId=jobId,exception=e)
                    return
                else:
                    return
            else:
                #  copy tarball and remote.sh to local tmp for local test
                import shutil
                for localFile,remoteFile in [[local_sh,sP.remoteSh],[str(sP.local_tarball),remote_tarball]]:
                    try:
                        shutil.copyfile(localFile,remoteFile)
                    except Exception as e:
                        print('ERROR copying file for local test',localFile,'to',remoteFile)
                    else:
                        print('Copy for local test:',localFile,'to',remoteFile)
        else:
            sP.local_tarball = os.path.join( sP.projectDirectory,'CCP4_TMP','job_'+sP.jobNumber+'_setup.ccp4db.zip')
            sP.local_report_file = os.path.join( sP.projectDirectory,'CCP4_JOBS','job_'+sP.jobNumber,'program.xml')
            sP.local_finished_tarball=os.path.join( sP.projectDirectory,'CCP4_TMP','job_'+sP.jobNumber+'_finished.ccp4db.zip')
            try:
                self.customHandler(jobId).setup(jobId)
            except Exception as e:
                self.handleServerFail(jobId=jobId,exception=e)
        self.runOnServer2(jobId)
    
    def runOnServer2(self,jobId):
        sP = self.serverParams(jobId)
        if sP is None: return
        self._transportFilesFinished(jobId)
        from dbapi import CCP4DbApi
        # Use paramiko to get ssh connection - beware failing and not resetting job status
        if sP.mechanism in [ 'ssh_shared', 'ssh' ]:
            try:
                client = self.runInSSH(jobId,sP.remoteSh)
            except CException as e:
                self.handleServerFail(jobId=jobId,exception=e)
            else:
                self.db.updateJobStatus(jobId=jobId,status=CCP4DbApi.JOB_STATUS_REMOTE)
        elif sP.mechanism in [ 'test' ]:
            try:
                self.runLocalTest(jobId,sP.remoteSh)
            except CException as e:
                self.handleServerFail(jobId=jobId,exception=e)
            else:
                self.db.updateJobStatus(jobId=jobId,status=CCP4DbApi.JOB_STATUS_REMOTE)
        elif sP.mechanism in [ 'qsub_local','qsub_shared','qsub_remote' ]:
            # Fire off qsub
            optionsFile = self.qsubOptionsFile()
            try:
                self.runInQsub(jobId,sP.remoteSh,optionsFile=optionsFile,mechanism=sP.mechanism)
                #print 'after runInQsub serverParams',sP
            except Exception as e:
                self.handleServerFail(jobId=jobId,exception=e)
            else:
                self.db.updateJobStatus(jobId=jobId,status=CCP4DbApi.JOB_STATUS_REMOTE)
        elif sP.mechanism in [ 'slurm_remote' ]:
            # Fire off sbatch
            optionsFile = self.slurmOptionsFile()
            try:
                self.runInSlurm(jobId,sP.remoteSh,optionsFile=optionsFile,mechanism=sP.mechanism)
                #print 'after runInSlurm serverParams',sP
            except Exception as e:
                self.handleServerFail(jobId=jobId,exception=e)
            else:
                self.db.updateJobStatus(jobId=jobId,status=CCP4DbApi.JOB_STATUS_REMOTE)
        elif sP.mechanism in [ 'custom' ]:
            try:
                self.customHandler(jobId).openConnection(jobId)
            except Exception as e:
                print('openConnection fail',e)
                self.handleServerFail(jobId=jobId,exception=e)
            else:
                try:
                    self.customHandler(jobId).submit(jobId)
                except Exception as e:
                    self.handleServerFail(jobId=jobId,exception=e)
                else:
                    print('runOnServer2 setServerParam pollReport')
                    self.setServerParam(jobId,'pollFinish',3)
                    self.setServerParam(jobId,'pollReport',True)
                    self.db.updateJobStatus(jobId=jobId,status=CCP4DbApi.JOB_STATUS_REMOTE)

    def handleServerFail(self,jobId,exception=None):
        self.handleFail(jobId,exception)
    
    @QtCore.Slot(tuple)
    def handleFailTuple(self,jobId_exception):
        jobId,exception = jobId_exception
        self.handleFail(jobId,exception)

    def handleFail(self,jobId=None,exception=None):
        from dbapi import CCP4DbApi
        if jobId is None:
            return
        self.db.updateJobStatus(jobId=jobId,status=CCP4DbApi.JOB_STATUS_PENDING)
        sP = self.serverParams(jobId)
        projectId = copy.deepcopy(sP.projectId)
        self.deleteServerParams(jobId)
        if not isinstance(exception,CErrorReport):
            err = CErrorReport(self.__class__,115,details = str(exception))
        else:
            err = exception
        print('Remote run fail:',str(exception))
        traceback.print_exc()
        self.serverJobFailed.emit((jobId,projectId,err))
        
    @QtCore.Slot(str,int,int)
    def handleFinish(self,jobId,exitCode=0,exitStatus=0):
        #print 'CJobController.handleFinish',jobId,status,exitStatus
        #if (exitStatus == QtCore.QProcess.CrashExit or status > 0) and jobId is not None:
        #SJM 15/1/2019 - removed the QProcess.CrashExit test as it seems to be random and cause unneccessary problems with bucref.
        if (exitCode > 0) and jobId is not None:
            from dbapi import CCP4DbApi
            self.db.updateJobStatus(jobId,CCP4DbApi.JOB_STATUS_FAILED)
        #from dbapi import CCP4DbUtils
        #CCP4DbUtils.makeJobBackup(jobId=jobId,db=self.db)

    def killJobProcess(self,jobId):
        err = CErrorReport()
        if self.serverParams(jobId) is not None:
            self.killRemoteJob(jobId)
            self.deleteServerParams(jobId)
            return
        pid = self.db.getJobInfo(jobId,'processId')
        #print 'CJobController.killJobProcess',jobId,pid
        if pid is not None:
            import signal
            try:
                self.killChildProcesses(pid)
            except Exception as e:
                err.append(self.__class__,114,details=str(e))
            if sys.platform.startswith('win'):
                try:
                    """
                    os.kill(pid,signal.CTRL_C_EVENT)
                    """
                    os.popen("taskkill /pid "+str(pid)+" /F")
                except Exception as e:
                    err.append(self.__class__,112,details=str(e))
            else:
                # SJM 23/09/2014 - SIGQUIT generates coredump which on OS X causes
                #                  unwanted stack trace window to appear. Is there
                #                  a reason to use QUIT instead of INT (Ctrl-C)?
                #                  os.kill(process,signal.SIGQUIT)
                try:
                    os.kill(pid,signal.SIGINT)
                except Exception as e:
                    err.append(self.__class__,112,details=str(e))
            from dbapi import CCP4DbApi
            #print 'to updateJobStatus',jobId,type(jobId)
            self.db.updateJobStatus(jobId,CCP4DbApi.JOB_STATUS_FAILED)
        else:
            err.append(self.__class__,113)
        if err.maxSeverity()>SEVERITY_WARNING:
            print('Error killing job process',err.report())
        return err

    def killChildProcesses(self,parent_pid):
        import signal
        import psutil
        # BEWARE this is picking up some other psutil!!!
        # See https://github.com/giampaolo/psutil
        # http://stackoverflow.com/questions/3332043/obtaining-pid-of-child-process
        if sys.platform == 'win32':
            sig = signal.CTRL_C_EVENT
        else:
            sig = signal.SIGQUIT
        try:
            proc = psutil.Process(parent_pid)
        except:
            print('killChildProcesses no child process found',parent_pid)
            return
        #print 'killChildProcesses proc',proc,type(proc),dir(proc)
        child_proc = proc.children(recursive=True)
        print('killChildProcesses child_proc',child_proc)
        for proc in child_proc:
            try:
                if sys.platform.startswith('win'):
                    os.popen("taskkill /pid "+str(proc.pid)+" /F")
                else:
                    os.kill(proc.pid, sig)
            except Exception as e:
                print('Failed to kill process',proc)
                print(e)

    def Exit(self):
        self.timer.stop()
        #sys.__stdout__.write('CJobController blockExit '+str(self.blockExit)+'\n');sys.__stdout__.flush()
        CCP4JobServer.CJobServer.Exit(self)

    def loadRemoteRun(self,jobId=None,compressedFile=None,xmlDbFile=None):
        from dbapi import CCP4DbApi
        from core import CCP4Modules
        #print 'loadRemoteRun',jobId,compressedFile,xmlDbFile
        # Extract database xmlfile
        if compressedFile is not None:
            #try:
            xmlDbFile = CCP4Modules.PROJECTSMANAGER().extractDatabaseXml(compressedFile)
            '''
            except CException as e:
              e.warningMessage('Reload exported job','Failed extracting database XML file from compressed file')
              return False
            except Exception as e:
              print 'Failed unpacking compressed file',str(e)
            '''
        else:
            xmlDbFile = os.path.join(os.path.split(xmlDbFile)[0],'DATABASE_final.db.xml')
    
        if not os.path.exists(xmlDbFile):
            try:
                CCP4Modules.PROJECTSMANAGER().updateJobStatus(jobId,CCP4DbApi.JOB_STATUS_FAILED)
            except:
                print('Failed updating remote job status for:',jobId,'Probably deleted job')
                pass
            return False

        try:
            projectId = CCP4Modules.PROJECTSMANAGER().db().getJobInfo(jobId=jobId,mode='projectid')
        except Exception as e:
            projectId = None
            print('Error in loadRemoteRun',jobId,compressedFile)
            print(e)
            return
        # Read dbxml file into  a CDbXml and check that it is for this project
        dbImport = CCP4DbApi.CDbXml(db=CCP4Modules.PROJECTSMANAGER().db(),xmlFile=xmlDbFile)
        importProjectInfo = dbImport.loadProjectInfo()
        if projectId is not None and dbImport.projectId != projectId:
            return
        dbImport.createTempTables()
        dbImport.loadTempTable()
        # If loading jobs to an existing project flag up jobs in temp tables that are already in db
        #dbImport._diagnostic = True
        dbImport.setExclInTempTables()

        # Extract job files from the compressed file
        if compressedFile is not None:
            projectDir =  CCP4Modules.PROJECTSMANAGER().db().getProjectInfo(projectId=projectId,mode='projectdirectory')
            from qtcore import CCP4Export
            # Despit the name this is not a separate thread!
            importThread = CCP4Export.ImportProjectThread(self,projectDir=projectDir,compressedFile=compressedFile)
            importThread.extractJobs(importThread.compressedFile,importThread.projectDir,dbImport=dbImport)

        # Update the database
        dbImport.cleanupTempTables()
        stats = dbImport.importStats()
        #print 'loadRemoteRun stats',stats
        if compressedFile is None:
            dbImport.importTempTables(includeCode=0)
        else:
            dbImport.importTempTables()
        dbImport.removeTempTables()
        #print 'loadRemoteRun DONE'

        dbImport.db.projectReset.emit({'projectId':dbImport.projectId})
        status = dbImport.db.getJobInfo(jobId,'status')
        dbImport.db.jobFinished.emit({'jobId':jobId,'projectId':dbImport.projectId,'status':status})
        CCP4Modules.PROJECTSMANAGER().backupDB()

    def updateReports(self):
        from qtgui import CCP4ProjectViewer
        currentlyOpenJobs = CCP4ProjectViewer.currentlyOpenJobs()
        #print 'updateReports',self.getJobsToPollReport(),'currentlyOpenJobs',currentlyOpenJobs
        
        for jobId in self.getJobsToPollReport():
            if jobId in currentlyOpenJobs:
                sP = self.serverParams(jobId)
                #print 'Polling for report update',jobId
                try:
                    if sP.mechanism in [ 'ssh_shared', 'qsub_local', 'qsub_shared' ]:
                        if not os.path.exists(sP.local_report_file):
                            #print 'updateReports program.xml not there'
                            pass
                        else:
                            # QFileSystemWatcher does not see changes made by another machine so must watch the
                            # local program.xml here
                            size = os.stat(sP.local_report_file).st_size
                            #print 'updateReports program.xml size',size
                            if not hasattr(sP,'local_report_size') or size>sP.local_report_size:
                                sP.local_report_size = size
                                self.remoteJobUpdated.emit((jobId,sP.local_report_file))
                    elif sP.mechanism in [ 'custom' ]:
                        self.customHandler(jobId).transportFiles(jobId,copyList= [[sP.local_report_file,sP.remote_report_file]] , mode='get')
                        self.remoteJobUpdated.emit((jobId,sP.local_report_file))
                    else:
                        #print 'CJobController.updateReports recover',jobId
                        try:
                            self.transportFiles(jobId,copyList= [[sP.local_report_file,sP.remote_report_file]] , mode='get', failSignal=False, diagnostic=False)
                        except Exception as e:
                            print('Failed copy program.xml')
                            pass
                        else:
                            pass
                except:
                    pass

    def listLocalProcesses(self,containsList = []):
        import psutil

        def contains(exe,containsList):
            if exe is None: return False
            for item in containsList:
                if exe.count(item): return True
            return False

        pInfoDict = {}
        me = CCP4Utils.getUserId()
        for proc in psutil.process_iter():
            try:
                pinfo = proc.as_dict(attrs=['pid', 'name', 'username','exe','create_time'])
            except psutil.NoSuchProcess:
                pass
            else:
                if pinfo['username'] == me and (len(containsList)==0 or contains(pinfo['exe'],containsList)):
                    try:
                        pinfo['parent'] = proc.parent().as_dict(attrs=['pid'])['pid']
                    except:
                        pinfo['parent'] = -1
                    pinfo['children'] = []
                    for p in proc.children(recursive=True):
                        try:
                            pinfo['children'].append(p.pid)
                        except:
                            pass
                    #print(pinfo)
                    pInfoDict[pinfo['pid']] = pinfo
        return pInfoDict


#===========================================================================================================
import unittest

class testController(unittest.TestCase):

    def setUp(self):
        self.controller = CJobController()

    def testRunFreerflag(self):
        from core import CCP4Container
        controlFile = os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','freerflag','test_data','test1.data.xml')
        stdoutFile = os.path.join(CCP4Utils.getTestTmpDir(),'CJobController_test.stdout')
        c = CCP4Container.CContainer()
        c.loadDataFromXml(controlFile)
        output = str(c.outputData.HKLOUT)
        print('Output file:',output,'Stdout:',stdoutFile)
        if os.path.exists(output): os.remove(output)
        self.controller.runTask(fileName=controlFile,wait=3000)
        self.assertTrue(os.path.exists(output),'No output file created')


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testController)
    return suite
    
def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
