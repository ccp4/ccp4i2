"""
Liz Potterton Apr 2016 - Separate 'remote' server code out from CCP4JobController
"""

import functools
import inspect
import os
import re
import socket
import subprocess
import sys

import UtilityThread
from lxml import etree
from PySide2 import QtCore
import paramiko

from . import CCP4Utils
from .CCP4ErrorHandling import CException


PARAMIKO_PORT=22

class CServerParams:
    SAVELIST = ['machine', 'username', 'ccp4Dir', 'tempDir', 'mechanism', 'keyFilename', 'serverProcessId',
                'remotePath', 'validate', 'customCodeFile', 'queueOptionsFile', 'sge_root', 'serverGroup']

    def __init__(self, **kw):
        self.jobId = kw.get('jobId', None)
        self.password = kw.get('password', None)
        self.jobNumber = None
        self.projectId = None
        self.projectName = None
        self.projectDirectory = None
        # pollFinish flags: 0=no poll (eg ssh expecting a signal) 1=poll for FINISHED file 2=poll qsub
        #                   3=custom poll (started in current i2 session) 4 = custom poll (from previous i2 session)
        self.pollFinish = 0
        self.pollReport = 0
        for item in CServerParams.SAVELIST:
            setattr(self, item, kw.get(item, None))
        self.serverProcessId = None

    def __str__(self):
        tx = ''
        for item in CServerParams.SAVELIST + ['jobNumber', 'projectName', 'projectDirectory', 'pollFinish', 'pollReport']:
            tx = tx + item + ': ' + str(getattr(self, item)) + '\n'
        return tx

    @property
    def remote_finished_tarball(self):
        return self.remotePath + 'finished.ccp4db.zip'

    @property
    def remote_report_file(self):
        if self.mechanism in ['ssh_shared', 'qsub_local', 'qsub_shared']:
            return None
        else:
            if self.remotePath is None:
                return 'work/project/CCP4_JOBS/job_' + self.jobNumber + '/program.xml'
            else:
                return self.remotePath + 'work/project/CCP4_JOBS/job_' + self.jobNumber + '/program.xml'

    @property
    def local_tarball(self):
        return os.path.join(self.projectDirectory, 'CCP4_TMP', 'job_' + self.jobNumber + '_setup.ccp4db.zip')

    @property
    def local_report_file(self):
        return os.path.join(self.projectDirectory, 'CCP4_JOBS', 'job_' + self.jobNumber, 'program.xml')

    @property
    def local_finished_tarball(self):
        return os.path.join(self.projectDirectory, 'CCP4_TMP', 'job_' + self.jobNumber + '_finished.ccp4db.zip')

    @property
    def dbXml(self):
        if self.mechanism not in ['ssh_shared', 'test', 'qsub_local', 'qsub_shared']:
            return os.path.join(self.projectDirectory, 'CCP4_TMP', 'DATABASE_job_' + self.jobNumber + '.db.xml')
        else:
            return os.path.join(self.projectDirectory, 'CCP4_JOBS', 'job_' + self.jobNumber, 'DATABASE.db.xml' )

    @property
    def finishedFile(self):
        if self.mechanism in ['ssh_shared', 'qsub_shared', 'qsub_local']:
            return os.path.join(self.projectDirectory, 'CCP4_JOBS', 'job_' + self.jobNumber, 'FINISHED')
        else:
            return self.remotePath[0:-1]+'.FINISHED'

    '''
    @property
    def jobNumber(self):
        """Get the jobNumber"""
        print 'jobNumber',self.jobId
        if self._jobNumber is None:
            jobInfo = self.db.getJobInfo(jobId=self.jobId,mode=['jobnumber','projectid'])
            self._jobNumber = jobInfo['jobnumber']
            self._projectId = jobInfo['projectid']
        return self._jobNumber

    @property
    def projectName(self):
        """Get the projectName"""
        if self._projectName is None:
            if self._projectId is None: self.jobNumber()
            projectInfo = self.db.getProjectInfo(projectId=self._projectId)
            self._projectName = projectInfo['projectname']
            self._projectDirectory = projectInfo['projectdirectory']
        return self._projectName

    @property
    def projectDirectory(self):
        """Get the projectDirectory"""
        print 'projectDirectory',self._projectId
        if self._projectDirectory is None:
            if self._projectId is None: self.jobNumber()
            projectInfo = self.db.getProjectInfo(projectId=self._projectId)
            self._projectName = projectInfo['projectname']
            self._projectDirectory = projectInfo['projectdirectory']
        return self._projectDirectory
    '''

    def getEtree(self):
        ele = etree.Element('serverParams')
        ele.set('jobId', str(self.jobId))
        for key in CServerParams.SAVELIST:
            if getattr(self, key, None) is not None:
                e = etree.Element(key)
                e.text = str(getattr(self,key))
                ele.append(e)
        return ele

    def setEtree(self,ele):
        self.jobId = ele.get('jobId')
        for e in ele.iterchildren():
            setattr(self, str(e.tag), str(e.text))
        self.setDbParams()

    def setDbParams(self):
        from .CCP4ProjectsManager import PROJECTSMANAGER
        jobInfo = PROJECTSMANAGER().db().getJobInfo(self.jobId, ['jobnumber', 'projectid', 'projectname'])
        self.jobNumber = jobInfo['jobnumber']
        self.projectId = jobInfo['projectid']
        self.projectName = jobInfo['projectname']
        self.projectDirectory = PROJECTSMANAGER().db().getProjectInfo(self.projectId, 'projectdirectory')


class CJobServer(QtCore.QObject):

    failedOpenConnection = QtCore.Signal(tuple)
    failedRemoteCommand = QtCore.Signal(tuple)
    testMessage = QtCore.Signal(str)
    remoteJobMessage = QtCore.Signal(str,str)
    testRemoteFilesSignal = QtCore.Signal(tuple)
    remoteProcessesList = QtCore.Signal(tuple)

    DB_KEYS = ['jobId', 'machine', 'username', 'mechanism', 'remotePath', 'customCodeFile', 'validate',
               'keyFilename', 'serverProcessId', 'serverGroup']
    ERROR_CODES = {301 : {'description' : 'Failed opening SSH connection'},
                   302 : {'description' : 'Failed opening SSH connection'},
                   303 : {'description' : 'Failed opening SSH connection'},
                   304 : {'description' : 'Failed opening SSH connection'},
                   320 : {'description' : 'Failed submitting job'},
                   330 : {'description' : 'Failed to open connection to copy files'},
                   331 : {'description' : 'Failed copying file to server'},
                   332 : {'description' : 'Failed copying file from server'},
                   340 : {'description' : 'Failed killing remote job - can not recover job id'},
                   350 : {'description' : 'Failed running job startup script on server'},
                   360 : {'description' : 'Failed to kill job - job id unknown'},
                   361 : {'description' : 'Failed to kill job - possible communication fail'},
                   362 : {'description' : 'Failed to find job status'},
                   363 : {'description' : 'Failed to recover remote process id'},
                   364 : {'description' : 'Failed to recover remote process status'},
                   365 : {'description' : 'Failed killing remote process - no information on the remote process'}}

    def __init__(self):
        self._diagnostic = False
        self._serverParams = {}
        self.runInSSHThreads = []
        # Lists of outstanding jobs when i2 restarts - check if they are finished and
        # then poll in CJobController.doChecks()
        paramiko.util.log_to_file(os.path.join(CCP4Utils.getDotDirectory(), 'status', 'paramiko.log'))

    def restoreRunningJobs(self):
        #self.loadParams()
        self.restoreFromDb()

    def getJobsToPollFinish(self, pollFinish=[]):
        ret = []
        for jobId,sP in list(self._serverParams.items()):
            if sP.pollFinish in pollFinish:
                ret.append(jobId)
        return ret

    def getJobsToPollReport(self):
        jobIdList = []
        for jobId,sP in list(self._serverParams.items()):
            if sP.pollReport > 0:
                jobIdList.append(jobId)
        return jobIdList

    def Exit(self):
        if len(self._serverParams) > 0:
            self.saveParams()
        elif os.path.exists(self.defaultParamsFileName()):
            os.remove(self.defaultParamsFileName())

    def serverParams(self, jobId):
        try:
            return self._serverParams[jobId]
        except:
            return None

    def saveParams(self, fileName=None):
        from . import CCP4File
        if fileName is None:
            fileName = self.defaultParamsFileName()
        fileObj = CCP4File.CI2XmlDataFile(fileName)
        fileObj.header.setCurrent()
        fileObj.header.function.set('JOBSERVERSTATUS')
        bodyEtree = etree.Element('serverParamsList')
        for jobId, sP in list(self._serverParams.items()):
            bodyEtree.append(sP.getEtree())
        fileObj.saveFile(bodyEtree)

    def defaultParamsFileName(self):
        return os.path.join(CCP4Utils.getDotDirectory(), 'status', 'jobServer.params.xml')

    def restoreFromDb(self):
        from .CCP4ProjectsManager import PROJECTSMANAGER
        jobInfoList = PROJECTSMANAGER().db().getServerJobs()
        #print 'restoreFromDb',jobInfoList
        for jobInfo in jobInfoList:
            sP = self._serverParams[jobInfo[0]] = CServerParams()
            ii = 0
            for key in self.DB_KEYS:
                setattr(sP, key, jobInfo[ii])
                ii += 1
            sP.setDbParams()
            self.setPollStatus(sP)

    def setPollStatus(self, sP):
        if sP.mechanism in ['qsub']:
            sP.pollFinish = 2
        elif sP.mechanism in ['qsub_remote', 'ssh', 'ssh_shared', 'slurm_remote']:
            if sP.validate in ['password', 'pass_key_filename']:
                sP.pollFinish = -1
            else:
                sP.pollFinish = 1
        elif sP.mechanism in ['custom']:
            sP.pollFinish = 4
        # And set to poll for report update
        if sP.remote_report_file is not None:
            sP.pollReport = 1

    def loadParams(self, fileName=None):
        from . import CCP4File
        if fileName is None:
            fileName = self.defaultParamsFileName()
        if not os.path.exists(fileName):
            return
        fileObj = CCP4File.CI2XmlDataFile(fileName)
        body = fileObj.getBodyEtree()
        #print 'CJobServer.loadParams body',body.tag
        # pollFinish flags: 0=no poll (eg ssh expecting a signal) 1=poll for FINISHED file 2=poll qsub
        #                   3=custom poll (started in current i2 session) 4 = custom poll (from previous i2 session)
        for ele in body.iterchildren():
            jobId = ele.get('jobId')
            if len(ele) > 0:
                sP = self._serverParams[jobId] = CServerParams()
                sP.setEtree(ele)
                sP.setDbParams()
                self.setPollStatus(sP)

    def createServerParams(self, jobId, params):
        self._serverParams[jobId] = params
        self._serverParams[jobId].jobId = jobId
        jobInfo = self.db.getJobInfo(jobId=jobId, mode=['jobnumber', 'projectid'])
        self._serverParams[jobId].jobNumber = jobInfo['jobnumber']
        self._serverParams[jobId].projectId = jobInfo['projectid']
        projectInfo = self.db.getProjectInfo(projectId=self._serverParams[jobId].projectId)
        self._serverParams[jobId].projectName = CCP4Utils.safeOneWord(projectInfo['projectname'])
        self._serverParams[jobId].projectDirectory = projectInfo['projectdirectory']
        #print 'setServerParams',jobId, self._serverParams[jobId]
        self.db.createServerJob(jobId, self._serverParams[jobId])

    def setServerParam(self, jobId, key, value):
#Trying to set the keys below is bad on Python3 and I guess not neccessary as surely property method should override?
        if jobId not in self._serverParams or key in ("remote_finished_tarball","remote_report_file","local_tarball","local_report_file","local_finished_tarball","dbXml","finishedFile"):
            pass
        else:
            setattr(self._serverParams[jobId] ,key, value)
            if key in self.DB_KEYS:
                self.db.updateServerJob(jobId, key, value)

    def getServerParam(self, jobId, key):
        from .CCP4ServerSetup import SERVERSETUP
        if jobId not in self._serverParams:
            return None
        rv = getattr(self._serverParams[jobId], key, None)
        if rv is not None:
            return rv
        if key in ['ccp4Dir']:
            #Try recovering info from the server setup data
            #This is in a def file (.CCP4I2/configs/serverSetup.def.xml) and loaded into
            # a container with each serverGroup labelled SERVERGROUPn
            serverGroup = getattr(self._serverParams[jobId], 'serverGroup', None)
            if serverGroup is None:
                return None
            if getattr(self, 'setupContainer', None) is None:
                self.setupContainer = SERVERSETUP()
                self.setupContainer.load()
            for dataObjName in self.setupContainer.dataOrder():
                if self.setupContainer.get(dataObjName).name == serverGroup:
                    value = self.setupContainer.get(dataObjName).get(key)
                    setattr(self._serverParams[jobId], key, value)
                    return value
        return None
 
    def deleteServerParams(self, jobId):
        try:
            self.db.deleteServerJob(jobId)
        except Exception as e:
            print('ERROR removing server job record from database\n',e)
        if jobId in self._serverParams:
            del self._serverParams[jobId]
        else:
            return
  
    def deleteServerParam(self, jobId, key):
        try:
            delattr(self._serverParams[jobId], key)
        except:
            pass
  
    def patchRemotePasswords(self, jobId,password):
        sP = self.serverParams(jobId)
        count = 0
        for jI in list(self._serverParams.keys()):
            if jI == jobId or (self._serverParams[jI].pollFinish < 0 and \
                               self._serverParams[jI].machine == sP.machine and self._serverParams[jI].username == sP.username):
                count += 1
                self._serverParams[jI].password = password
                self._serverParams[jI].pollFinish = 1
                self._serverParams[jI].pollReport = 1
        return count

    def checkForRemotePasswords(self, projectId):
        requirePass = []
        for jobId in self.getJobsToPollFinish([-1, -2, -3, -4]):
            sP = self.serverParams(jobId)
            if sP.pollFinish < 0 and sP.projectId == projectId:
                requirePass.append({'jobId' : jobId , 'machine' : sP.machine, 'username' : sP.username})
        return requirePass

    def getPKey(self,jobId):
        sP = self.serverParams(jobId)
        if sys.platform != "win32":
            keyFilename = re.sub(r'\$HOME', CCP4Utils.getHOME(), sP.keyFilename)
        else:
            keyFilename = sP.keyFilename
        if sP.validate == 'pass_key_filename': 
            pkey = paramiko.RSAKey.from_private_key_file(keyFilename, password=sP.password)
        else:
            pkey = paramiko.RSAKey.from_private_key_file(keyFilename)
        #print 'CJobServer.getPKey', keyFilename
        return pkey

    def pollQsubStat(self, mode='finished'):
        procCheck = subprocess.Popen(["qstat", "-xml"], stdout=subprocess.PIPE, universal_newlines=True, stderr=subprocess.STDOUT)
        outCheck,errCheck = procCheck.communicate()
        if self._diagnostic:
            print('pollQsubStat', outCheck, errCheck)
        if mode == 'finished':
            qsubIdDict= self.parseQsubStat(outCheck)
            for jobId in self.getJobsToPollFinish([2]):
                if self.serverParams(jobId).serverProcessId not in qsubIdDict:
                    self.handleFinishedServerJob(jobId)
        elif mode == 'status':
            return self.parseQsubStat(outCheck, mode)

    def handleFinishedServerJob(self, jobId):
        # Should be reimplemented in application
        # Do whatever necessary to retrieve and use job
        # Remember to delete from self.serverParams
        self.setServerParam(jobId, 'pollFinish', 0)
        self.setServerParam(jobId, 'pollReport', 0)

    def parseQsubStat(self, returnString, mode='finished'):
        parser = etree.XMLParser()
        tree = etree.fromstring(returnString, parser)
        if mode == 'finished':
            qsubIdDict = {}
            jobs = tree.xpath("job_info|queue_info")
            for job in jobs:
                for jobl in job.xpath("job_list"):
                    qsubIdDict[str(jobl.xpath("JB_job_number")[0].text)] = str(jobl.attrib["state"])
            #print 'parseQsubStat',qsubIdDict
            return qsubIdDict
        elif mode == 'status':
            return {}


    def pidFile(self, jobId, local=False):
        sP = self.serverParams(jobId)
        if local or sP.mechanism == 'ssh_shared':
            return  os.path.join(sP.projectDirectory, 'CCP4_JOBS', 'job_' + sP.jobNumber, 'pidfile.txt')
        elif sP.mechanism == 'ssh':
            return sP.remotePath + 'pidfile.txt'
        else:
            return None

    def transportFiles(self, jobId, copyList=[], mode='put', finishHandler=None, failSignal=True, diagnostic=True):
        sP = self.serverParams(jobId)
        #print  'transport files',mode,copyList
        if getattr(sP, 'sshThreadTransport', None) is not None:
            return False
        if finishHandler is None:
            finishHandler = self._transportFilesFinished
        """
        sP.sshThreadTransport = UtilityThread.UtilityThread(functools.partial(self._transportFiles, jobId, copyList, mode, failSignal, diagnostic))
        sP.sshThreadTransport.finished.connect(functools.partial(finishHandler, jobId))
        self.runInSSHThreads.append(sP.sshThreadTransport)
        sP.sshThreadTransport.start()
        """
        self._transportFiles( jobId, copyList, mode, failSignal, diagnostic)
        finishHandler(jobId)
        return True

    def _transportFiles(self, jobId, copyList=[], mode='put', failSignal=True, diagnostic=True):
        sP = self.serverParams(jobId)
        try:
            if len(sP.machine.split(":")) > 1:
                mach,port = sP.machine.split(":")
                port = int(port)
            else:
                mach = sP.machine
                port = PARAMIKO_PORT

            """
            transport = paramiko.Transport((mach, port))
            if sP.validate in ['key_filename', 'pass_key_filename'] and sP.keyFilename is not None and len(sP.keyFilename) > 0:
                transport.connect(username=sP.username, password=sP.password, pkey=self.getPKey(jobId))
            else:
                transport.connect(username=sP.username, password=sP.password)
            sftp = paramiko.SFTPClient.from_transport(transport)
            """
            ssh = paramiko.SSHClient()
            ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
            if sP.validate in ['key_filename', 'pass_key_filename'] and sP.keyFilename is not None and len(sP.keyFilename) > 0:
                ssh.connect(mach,username=sP.username, pkey=self.getPKey(jobId), auth_timeout=10, banner_timeout=10, timeout=10)
            else:
                ssh.connect(mach,username=sP.username, password=sP.password, auth_timeout=10, banner_timeout=10, timeout=10)
            sftp = ssh.open_sftp()
        except Exception as e:
            print('ERROR setting up paramiko FTP file transfer')
            print(e)
            err = CException(self.__class__, 330, 'Machine:' + str(self.getServerParam(jobId, 'machine')) + '\n' + str(e))
            if failSignal:
                self.failedOpenConnection.emit((jobId, err))
            raise err
        err = CException()
        for localName,remoteName in copyList:
            try:
                if mode == 'put':
                    sftp.put(localName, remoteName)
                    if diagnostic:
                        print("Copied", localName, 'to', remoteName)
                else:
                    sftp.get(remoteName, localName )
                    if diagnostic:
                        print("Copied", remoteName, 'to', localName)
            except Exception as e:
                if mode == 'put':
                    if diagnostic:
                        print('ERROR copying files', str(localName), 'to', str(remoteName))
                    err.append(self.__class__, 331, 'Local:' + str(localName) + ' Remote:' + str(remoteName) + '\n' + str(e))
                else:
                    if diagnostic:
                        print('ERROR copying files',str(localName),'from',str(remoteName))
                    err.append(self.__class__, 332, 'Local:' + str(localName) + ' Remote:' + str(remoteName) + '\n' + str(e))
        sftp.close()
        """
        transport.close()
        """
        if len(err) > 0:
            if failSignal:
                self.failedOpenConnection.emit((jobId, err))
            #raise err

    def _transportFilesFinished(self,jobId):
        #print '_transportFilesFinished',jobId
        self.deleteServerParam(jobId,'sshThreadTransport')

    def openSSHConnection(self, jobId, machine, username, password, keyFilename=None, timeout=None, maxTries=2, emitFail=True):
        client = paramiko.SSHClient()
        client.set_missing_host_key_policy(paramiko.AutoAddPolicy())
        connected= False
        ntries = 0
        #print 'openSSHConnection',machine,username,keyFilename,timeout,maxTries
        self.testMessage.emit('Connecting to ' + machine + ' with timeout ' + str(timeout))
        while ntries < maxTries and not connected:
            failCode = 0
            try:
                if len(machine.split(":")) > 1:
                    mach,port = machine.split(":")
                    port = int(port)
                else:
                    mach = machine
                    port = PARAMIKO_PORT
                if self._diagnostic:
                    print('TO client.connect', ntries)
                if keyFilename is not None and keyFilename != 'None' and len(keyFilename) > 0: 
                    client.connect(mach, port=port, username=username, password=password, key_filename=keyFilename, timeout=timeout)
                else:
                    client.connect(mach, port=port, username=username, password=password, timeout=timeout)
                #print 'DONE client.connect'
                connected=True
            except paramiko.ssh_exception.AuthenticationException:
                failCode = 1
                ntries = maxTries
                mess = 'Authentication failed'
            except paramiko.ssh_exception.BadAuthenticationType:
                failCode = 2
                mess = 'Authentication type not supported'
            except socket.gaierror:
                failCode = 3
                mess = 'Socket error'
            except Exception as e:
                failCode = 4
                mess = str(e)
                print(e)
            #print 'runOnServer try', ntries, 'failCode', failCode
            if not connected:
                self.testMessage.emit('   Try: ' + str(ntries) + ' failed: ' + str(mess))
            else:
                self.testMessage.emit('   Try: ' + str(ntries) + ' succeeded')
            ntries = ntries + 1
        if not connected:
            print("Some error submitting job")
            exc_type, exc_value, exc_tb = sys.exc_info()[:3]
            print("Error running job using ssh" + str(exc_type) + "\n" + str(exc_value))
            err = CException(self.__class__, 300 + failCode, 'Machine:' + machine + '\n' + str(exc_type) + "\n" + str(exc_value) + '\n')
            if emitFail:
                self.failedOpenConnection.emit((jobId, err))
            sys.stderr.write(mach)
            sys.stderr.write(port)
            raise err
        else:
            return client

    def testRemoteFiles(self, jobId=None, machine=None, username=None, password=None, keyFilename=None, remoteFileList=[], timeout=None, maxTries=2):
        if self.getServerParam(jobId, 'sshThreadRemoteFiles') is not None:
            return
        sP = self.serverParams(jobId)
        sP.sshThreadRemoteFiles = UtilityThread.UtilityThread(functools.partial(self._testRemoteFiles, jobId, machine, username, password,
                                                                                keyFilename, remoteFileList, timeout, maxTries))
        sP.sshThreadRemoteFiles.finished.connect(functools.partial(self._sendSSHFinished, jobId, None, 'sshThreadRemoteFiles'))
        self.runInSSHThreads.append(sP.sshThreadRemoteFiles)
        sP.sshThreadRemoteFiles.start()

    def _testRemoteFiles(self, jobId=None, machine=None, username=None, password=None, keyFilename=None, remoteFileList=[], timeout=None, maxTries=2):
        #print 'CJobServer._testRemoteFiles',jobId,machine,username,remoteFileList
        ret = []
        if isinstance(jobId, list):
            jI = jobId[0]
        else:
            jI = jobId
        try:
            client = self.openSSHConnection(jobId=jI, machine=machine, username=username, password=password,
                                            keyFilename=keyFilename, timeout=timeout, maxTries=maxTries)
        except Exception as e:
            print('testRemoteFiles fail\n',e)
            for item in remoteFileList:
                ret.append(False)
            return ret
        for remoteFile in remoteFileList:
            stdin, stdout, stderr = client.exec_command('[ -e ' + remoteFile + ' ] && echo "Found" || echo "Lost"')
            out = ""
            err = ""
            for line in stdout:
                out += line
            for line in stderr:
                err += line
            if self._diagnostic:
                print('testRemoteFiles', remoteFile, out, '*', err)
            if out.count('Found') > 0:
                self.testMessage.emit('   File: ' + str(remoteFile) + ' exists')
            else:
                self.testMessage.emit('   File: ' + str(remoteFile) + ' not found')
            ret.append(out.count('Found') > 0)
        client.close()
        self.deleteServerParam(jobId, 'sshThreadRemoteFiles')
        self.testRemoteFilesSignal.emit((jobId, ret))

    def runOnServer(self, jobId=None):
        '''Run ssh or qsub with or without shared file system
         reimplement in application'''
        pass

    def runInQsub(self, jobId, remoteSh, optionsFile=None, mechanism='qsub_local'):
        path,name = os.path.split(remoteSh)
        name = name[0:-9]
        stdout = os.path.join(path, name + 'qsub_stdout.txt')
        stderr = os.path.join(path, name + 'qsub_stderr.txt')
        sP = self.serverParams(jobId)
        comList = ["qsub", '-o', stdout, '-e', stderr, '-S', '/bin/sh', remoteSh]
        if optionsFile is not None:
            comList.insert(1, '@')
            comList.insert(2, optionsFile)
        if mechanism == 'qsub_local':
            try:
                proc = subprocess.Popen(comList, stdout=subprocess.PIPE, universal_newlines=True, stderr=subprocess.STDOUT)
                out,err = proc.communicate()
                if self._diagnostic:
                    print('qsub start', out, err)
            except:
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                print("Error running job using remote queue command 'qsub'\n\n" + str(exc_type) + "\n\n" + str(exc_value))
                raise CException(self.__class__, 320)
            else:
                try:
                    self.setServerParam(jobId, 'serverProcessId', out[:out.find('(')].rstrip().split()[-1])
                    self.setServerParam(jobId, 'pollFinish',1)
                    self.setServerParam(jobId, 'pollReport',1)
                except:
                    pass
        else:
            # Use ssh to submit qsub on a different machine
            if sP.sge_root is not None:
                comLine = 'SGE_ROOT=' + str(sP.sge_root) + ' ' + comList[0]
            else:
                comLine = comList[0]
            for item in comList[1:]:
                comLine = comLine + ' ' + item
            self.sendSSHCommand(jobId, comLine=comLine, retHandler=self.handleRemoteQsubStart)
            self.setServerParam(jobId, 'pollFinish', 1)
            self.setServerParam(jobId, 'pollReport', 1)

    def runInSlurm(self, jobId, remoteSh, optionsFile=None, mechanism='slurm_remote'):
        path,name = os.path.split(remoteSh)
        name = name[0:-9]
        stdout = os.path.join(path, name + 'slurm_stdout.txt')
        stderr = os.path.join(path, name + 'slurm_stderr.txt')
        sP = self.serverParams(jobId)
        comList = ["sbatch", '-o', stdout, '-e', stderr, remoteSh]
        if optionsFile is not None:
            comList.insert(1, '@')
            comList.insert(2, optionsFile)
        if mechanism == 'slurm_local':
            try:
                proc = subprocess.Popen(comList, stdout=subprocess.PIPE, universal_newlines=True, stderr=subprocess.STDOUT)
                out,err = proc.communicate()
                if self._diagnostic:
                    print('slurm start', out, err)
            except:
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                print("Error running job using command 'sbatch'\n\n" + str(exc_type) + "\n\n" + str(exc_value))
                raise CException(self.__class__, 320)
            else:
                try:
                    match = re.search(r'[0-9]+',out)
                    pid = -999
                    if match:
                        pid = match.group(0)
                    self.setServerParam(jobId, 'serverProcessId', pid)
                    self.setServerParam(jobId, 'pollFinish',1)
                    self.setServerParam(jobId, 'pollReport',1)
                except:
                    pass
        else:
            # Use ssh to submit to Slurm on a different machine
            comLine = comList[0]
            for item in comList[1:]:
                comLine = comLine + ' ' + item
            self.sendSSHCommand(jobId, comLine=comLine, retHandler=self.handleRemoteSbatchStart)
            self.setServerParam(jobId, 'pollFinish', 1)
            self.setServerParam(jobId, 'pollReport', 1)
    

    def handleRemoteQsubStart(self, jobId, out, err):
        self.setServerParam(jobId,'serverProcessId', out[:out.find('(')].rstrip().split()[-1])
        if self._diagnostic:
            print('handleRemoteQsubStart', self.serverParams(jobId).serverProcessId) 

    def handleRemoteSbatchStart(self, jobId, out, err):
        # parse the output from the sbatch command
        pid = -999
        match = re.search(r'[0-9]+',out)
        if match:
            pid = match.group(0)
        # TODO: add some error handling
        self.setServerParam(jobId,'serverProcessId', pid)
        if self._diagnostic:
            print('handleRemoteSbatchStart', self.serverParams(jobId).serverProcessId)

    def sendSSHCommand(self, jobId, comLine, finishHandler=None, retHandler=None, emitFail=True):
        if self._diagnostic:
            print('Remote server command:', comLine)
        sP = self.serverParams(jobId)
    
        sP.sshThread = UtilityThread.UtilityThread(functools.partial(self._sendSSHCommand, jobId, comLine, emitFail, retHandler))
        if finishHandler is not None:
            sP.sshThread.finished.connect(functools.partial(finishHandler, jobId))
        self.runInSSHThreads.append(sP.sshThread)
        sP.sshThread.start()
        #print 'runInSSH DONE'

    def _sendSSHCommand(self, jobId, comLine, emitFail=True, retHandler=None):
        sP = self.serverParams(jobId)
        if sP.keyFilename is not None and len(sP.keyFilename) > 0:
            if sP.keyFilename.count('HOME'):
                keyFilename = re.sub(r'\$HOME', CCP4Utils.getHOME(), sP.keyFilename)
            else:
                keyFilename = sP.keyFilename
        else:
            keyFilename = None
        sP.sshClient = self.openSSHConnection(jobId, sP.machine, sP.username, sP.password, keyFilename, emitFail=emitFail)
        stdin, stdout, stderr = sP.sshClient.exec_command(str(comLine))
        out = ""
        err = ""
        for line in stdout:
            out += line
        for line in stderr:
            err += line
        if self._diagnostic:
            print('Remote start stdout:', out)
            print('Remote start stderr:', err)
        if len(err) > 0:
            exc = CException(self.__class__, 350, err)
            if emitFail:
                self.failedRemoteCommand.emit((jobId, exc))
        if retHandler:
            try:
                retHandler(jobId, out, err)
            except Exception as e:
                print('ERROR handling return output from SSH command')
                print(e)
        return out, err

    @QtCore.Slot(str,str,str)
    def _sendSSHFinished(self, jobId, clientParam='sshClient', threadParam='sshThread'):
        #print "Finished remote job",jobId,threadParam
        #print '_runInSSHFinished',self.serverParams[jobId].mechanism
        sP = self.serverParams(jobId)
        if sP is None:
            return
        if clientParam is not None and hasattr(sP, clientParam):
            getattr(sP, clientParam).close()
            #getattr(sP, clientParam).__del__() - daft but how do we be show it is zapped?
        if hasattr(sP, threadParam):
            if getattr(sP, threadParam).retval is not None:
                stdout, stderr = getattr(sP,threadParam).retval
                if self._diagnostic:
                    print('Remote finish stdout:', stdout)
                    print('Remote finish stderr:', stderr)
            else:
                print('Failed analysing return from remote command')
        if hasattr(sP, threadParam):
            self.runInSSHThreads.remove(getattr(sP, threadParam))
            #getattr(sP,threadParam).__del__() -- this is nonsense for a thread

    def runInSSH(self, jobId, remoteSh):
        comLine = 'sh ' + remoteSh
        self.sendSSHCommand(jobId, comLine, finishHandler=self._sendSSHFinished, retHandler=self.handleRunInSSHDone)
        #retHandler=self.handleRunInSSHDone)
        self.setServerParam(jobId, 'pollReport', 1)
        self.setServerParam(jobId, 'pollFinish', 1)
        
    def _runInSSHFinished(self, jobId):
        self._sendSSHFinished(jobId)
        self.handleFinishedServerJob(jobId)

    def handleRunInSSHDone(self, jobId, out, err):
        if self._diagnostic: 
            print('handleRunInSSHDone', jobId, 'out', out, 'err', err)
        for line in out.split('\n'):
            if line.startswith('PID='):
                self.setServerParam(jobId, 'serverProcessId', int(line[4:]))
                break

    def handleRemoteSSHStatus(self, jobId, out, err):
        if self._diagnostic:
            print('handleRemoteSSHStatus', jobId, out, err)
        sP = self.serverParams(jobId)
        self.remoteJobMessage.emit(jobId, 'The process id for this job: ' + str(sP.serverProcessId), out, sP.machine)

    def checkLocalQsubStatus(self, jobId):
        try:      # KJS : Problem here. This func looks broken.
            proc = subprocess.Popen(comList='qstat', stdout=subprocess.PIPE, universal_newlines=True, stderr=subprocess.STDOUT)
            out,err = proc.communicate()
            if self._diagnostic:
                print('checkLocalQsubStatus', out, err)
        except:
            exc_type, exc_value,exc_tb = sys.exc_info()[:3]
            print("Error checking qsub status\n" + str(exc_type) + "\n" + str(exc_value))
            raise CException(self.__class__, 362)
        else:
            self.handleRemoteQsubStatus(jobId, out, err)

    def handleRemoteQsubStatus(self, jobId, out, err):
        # VU: It seems this function is not used at all.
        sP = self.serverParams(jobId)
        # VU: this check makes no sense, if the cluster runs other jobs, too
        if len(out) < 10:
            out = 'No jobs in queue (imples all finished/failed)'
        self.remoteJobMessage.emit(jobId, 'QSub status on ' + sP.machine, out, sP.machine)
        #print('handleRemoteQsubStatus: '+out)

    def getPidFileContent(self, jobId):
        sP = self.serverParams(jobId)
        if sP.mechanism == 'ssh':
            self.transportFiles(jobId,[[self.pidFile(jobId, local=True), self.pidFile(jobId)]], 'get', finishHandler=self.getPidFileContent2, failSignal=False)
        else:
            self.getPidFileContent2(jobId)

    @QtCore.Slot(str)
    def getPidFileContent2(self,jobId):
        self._transportFilesFinished(jobId)
        ret = CCP4Utils.readFile(self.pidFile(jobId, local=True))
        try:
            ret= int(ret.strip())
        except:
            ret = None
        else:
            self.setServerParam(jobId, 'serverProcessId', ret)
            #comLine = 'pgrep -P ' + str(self.serverParams[jobId].serverProcessId)+';ps -ef'
            comLine = 'ps -efu ' + self.serverParams(jobId).username
            self.sendSSHCommand(jobId, comLine=comLine, retHandler=self.handleRemoteSSHStatus, emitFail=False)
        return ret

    def killSSHJob(self, jobId):
        sP =  self.serverParams(jobId)
        '''
        if sP.serverProcessId is None: self.getPidFileContent(jobId)
        if sP.serverProcessId is None:
           raise CException(self.__class__,340)
        '''
        self.sendSSHCommand(jobId, comLine='pkill -F ' + self.pidFile(jobId), retHandler=self.handleRemoteDel)

    def killQsubJob(self, jobId):
        sP = self.serverParams(jobId)
        if sP.serverProcessId is None:
            raise CException(self.__class__, 360)
        comList = ['qdel', sP.serverProcessId]
        if sP.mechanism in ['qsub_local']:
            try:    # KJS Again a problem.
                proc = subprocess.Popen(comList, stdout=subprocess.PIPE,universal_newlines=True, stderr=subprocess.STDOUT)
                out,err = proc.communicate()
                print('kill qsub',out,err)
            except:
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                print("Error killing job\n"+str(exc_type)+"\n"+str(exc_value))
                raise CException(self.__class__,361)
            else:
                from ..dbapi import CCP4DbApi
                from .CCP4ProjectsManager import PROJECTSMANAGER
                self.deleteServerParams(jobId)
                PROJECTSMANAGER().db().updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_FAILED)
        elif sP.mechanism in ['qsub_remote']:
            comline = 'qdel ' + sP.serverProcessId
            self.sendSSHCommand(jobId, comLine='qdel ' + sP.serverProcessId, retHandler=self.handleRemoteDel)

    def killSlurmJob(self, jobId):
        sP = self.serverParams(jobId)
        if sP.serverProcessId is None:
            raise CException(self.__class__, 360)
        comList = ['scancel', sP.serverProcessId]
        if sP.mechanism in ['slurm_local']:
            try:
                proc = subprocess.Popen(comList, stdout=subprocess.PIPE,universal_newlines=True, stderr=subprocess.STDOUT)
                out,err = proc.communicate()
                print('kill Slurm job',out,err)
            except:
                exc_type, exc_value,exc_tb = sys.exc_info()[:3]
                print("Error killing job\n"+str(exc_type)+"\n"+str(exc_value))
                raise CException(self.__class__,361)
            else:
                from ..dbapi import CCP4DbApi
                from .CCP4ProjectsManager import PROJECTSMANAGER
                self.deleteServerParams(jobId)
                PROJECTSMANAGER().db().updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_FAILED)
        elif sP.mechanism in ['slurm_remote']:
            self.sendSSHCommand(jobId, comLine='scancel ' + sP.serverProcessId, retHandler=self.handleRemoteDel)

    def handleRemoteDel(self, jobId, out, err):
        from ..dbapi import CCP4DbApi
        from .CCP4ProjectsManager import PROJECTSMANAGER
        if self._diagnostic:
            print('handleRemoteQsubDel', jobId, out, err)
        self.deleteServerParams(jobId)
        PROJECTSMANAGER().db().updateJobStatus(jobId=jobId, status=CCP4DbApi.JOB_STATUS_FAILED)

    def killRemoteJob(self, jobId):
        sP = self.serverParams(jobId)
        if sP is None:
            raise CException(self.__class__, 365)
        if sP.mechanism in ['ssh', 'ssh_shared']:
            self.killSSHJob(jobId)
        elif sP.mechanism in ['qsub_local', 'qsub_shared', 'qsub_remote']:
            self.killQsubJob(jobId)
        elif sP.mechanism in ['slurm_remote']:
            self.killSlurmJob(jobId)
        elif sP.mechanism in ['custom']:
            self.customHandler(jobid).killJob(jobId)

    def runLocalTest(self, jobId, remoteSh):
        from .CCP4ProcessManager import PROCESSMANAGER
        PROCESSMANAGER().startProcess('sh', [remoteSh], handler=[self.localTestFinished, {'jobId':jobId}], ifAsync=True)

    def localTestFinished(self, pid, jobId=None):
        #print 'localTestFinished',pid,args
        self.handleFinishedServerJob(jobId)

    def pollForFinishedFlagFile(self):
        if getattr(self,'pollFinishConnected',None) is None:
            self.testRemoteFilesSignal.connect(self.pollForFinishedFlagFile1)
        pollByMachine = {}
        for jobId in self.getJobsToPollFinish([1]):
            sP = self.serverParams(jobId)
            if sP.mechanism in ['ssh', 'ssh_shared', 'qsub_remote', 'qsub_shared', 'qsub_local', 'slurm_remote']:
                if sP.mechanism in ['ssh_shared', 'qsub_shared', 'qsub_local']:
                    if os.path.exists(sP.finishedFile):
                        self.handleFinishedServerJob(jobId)
                else:
                    if sP.machine not in pollByMachine:
                        pollByMachine[sP.machine] = []
                    pollByMachine[sP.machine].append(jobId)
        #print 'pollForFinishedFlagFile pollByMachine',pollByMachine
        for machine, jobIdList in list(pollByMachine.items()):
            finfileList = []
            for jobId in jobIdList:
                sP = self.serverParams(jobId)
                if sP is not None:
                    finfileList.append (sP.finishedFile)
                #print 'pollForFinishedFlagFile finfileList',finfileList
            if len(finfileList) > 0:
                self.testRemoteFiles(jobId=jobId, machine=sP.machine, username=sP.username, password=sP.password, keyFilename=sP.keyFilename, remoteFileList=finfileList)
            
              
    @QtCore.Slot(tuple)
    def pollForFinishedFlagFile1(self, jobId_fileStatus):
        jobId, fileStatus = jobId_fileStatus
        if self._diagnostic:
            print('pollForFinishedFlagFile', jobId, fileStatus)
        if not isinstance(jobId,list):
            jobId = [jobId]
        for ii in range(len(jobId)):
            if fileStatus[ii]:
                self.handleFinishedServerJob(jobId[ii])

    def customHandler(self, jobId=None, customCodeFile=None):
        # This is creating one instance of the custom handler for each session
        # The alternative of one per job may be better for some cases???
        customCodeFile = str(customCodeFile)
        if not hasattr(self,'_customHandler'):
            self._customHandler = {}
        if jobId is not None:
            customCodeFile  = self.getServerParam(jobId, 'customCodeFile')
        if customCodeFile in self._customHandler:
            return self._customHandler[customCodeFile]
        if customCodeFile is None:
            print('ERROR attempting to access server custom code when no file provided')
            return None
        elif not os.path.exists(str(customCodeFile) ):
            print('ERROR server custom code file does not exist:' + customCodeFile)
            return None
        else:
            # If a startup.py file exists exec it first
            startupFile = os.path.join(os.path.split(customCodeFile)[0], 'startup.py')
            if os.path.exists(startupFile):
                try:
                    exec(compile(open(startupFile).read(), startupFile, 'exec'))
                except Exception as e:
                    print('ERROR execing custom server interface startup code:', startupFile)
                    print(e)
                    return None
            sys.path.append(os.path.split(customCodeFile)[0])
            module, err = CCP4Utils.importFileModule(customCodeFile)
            if err is not None:
                print('ERROR loading custom code file:' + customCodeFile)
                print(err)
                return None
            clsList = inspect.getmembers(module, inspect.isclass)
            if len(clsList) == 0:
                print('ERROR no classes found in custom code file:' + customCodeFile)
                return None
            if self._diagnostic:
                print('CServerParams.customHandler instantiating customHandler', clsList[0][0], self._serverParams)
            try:
                self._customHandler[customCodeFile] = clsList[0][1](self.serverParams)
            except Exception as e:
                print('ERROR instantiating custom handler class')
                print(e)
                return None
            return self._customHandler[customCodeFile]

    def listRemoteProcesses(self):
        runningJobs = self.getJobsToPollFinish(pollFinish=[0, 1])
        jobsByMachine = {}
        for jobId in runningJobs:
            sP = self.serverParams(jobId)
            if sP.machine not in jobsByMachine:
                jobsByMachine[sP.machine] = {}
            jobsByMachine[sP.machine][jobId] = sP.serverProcessId
        for machine, jobDict in list(jobsByMachine.items()):
            comLine=self.getServerParam(jobId, 'ccp4Dir') + '/bin/ccp4-python ' + \
                     self.getServerParam(jobId,'ccp4Dir') + '/share/ccp4i2-devel/bin/listProcesses.py'
            self.sendSSHCommand(list(jobDict.keys())[0], comLine=comLine,
                                retHandler=functools.partial(self.handleListRemoteProcesses, machine, jobDict), emitFail=False)
        # get qsub status
        runningJobs = self.getJobsToPollFinish(pollFinish=[2])

    def handleListRemoteProcesses(self, machine, jobDict, jobId, out, err):
        #print 'handleListRemoteProcesses out',out
        #print 'handleListRemoteProcesses err',err
        # Extract the processes and time from the remote job stdout
        for line in out.split('\n'):
            if line.startswith('processes') or line.startswith('atTime'):
                #print 'handleListRemoteProcesses line',line
                try:
                    exec(line)
                except Exception as e:
                    print('Failed evaling the return from listRemoteProcesses')
                    print('Error:', e)
        #print 'processes',processes
        #print 'atTime',atTime
        try:
            self.remoteProcessesList.emit((machine, jobDict, processes, atTime)) # KJS : Another basic error here by the looks of it.
        except:
            print('ERROR in handling listing of remote processes')
            print('execing line', line)
