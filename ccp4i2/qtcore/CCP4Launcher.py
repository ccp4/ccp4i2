"""
Copyright (C) 2011 University of York
Liz Potterton June 2011 - Class to launch other viewers
"""

import functools
import os
import re
import shutil
import socket
import sys

from PySide2 import QtCore, QtGui, QtWidgets

from ..core import CCP4Utils
from ..core.CCP4Modules import PREFERENCES, PROJECTSMANAGER, TASKMANAGER
from ..utils.QApp import QTAPPLICATION


class NamedSocket(socket.socket):
    def __init__(self,domain,protocol):
        socket.socket.__init__(self,domain,protocol)
        self._myName = 'Unknown'
        self._isActive = True

    def setActive(self,Active):
        self._isActive = Active

    def setName(self,name):
        self._myName = name

    def getName(self):
        return self._myName


def LAUNCHER():
    if CLauncher.insts is None:
        parent = QTAPPLICATION()
        CLauncher(parent)
    return CLauncher.insts

class CLauncher(QtCore.QObject):
    insts = None

    def __init__(self,parent=None):
        if parent is None:
            parent = QTAPPLICATION()
        QtCore.QObject.__init__(self,parent)
        if CLauncher.insts is None: CLauncher.insts = self
        self.hostPorts = {'CCP4mg' : {'hostname': 'localhost', 'port' : 9000}}
        self.sockets = {}
        self.launchTries = {}
        self.blockExit = False

    def Exit(self):
        pass
        #sys.__stdout__.write('CLauncher.Exit blockExit '+str(self.blockExit)+'\n');sys.__stdout__.flush()

    @QtCore.Slot(str,int,str,str)
    def openSocket(self, hostname, port, name='Unknown', command=None):
        #print 'openSocket',self.launchTries,hostname,port,name,command
        self.blockExit = True
        csocket = NamedSocket(socket.AF_INET, socket.SOCK_STREAM)
        try:
            csocket.connect((hostname,port))
            print("Successful connection to "+str(hostname)+':'+str(port))
            csocket.setName(name)
            self.sockets[name] = csocket
            if name in self.launchTries: del self.launchTries[name]
        except:
            print("Failed connection to "+str(hostname)+':'+str(port))  
            
            if name not in self.launchTries:
                print('openSocket launching:',name)
                p = self.launch(name)
                if p is not None:
                    self.launchTries[name] = 0
                else:
                    self.launchTries[name] = -1
            print('launchTries',self.launchTries)
                
            if self.launchTries[name]>=0:
                self.launchTries[name] = self.launchTries[name] + 1
                if self.launchTries[name]>10:
                    print('Failed to launch program',name)
                    self.launchTries[name] = -1          
                else:
                    timer = QtCore.QTimer(self)
                    timer.setSingleShot(True)
                    timer.timeout.connect(functools.partial(self.openSocket, hostname, port, name, command))
                    timer.start(1000)
        if name in self.sockets and command is not None:
            print('openSocket trying send')
            self.sockets[name].send(command)
        self.blockExit = False

    def getExecutable(self, viewer, guiParent=None):
        viewer = viewer.lower()
        path = None
        if viewer == 'ccp4mg':
            path = str(PREFERENCES().CCP4MG_EXECUTABLE)
            if  path is not None and os.path.isfile(path) and os.access(path, os.X_OK):
                return path
            if shutil.which(viewer) is not None:
                if sys.platform[0:3] == 'win':
                    whichExe = shutil.which(viewer)
                    if not viewer.lower().endswith(".bat") and whichExe.lower().endswith(".bat"):
                        viewer = whichExe
                return viewer
            if sys.platform[0:3] == 'win':
                #path = os.path.join(os.environ['CCP4'],'bin','winccp4mg.exe')
                # FIXME : This is a temporary hack whilst MG's ccp4mg.bat is broken.
                path = os.path.join(CCP4Utils.getCCP4I2Dir(),'bin','ccp4mg.bat')
            else:
                path = os.path.join(os.environ['CCP4'],'bin','ccp4mg')
            if  os.path.exists(path): return path
            if guiParent is not None:
                return self.queryExecutable(viewer=viewer,guiParent=guiParent)
        elif  viewer == 'coot':
            path = str(PREFERENCES().COOT_EXECUTABLE)
            #print 'CLauncher.getExecutable  prfrences coot',path
            if path is not None and os.path.isfile(path) and os.access(path, os.X_OK):
                if sys.platform == 'win32':
                    altpath = self.modifyCootBat()
                    if altpath is not None: return altpath
                return path
            if shutil.which(viewer) is not None: return viewer
            if guiParent is not None:
                return self.queryExecutable(viewer=viewer,guiParent=guiParent)
        elif  viewer == 'lidia':
            path = str(PREFERENCES().COOT_EXECUTABLE)
            #print 'CLauncher.getExecutable  prfrences coot',path
            if  path is not None and os.path.isfile(path) and os.access(path, os.X_OK):
                if sys.platform == 'win32':
                    altpath = self.modifyLidiaBat()
                    if altpath is not None: return altpath
                return path
            if shutil.which(viewer) is not None: return viewer
            if guiParent is not None:
                return self.queryExecutable(viewer=viewer,guiParent=guiParent)
        return path

    def launch(self,viewer=None,argList=[],envEdit=[],projectId=None,logFile=None,callBack=None,guiParent=None):
        #Bother - launch() gets called directly from CCP4ReportWidgets
        # Beware multiple calls to launch (from coot_model_building) seems to retain previous envEdit
        #print "CLauncher launch"
        editEnv = []
        editEnv.extend(envEdit)
        pwdSet = False
        for item in editEnv:
            if item[0]=='PWD': pwdSet = True
        if viewer in ['ccp4mg','coot','moorhen']:
            if not pwdSet and projectId is not None:
                cootD = self.cootDir(projectId)
                if cootD is not None:
                    editEnv.append(['PWD',self.cootDir(projectId)])
            editEnv.append( ['PYTHONHOME'] )
        # if viewer == 'moorhen':
        #     argList.insert(0,os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','coot_model_building','test_data','moorhen.py')))
        #     exe = os.path.join(CCP4Utils.getCCP4I2Dir(),'bin','Python')
        if viewer == 'loggraph':
            argList.insert(0,os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'pimple','MGQTmatplotlib.py')))
            if sys.platform == "win32" or sys.platform == "win64":
                exe = sys.executable
            else:
                exe = CCP4Utils.pythonExecutable()
        elif viewer == 'logview':
            exe = CCP4Utils.pythonExecutable()
            argList.insert(0,os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'smartie','qtrview.py')))
        elif viewer == 'PdbView':
            argList.insert(0,os.path.normpath(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','pdbview_edit','script','PdbViewView.py')))
            if sys.platform == "win32" or sys.platform == "win64":
                exe = sys.executable
            else:
                exe = CCP4Utils.pythonExecutable()
        elif viewer == 'viewhkl' and sys.platform == "darwin" and os.path.exists(os.path.join(os.environ["CCP4"],"ViewHKL.app")):
            exe = "open"
            argList.insert(0,"--args") 
            argList.insert(0,os.path.join(os.environ["CCP4"],"ViewHKL.app"))
            argList.insert(0,"-a")
        else:
            exe = self.getExecutable(viewer,guiParent=guiParent)
        if sys.platform[0:3] == "win" and viewer in ("ccp4mg", "lidia"):
            whichExe = shutil.which(exe)
            if not exe.lower().endswith(".bat") and whichExe.lower().endswith(".bat"):
                exe = whichExe
        if exe is None:
            print('Do not know executable path so calling as: ',viewer)
            exe = viewer
        if sys.platform[0:3] == "win" and viewer == "logview":
            exe += ".bat"

        qArgList = []
        for item in argList:
            if item is not None:
                qArgList.append(item)
        p = QtCore.QProcess(self)
        if len(editEnv) > 0:
            #if len(editEnv)>1 or editEnv[0][0] not in  ['pwd','PWD']:
            processEnvironment = QtCore.QProcessEnvironment.systemEnvironment()
            for editItem in editEnv:
                if len(editItem)==1:
                    processEnvironment.remove(editItem[0])
                else:
                    if editItem[0] in ['pwd','PWD']:
                        p.setWorkingDirectory(editItem[1])
                        print('1 Setting process working directory',editItem[1])
                    #else:
                    processEnvironment.insert(editItem[0],editItem[1])
            p.setProcessEnvironment(processEnvironment)
            #else:
            #  p.setWorkingDirectory(editEnv[0][1])
            #  print '2 Setting process working directory',editEnv[0][1]
        print('launch', viewer, type(viewer), exe, argList)
        if logFile is not None:
            p.setStandardOutputFile(logFile)
        if callBack is not None:
            p.finished.connect(callBack)
        if  exe is None:
            p.start(viewer,qArgList)
        else:
            #print('\n\nexe',exe,qArgList)
            p.start(exe,qArgList)
        return p

    def modifyCootBat(self):
        cootBat = PREFERENCES().COOT_EXECUTABLE.__str__()
        if not os.path.splitext(cootBat)[1] == '.bat' or not os.path.exists(cootBat):
            return None
        text = CCP4Utils.readFile(cootBat)
        text0 = re.sub('start .* coot-bin.exe', 'coot-bin.exe', text)
        modFile = os.path.join(CCP4Utils.getDotDirectory(), 'runwincoot.bat')
        CCP4Utils.saveFile(modFile, text0, overwrite=True)
        return modFile

    def modifyLidiaBat(self):
        cootBat = PREFERENCES().COOT_EXECUTABLE.__str__()
        if not os.path.splitext(cootBat)[1] == '.bat' or not os.path.exists(cootBat):
            return None
        text = CCP4Utils.readFile(cootBat)
        text0 = re.sub('start .* coot-bin.exe','lidia.exe',text)
        modFile = os.path.join(CCP4Utils.getDotDirectory(),'runlidia.bat')
        CCP4Utils.saveFile(modFile,text0,overwrite=True)
        return modFile

    def cootDir(self,projectId):
        if projectId is None:
            return None
        else:
            cootDir = os.path.join(PROJECTSMANAGER().getProjectDirectory(projectId=projectId),'CCP4_COOT')
            if not os.path.exists(cootDir):
                try:
                    os.mkdir(cootDir)
                except:
                    return None
        #print '  cootDir',projectId,cootDir
        return cootDir
    
    def openInViewer(self,viewer=None,fileName=None,jobId=None,projectId=None,style=None,guiParent=None,fileType='chemical/x-pdb'):
        if projectId is None and jobId is not None:
            try:
                projectId = PROJECTSMANAGER().db().getJobInfo(jobId,mode='projectid')
            except:
                pass

        if viewer == '':
            root, ext = os.path.splitext(str(fileName))
            if ext == '.mtz':
                viewer = 'viewhkl'

        #print 'LAUNCH.openInViewer viewer',viewer,jobId,projectId,fileName
        if viewer.lower() == 'ccp4mg':
            argList = self.mgComLine(fileName=fileName,jobId=jobId,style=style)
            self.launch(viewer='ccp4mg',argList=argList,projectId=projectId,guiParent=guiParent)
        elif viewer.lower() == 'pdbview':
            argList = [fileName]
            self.launch(viewer='PdbView',argList=argList,projectId=projectId,guiParent=guiParent)
        elif viewer.lower() == 'coot0':
            argList = self.cootComLine(fileName=fileName,jobId=jobId)
            self.launch(viewer='coot',argList=argList,projectId=projectId,guiParent=guiParent)
            #else:
            #  self.launch(viewer='moorhen',argList=argList)
        elif viewer.lower() == 'coot_job' or viewer.lower() == 'coot':
            self.runCootJob(contextJobId=jobId, projectId=projectId, fileName=fileName, fileType=fileType )
        elif viewer.lower() == 'lidia':
            cootExeDir = None
            if hasattr(PREFERENCES(), 'COOT_EXECUTABLE'):
                if os.path.isfile(str(PREFERENCES().COOT_EXECUTABLE)):
                    cootExeDir = str(PREFERENCES().COOT_EXECUTABLE)
            if cootExeDir is None:
                cootExeDir = shutil.which('coot')
            cootDir = os.path.normpath(os.path.dirname(os.path.dirname(cootExeDir)))
            envEdit = [['COOT_PREFIX', cootDir]]
            COOT_DATA_DIR = os.path.normpath(os.path.join(cootDir, 'share', 'coot'))
            envEdit.append(['COOT_DATA_DIR',COOT_DATA_DIR])
            self.launch(viewer='lidia', argList=[fileName], envEdit=envEdit, projectId=projectId, guiParent=guiParent)
        elif viewer.lower() == 'viewhkl':
            if not isinstance(fileName,(tuple,list)):
                fileNameList = [fileName]
            else:
                fileNameList = fileName
            self.launch('viewhkl',fileNameList)
        elif viewer.lower() == 'loggraph':
            print('openInViewer', fileName, jobId)
            if fileName is not None:
                self.launch(viewer='loggraph', argList=[fileName])
            elif jobId is not None:
                logfile = PROJECTSMANAGER().makeFileName(jobId=jobId, mode='REPORT')
                #print 'openInViewer loggraph',logfile
                self.launch(viewer='loggraph', argList=[logfile])
                #self.launch(viewer='loggraph')
        elif viewer.lower() == 'postscript':
            QtGui.QDesktopServices.openUrl(QtCore.QUrl.fromLocalFile(fileName))

    def mgComLine(self, fileName=None, jobId=None, style='Bonds:All_atoms'):
        # common styles: 'bonds:all_atoms' 'ribbons:colour_chains' 'ribbons:colour_blend_thru_chain' 'ribbons:secondary_structure'
        #print 'LAUNCHER.mgComLine jobId',jobId
        if style is None: style = 'Bonds:All_atoms'
        line = ['-norestore','-drawing_style',style]
        if fileName is None:
            pass
        elif isinstance(fileName, (tuple,list)):
            line.extend(fileName)
        else:
            line.append(fileName)
        if jobId is not None:
            scenefileList =PROJECTSMANAGER().getSceneFiles(jobId=jobId)
            if len(scenefileList) > 0:
                line.append(scenefileList[0])
            else:
                fileList = self.getFilesToDisplay(jobId)
                line.extend(fileList)
        return line

    def runCootJob(self,projectId=None,contextJobId=None,fileName=None,fileType='chemical/x-pdb'):
        from ..dbapi.CCP4DbApi import FILETYPELIST
        from ..dbapi.CCP4DbUtils import COpenJob
        if projectId is None:
            projectId = PROJECTSMANAGER().db().getJobInfo(contextJobId,'projectid')
        openJob = COpenJob(projectId=projectId)
        openJob.createJob(taskName='coot_rebuild',contextJobId=contextJobId)
        if fileName is not None:
            #need to establish what sort of file we are trying to view...difficult from here since we have only a fileName
            #and a projectId
            db = PROJECTSMANAGER().db()
            
            #SELECT Files.JobId,Files.FileID,ImportFiles.ImportId,Files.FileTypeId,Files.Filename
            
            selectedInfos = []
            fileInfos=db.getProjectFiles(projectId=projectId, topLevelOnly=False)
            selectedInfos = [anInfo for anInfo in fileInfos if db.getFullPath(anInfo[1]) == fileName]
            #print 'runCootJob selectedInfos',selectedInfos
            if len(selectedInfos) > 0:
                partFileInfo = selectedInfos[0]
                fileInfo = db.getFileInfo(fileId=partFileInfo[1],mode=['jobid','filename','relpath','projectid','annotation','filecontent','filesubtype'])
                dobj = None
                if FILETYPELIST[int(partFileInfo[3])][1] == 'chemical/x-pdb':
                    dobj = openJob.container.inputData.XYZIN_LIST
                elif FILETYPELIST[int(partFileInfo[3])][1] == 'application/CCP4-mtz-map' and fileInfo['filesubtype']== 1:
                    dobj = openJob.container.inputData.FPHIIN_LIST
                elif FILETYPELIST[int(partFileInfo[3])][1] == 'application/CCP4-mtz-map' and fileInfo['filesubtype']== 2:
                    dobj = openJob.container.inputData.DELFPHIIN_LIST
                if dobj is not None:
                    dobj.addItem()
                    dobj[-1].set( { 'baseName' : fileInfo['filename'],
                                  'relPath' :  fileInfo['relpath'],
                                  'project' : projectId,
                                  'annotation' : fileInfo['annotation'],
                                  'dbFileId' :selectedInfos[0][1],
                                  'contentFlag' :fileInfo['filecontent'],
                                  'subType' :fileInfo['filesubtype'] } )
            else:
                dobj = openJob.container.inputData.XYZIN_LIST
                dobj.addItem()
                dobj[-1].setFullPath(fileName)
        openJob.runJob()

    def cootComLine(self,fileName=None,jobId=None):
        line = ['--no-state-script']
        if fileName is not None:
            fileNameList = []
            if isinstance(fileName,(tuple,list)):
                fileNameList.extend(fileName)
            else:
                fileNameList.append(fileName)
            for f in  fileNameList:
                ext = os.path.splitext(f)[1]
                if ext == '.pdb':
                    line.extend(['--pdb',f])
                elif ext == '.mtz':
                    #line.extend(['--auto',f])
                    line.extend(['--data',f])
                elif ext == '.map':
                    line.extend(['--map',f])
            return line
        elif jobId is not None:
            paramsFile, comLine = self.makeCootParamsFile(contextJobId=jobId)
            #line.extend(['--i2params',paramsFile])
            line.extend(comLine)
            return line

    def makeCootParamsFile(self,contextJobId=None,projectId=None,taskName='coot_rebuild'):
        from ..core import CCP4Container
        from ..core import CCP4File
        from ..core import CCP4ModelData
        from ..core.CCP4XtalData import CMapCoeffsDataFile
        comLine = []
        # Make container
        defFile = TASKMANAGER().lookupDefFile(name=taskName)
        container = CCP4Container.CContainer()
        container.loadContentsFromXml(fileName=defFile)
        #Loop over inputData files to pull best file of type from database
        db = PROJECTSMANAGER().db()
        for key in container.inputData.dataOrder():
            dobj = container.inputData.get(key)
            if isinstance(dobj,CCP4File.CDataFile):
                fileIdList = db.getFileByJobContext(contextJobId=contextJobId,fileType=dobj.qualifiers('mimeTypeName'),
                                                    subType=dobj.qualifiers('requiredSubType'),contentFlag=dobj.qualifiers('requiredContentFlag'),
                                                    projectId=projectId)
                #print 'makeCootParamsFile',contextJobId,key,dobj.qualifiers('mimeTypeName'),dobj.qualifiers('requiredSubType'),fileIdList
                if len(fileIdList)>0:
                    fileInfo = db.getFileInfo(fileId=fileIdList[0], mode=['jobid', 'filename', 'relpath', 'projectid',
                                                                          'annotation', 'filecontent', 'filesubtype'])
                    if projectId is None: projectId = fileInfo['projectid']
                    dobj.set({'baseName' : fileInfo['filename'], 'relPath' : fileInfo['relpath'],
                              'project' : projectId, 'annotation' : fileInfo['annotation'], 'dbFileId' :fileIdList[0],
                              'contentFlag' :fileInfo['filecontent'], 'subType' :fileInfo['filesubtype'] } )
                if isinstance(dobj,CCP4ModelData.CPdbDataFile):
                    comLine.extend(['--pdb',dobj.__str__()])
                elif isinstance(dobj,CMapCoeffsDataFile):
                    comLine.extend(['--data',dobj.__str__()])
        paramsFileName = os.path.join(db.getProjectDirectory(jobId=contextJobId),'CCP4_COOT','params.xml')
        #print 'makeCootParamsFile paramsFileName',paramsFileName
        container.header.projectId = projectId
        container.saveDataToXml(  fileName=paramsFileName )
        return paramsFileName,comLine

    def getFilesToDisplay(self,jobId):
        from ..core import CCP4Data
        taskName =  PROJECTSMANAGER().db().getJobInfo(jobId,'taskname')
        paramList = TASKMANAGER().getTaskAttribute(taskName,'MGDISPLAYFILES',None)
        # No DISPLAYFILES defined for task so Just give the output files
        if paramList is None:
            fileNameList = PROJECTSMANAGER().db().getJobFiles(jobId=jobId,mode='fullPath')
        else:
            container = PROJECTSMANAGER().db().getParamsContainer(jobId)
            fileNameList = []
            for param in paramList:
                obj = container.find(param)
                #print 'LAUNCHER.getFilesToDisplay',param,repr(obj)
                if obj is not None:
                    if isinstance(obj,(list,CCP4Data.CList)):
                        for item in obj:
                            fileName = str(item)
                            if os.path.exists(fileName):
                                fileNameList.append(fileName)
                    else:
                        fileName = str(obj)
                        if os.path.exists(fileName):
                            fileNameList.append(fileName)
        #print 'LAUNCHER.getFilesToDisplay',fileNameList
        return fileNameList

    def makeCommand(self,viewer='CCP4mg',command=None,data=None):
        #print 'makeCommand',viewer,command
        if viewer == 'CCP4mg':
            if command == 'openFile':
#FIXME PYQT - potential hairiness
                t = data
                return "<begin>self.postEventSignal.emit(,{"+ \
                                "\"type\":"+str(int(QtCore.QEvent.FileOpen))+ \
                                ",\"file\":"+'"'+t+'"'+ \
                                "})<end>"
        return None

    def queryExecutable(self,viewer,guiParent):
        from ..qtgui.CCP4FileBrowser import CFileDialog
        self.findViewerDialog = CFileDialog(guiParent,'Find '+viewer,filters=[' (*)'],projectCombo=False)
        label = QtWidgets.QLabel("""Sorry - failed to find """+viewer+""". Please enter the program executable and then view again.\n The executable can also be set in Preferences.""",self.findViewerDialog)
        label.setStyleSheet( "QLabel { font-weight: bold;  border: 2px solid} ")
        self.findViewerDialog.addWidget(label)
        self.findViewerDialog.selectFile.connect(functools.partial(self.handleFindViewer,viewer))
        self.findViewerDialog.show()
        self.findViewerDialog.raise_()

    @QtCore.Slot(str,str)
    def handleFindViewer(self,viewer,filePath):
        #print 'handleFindViewer',viewer,filePath
        if os.path.exists(filePath):
            if viewer in ['coot','coot_job']:
                PREFERENCES().COOT_EXECUTABLE.set(filePath)
            else:
                PREFERENCES().CCP4MG_EXECUTABLE.set(filePath)
            PREFERENCES().save()
