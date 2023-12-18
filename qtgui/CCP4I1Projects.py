from __future__ import print_function


"""
     CCP4I1Projects.py: CCP4 GUI Project
     Copyright (C) 2016   STFC

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
   Liz Potterton Jan 2016 - Classes for importing and displaying CCP4i projects
"""

"""
To interpret i1 project def files we need to access
  ccp4i/etc/modules.def
  ccp4i/tasks/*.def
"""

import os,re,time,sys,traceback,copy,functools
from core import CCP4Utils
from qtgui import CCP4ProjectWidget, CCP4WebBrowser, CCP4WebView, CCP4TextViewer
from core import CCP4File
from core.CCP4Modules import *
from core.CCP4TaskManager import TASKMANAGER
from core.CCP4ErrorHandling import *
from PySide2 import QtCore,QtGui, QtWidgets
from lxml import etree


# CI1ProjectManager -> CI1Project -> CI1Job ->CI1File is hierarchy of
# data classes to load and hold data from i1 directories.def and database.def files
# They have appropriate  'get' methods to support the CI1ProjectWidget (ie the Qt tree view)
# I1PROJECTMANAGER is function to access the hierarchy


def CI1PREFERENCES():
  if CI1Preferences.insts is None:
      p = CI1Preferences()
  return CI1Preferences.insts

def isAlive(qobj):
    import shiboken2
    return shiboken2.isValid(qobj)


def splitDefLine(line):
  
  if len(line) == 0 or line[0] in ['#','@']:
    return None,'',None
  key,value = line.split(None,1)
  if value[0] == '_':
    return line.split(None,2)
  else:
    return key,'',value

def readI1DefFile(fileName,loadTypes=False,diagnostic=False):
  metaData = {}
  params = {}
  err = CErrorReport()

  try:
    text = CCP4Utils.readFile(fileName)
  except:
    err.append(CI1ProjectViewer,102,details=fileName,stack=False)
    return metaData,params,err

  for line in text.split('\n'):
    try:
      if line.startswith('#CCP4I'):   
        metaData[line.split()[1]] = line.split(None,2)[2]
      elif line.startswith('#'):
        pass
      elif line.startswith('@'):
        #Beware included files on the script def files
        #@ [FileJoin [GetEnvPath CCP4I_TOP] tasks harvest.def]
        importFile = 'UNKNOWN'
        importTask= 'UNKNOWN'
        try:
          importTask = line.replace(']','').split()[-1]
          importFile =  os.path.join(CCP4Utils.getCCP4Dir(),'share','ccp4i','tasks',importTask)
          importMetaData,importParams = readI1DefFile(importFile)
          params.update(importParams)
        except Exception as e:
          print('Failed loading subsiduary task def file',line)
          err.append(CI1ProjectViewer,110,details='For task: '+importTask+' file: '+importFile,stack=False)
        
      else:
        # Extract param name
        words = line.split()
        #print 'words',words
        if len(words)>=2:
          if words[0].count(','):
            pname,idx = words[0].split(',')
            idx = int(idx)
          else:
            pname = words[0]
            idx = -1
          # Extract param value
          resplit = re.search(r'(.*)"(.*)"(.*)',line)
          if resplit is not None:
            value = resplit.groups()[1]
          else:
            value = words[-1]
          # Create params item
          if pname not in params:
            if idx>=0:
              params[pname] = [ None , [] ]
            else:
              params[pname] = [ None , None ]
          # Set param type if possible
          if params[pname][0] is None and len(words) > 2 and words[1].startswith('_'):
            params[pname][0] = words[1]
          # Set param value - beware items may be missing from list
          if idx >= 0:
            while len(params[pname][1])<idx+1:
              params[pname][1].append(None)
            params[pname][1][idx] = value
          else:
            params[pname][1] = value

    
    except Exception as e:
      print('Error reading def file line:',line)
      print(e)
      err.append(CI1ProjectViewer,111,details=line,stack=False)
    
  if loadTypes and metaData.get('TASKNAME',None) is not None:
    taskDefFile = os.path.join(CCP4Utils.getCCP4Dir(),'share','ccp4i','tasks',metaData['TASKNAME']+'.def')
    try:
      taskMeta,taskParams,err0 = readI1DefFile(taskDefFile)
    except Exception as e:
      print('Failed reading task def file',taskDefFile)
      err.append(CI1ProjectViewer,112,details=taskDefFile,stack=False)
    else:
      for key in list(params.keys()):
        if params[key][0] is None and key in taskParams:
          params[key][0] = taskParams[key][0]

  if diagnostic:
    for key,value in list(metaData.items()):
      print(key,value)
    for key,value in list(params.items()):
      print(key,value)

  return metaData,params,err



  
class CI1TreeItemFolder(CCP4ProjectWidget.CTreeItemFolder):

  def __init__(self,parent=None,info={},tree=None):
    if tree is not None:
      info['name'] = tree.get('name')
    CCP4ProjectWidget.CTreeItemFolder.__init__(self,parent=parent,info=info)
    if tree is not None:
      folderList = tree.findall('folder')
      for fEle in folderList:
        fItem = CI1TreeItemFolder(self,tree=fEle)
        self.appendChildFolder(fItem)
      projectList = tree.findall('project')
      for pEle in projectList:
        pItem = CI1TreeItemProject(self,tree=pEle)
        self.appendChildProject(pItem)

  def getEtree(self):
    ele = ET.Element('folder')
    ele.set('name',self.name)
    for fItem in self.childFolders:
      ele.append(fItem.getEtree())
    for pItem in self.childProjects:
      ele.append(pItem.getEtree())
    return ele
    
  def mimeData(self):
    from lxml import etree
    root = ET.Element('name')
    root.text = self.name
    encodedData = QtCore.QByteArray()
    encodedData.append(etree.tostring(root,pretty_print=False))
    mimeData = QtCore.QMimeData()
    mimeData.setData('I1Folder',encodedData)
    return mimeData

      

class CI1TreeItemProject(CCP4ProjectWidget.CTreeItemProject):

  ERROR_CODES = { 101 : { 'description' : 'No job data extracted from the project database file' },
                  102 : { 'description' : 'Failed saving the updated project database file' },
                  103 : { 'description' : 'Failed updating the project database file' }
                    }
  def __init__(self,parent=None,infoList=[],directory=None,lastTaskTime=None,broken=0,tree=None):
    if tree is not None:
      name = tree.get('name')
      infoList = [name,name]
    CCP4ProjectWidget.CTreeItemProject.__init__(self,parent=parent,infoList=infoList)
    #print 'CI1TreeItemProject.__init__',self,tree,broken
    self.refProject = infoList[0]
    self.refDbIndex = None
    self.directory = directory
    if lastTaskTime is not None:
      self.lastTaskTime = CCP4ProjectWidget.formatData(lastTaskTime)
      self.machineTime = lastTaskTime
    elif tree is not None and tree.find('lastTaskTime') is not None:
      lastTaskTime = int(tree.find('lastTaskTime').text)
      self.lastTaskTime = CCP4ProjectWidget.formatData(lastTaskTime)
      self.machineTime = lastTaskTime      
    else:
#FIXME PYQT - or maybe None? This used to set QVariant.
      self.lastTaskTime = ""
      self.machineTime = None
    self.broken = broken
    self.annotation = None
    self.tagList = []
    if tree is not None:
      if tree.find('annotation') is not None:
        self.annotation = tree.find('annotation').text
      if tree.find('tagList') is not None:
        #print 'CI1TreeItemProject.__init__',tree.find('tagList').text
        for ele in tree.find('tagList').findall('tag'):
          self.tagList.append(str(ele.text))
        

  def setDirectory(self,directory):
    self.directory = directory
    if os.path.exists(directory):
      self.broken = 1
    else:
      self.broken = -1
    index = self.model().modelIndexFromProject(self.refProject)
    if index is not None:
      self.model().dataChanged.emit(index,index)

  def hasAnnotation(self):
    return (len(self.tagList)>0 or self.annotation is not None)
        

  def getEtree(self):
    #print 'CI1TreeItemProject.getEtree',self.annotation,self.tagList
    ele = ET.Element('project')
    ele.set('name',self.refProject)
    if self.annotation is not None:
      e = ET.Element('annotation')
      e.text = str(self.annotation)
      ele.append(e)
    if len(self.tagList)>0:
      eL = ET.Element('tagList')
      for item in self.tagList:
        e = ET.Element('tag')
        e.text = item
        eL.append(e)
      ele.append(eL)
    if self.machineTime is not None:
      e = ET.Element('lastTaskTime')
      e.text = str(self.machineTime )
      ele.append(e)
        
    return ele

  def loadDummy(self):
    jItem = CI1TreeItemJob(parent=self, info = { 'jobid' : '0',
             'jobnumber' : '0',
             'taskname' : 'No jobs in project',
             'jobtitle' : 'dummy',
             'evaluation' : 'Unknown',
             'status' :  'Failed',
             'finishtime' : time.time(),
             'parentjobid' : None
             })
    self.appendChildJob(jItem)
    
  def loadDatabase(self):
    #print 'loadDatabase',self.directory
    dbFile = os.path.normpath(os.path.join(self.directory,'CCP4_DATABASE','database.def'))
    print('Loading database file:',dbFile)
    metaData,params,err = readI1DefFile(dbFile)
    if err.maxSeverity()>SEVERITY_WARNING:
      print(err.report())
      return err
    if params.get('NJOBS',None) is None:
     return CErrorReport(self.__class__,101,details=dbFile,stack=False)
    self.nJobs =int( params['NJOBS'][1] )
    #print 'loadDatabase nJobs',self.nJobs
    lastTime = 0
    self.childJobs = []
    for nJ in range(1,self.nJobs+1):
      if params['TASKNAME'][1][nJ] is None:
        pass
      else:
        lastTime = max(lastTime,CCP4Utils.safeFloat(params['DATE'][1][nJ],0))
        info = { 'jobid' : str(nJ),
             'jobnumber' : str(nJ),
             'taskname' : params['TASKNAME'][1][nJ],
             'jobtitle' : params['TITLE'][1][nJ],
             'evaluation' : 'Unknown',
             'status' :  CI1TreeItemJob.STATUSCONV.get(params['STATUS'][1][nJ],'Pending'),
             'finishtime' : CCP4Utils.safeFloat(params['DATE'][1][nJ]),
             'parentjobid' : None
             }
        jItem = CI1TreeItemJob(self,info=info)
        self.appendChildJob(jItem)
        # Load the job def file to try to sort out file types
        #jobFile = os.path.join(dirName,str(nJ)+'_'+params['TASKNAME'][nJ]+'.def')
        #jobMetaData,jobParams,err = readI1DefFile(jobFile)
        inputFiles = self.extractFileList(params['INPUT_FILES'][1][nJ],params['INPUT_FILES_DIR'][1][nJ],jItem,input=1)
        outputFiles = self.extractFileList(params['OUTPUT_FILES'][1][nJ],params['OUTPUT_FILES_DIR'][1][nJ],jItem)
        #print 'loadDatabase outputFiles',self.jobs[-1].outputFiles
    if lastTime>0:  self.lastTaskTime = CCP4ProjectWidget.formatDate(lastTime)
    self.broken = 2
    return CErrorReport()

  def editDatabaseDef(self,resetDir=True,movedFiles={}):
    err = CErrorReport()
    dbFile = os.path.normpath(os.path.join(self.directory,'CCP4_DATABASE','database.def'))
    #print 'editDatabaseDef dbFile',dbFile
    try:
      lineList = CCP4Utils.readFile(dbFile).split('\n')
    except:
      err.append(self.__class__,101,dbFile)
      return err
    done = False

    if resetDir:
      for ii in range(len(lineList)):
        if resetDir and lineList[ii].startswith('#CCP4I PROJECT'):
          words = lineList[ii].split(None,3)
          lineList[ii] = '#CCP4I PROJECT '+str(self.refDbIndex)+'  '+str(self.directory)
          done = True
          break

    if len(movedFiles)>0:
      for ii in range(len(lineList)):
        if lineList[ii].startswith('INPUT_FILES'):
          nJ = lineList[ii].split()[0][12:]
          if movedFiles.get(nJ,None) is not None:
            text = 'INPUT_FILES,'+nJ+'              "'+movedFiles[nJ][0]
            for f in movedFiles[nJ]: text = text + ' '+f
            lineList[ii] = text+'"'
  
    if done:
      CCP4Utils.backupFile(dbFile,delete=True)
      text = ''
      for line in lineList: text = text + line + '\n'
      #print 'editDatabaseDef',text
      try:
        CCP4Utils.saveFile(dbFile,text=text)
      except:
        err.append(self.__class__,102,dbFile)
    else:
      err.append(self.__class__,103)
    return err
        
  def extractFileList(self,fStr,dStr,parentJob,input=0):
    fileObjList = []
    if fStr is not None:
      fList = fStr.split()
    else:
      fList = []
    if dStr is not None:
      dList = dStr.split()
    else:
      dList = []
    #print 'extractFileList',type(fList),fList,type(dList),dList
    for nF in range(min(len(fList),len(dList))):
      f = fList[nF].strip('"')
      d = dList[nF].strip('"')
      if d == 'FULL_PATH':
        import_ = True
      else:
        import_= None

      ext = os.path.splitext(f)[1]
      if ext == '.pdb':
        fType = 2
      elif ext == '.mtz':
        fType = 4
      else:
        fType = 0

      broken = 0
      if d == 'FULL_PATH':
        fileName = f
      elif d == self.refProject:
        fileName = f
        if os.path.exists(os.path.join(self.directory,f)):
          broken = 1
        else:
          broken = -1
      else:
        pObj = self.model().getProject(d)
        if pObj is not None:
          fileName = os.path.join(pObj.directory,f)
        else:
          fileName = d+'/'+f
          broken = -1
      if broken == 0: 
        if fileName is None:
          self.broken = -1
        elif not os.path.exists(fileName):
          self.broken = -1
        else:
         self.broken = 1
         
      infoList = [ '0','0',import_,fType , fileName, fileName ]
      fItem = CI1TreeItemFile(parentJob,infoList,broken=broken,input=input)
      if input:
        parentJob.appendChildInputFile(fItem)
      else:
        parentJob.appendChildOutputFile(fItem)
    return fileObjList
  
  def data(self,column,role):
    if role == QtCore.Qt.DisplayRole:
      if column == 0:
        return self.projectName
      elif column == 1 and self.lastTaskTime is not None:
        return self.lastTaskTime
    elif role == QtCore.Qt.DecorationRole:
      if column == 0:
        return CCP4ProjectWidget.jobIcon('ccp4')
    elif role == QtCore.Qt.ToolTipRole:
      if column == 0:
        text = 'Directory: ' + str(self.directory)
        if self.annotation is not None:
          text = text + '\nAnnotation: '+self.annotation
        if len(self.tagList)>0:
          text = text + '\nTags: '+ self.tagList[0]
          for t in self.tagList[1:]: text = text + ', ' + t
        return text
    elif role == QtCore.Qt.UserRole:
      if column == 0:
        return self.identifier
    elif role == QtCore.Qt.FontRole:
      return CCP4ProjectWidget.CTreeItemJob.boldFont
    elif role ==  QtCore.Qt.ForegroundRole:
      if self.broken<0:
        return QtGui.QBrush(QtCore.Qt.red)
    else:
#FIXME PYQT - or maybe None? This used to return QVariant.
      return None

  
  def canFetchMore(self):
    #print 'CI1TreeItemProject.canFetchMore',self.broken,self.broken in [0,1]
    return (self.broken in [0,1])

  def fetchMore(self):
    self.loadDatabase()
    

  def mimeData(self):
    #print 'CI1TreeItemProject.mimeData'
    from lxml import etree
    root = ET.Element('name')
    root.text = self.getProjectName()
    encodedData = QtCore.QByteArray()
    encodedData.append(etree.tostring(root,pretty_print=False))
    mimeData = QtCore.QMimeData()
    mimeData.setData('I1Project',encodedData)
    return mimeData

class CI1TreeItemJob(CCP4ProjectWidget.CTreeItemJob):
  STATUSCONV = { 'FINISHED' : 'Finished',
                 'FAILED' : 'Failed' }
      
  def __init__(self,parent=None,info={}):
    CCP4ProjectWidget.CTreeItemJob. __init__(self,parent=parent,info=info)
    self.finished = info.get('status','blah') == 'Finished'
    
  def data(self,column,role):
    if role == QtCore.Qt.DisplayRole:
      if column == 0:
        return self.name
      elif column == 1:
        return self.dateTime
    elif role == QtCore.Qt.UserRole:
      if column == 0:
        return self.identifier
    elif role == QtCore.Qt.DecorationRole:
      if column == 0:
        return self.status
#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  def setName(self,jobTitle=None):
    try:
      self.name = self.jobNumber +' '+ TASKMANAGER().getI1TaskTitle(self.taskName)
    except:
      self.name = self.jobNumber
      print('Error in CI1TreeItemJob.setName',self.jobNumber,'*',self.taskName,'*',TASKMANAGER().getI1TaskTitle(self.taskName))

  def logFile(self):
    log = os.path.join(self.parent().directory,str(self.jobNumber)+'_'+str(self.taskName)+'.log')
    if os.path.exists(log+'.html'):
        return log+'.html'
    elif os.path.exists(log):
        return log
    else:
        return None


class CI1TreeItemFile(CCP4ProjectWidget.CTreeItemFile):

  def __init__(self,parent=None,infoList=[],displayColumn=0,broken=0,input=0):
    CCP4ProjectWidget.CTreeItemFile.__init__(self,parent=parent,infoList=infoList,displayColumn=displayColumn)
    self.broken = broken
    self.input = input
    
  def data(self,column,role):
    if role == QtCore.Qt.DisplayRole:
      if column == 0:
        return self.name
    elif role == QtCore.Qt.UserRole:
      if column == 0:
        return self.identifier
    elif role == QtCore.Qt.DecorationRole:
      if column == 0:
        return self.icon
    elif role ==  QtCore.Qt.ForegroundRole:
      if self.broken<0:
        return QtGui.QBrush(QtCore.Qt.red)
    elif role ==  QtCore.Qt.FontRole:
      if column == 0 and self.input>0:
        return CCP4ProjectWidget.CTreeItemJob.italicFont
#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  
  def setName(self,annotation,mimeType=None,maxChar=40,displayJobNumber=False):
    #print 'CI1TreeItemFile.setName',annotation
    self.name = annotation
  
  def filePath(self):
    projDir = self.parent().parent().directory
    #print 'CI1TreeItemFile.filePath',projDir
    f = os.path.join(projDir,self.fileName)
    if os.path.exists(f):
      return f
    else:
      return None
  
  def mimeData(self):
    fN = self.filePath()
    if fN is None: return None
    urlList = [QtCore.QUrl()]
    urlList[0].setPath( fN )
    mimeData = QtCore.QMimeData()
    mimeData.setUrls(urlList)
    #print 'CI1TreeItemFile.mimeData',fN
    return mimeData

    
class CI1ProjectModel(QtCore.QAbstractItemModel):
  COLUMNS = ['name','date']
  ERROR_CODES = { 101 : { 'description' : 'Failed reading directories.def file'},
                  102 : { 'description' : 'Failed to find and update directory in directories.def file'},
                  103 : {  'description' : 'Failed writing directories.def file'}
                  }
                  
  def __init__(self):
    QtCore.QAbstractItemModel.__init__(self)
    self.rootItem=CI1TreeItemFolder(self,{ 'name' : 'root' } )
    self.loadErrorReport = CErrorReport()
    self.initialised = False
    self.sourceFile = None

  def resetAll(self,args):
    self.beginResetModel()
    self.rootItem =  CCP4ProjectWidget.CTreeItemProject(self,['-1','root']) 
# ???
    #self.setupModel()
    self.endResetModel()

  def loadSupplement(self,fileName=None):
    if fileName is None:
      fileName = os.path.join(CCP4Utils.getDotDirectory(),'i1supplement','temporary.xml')
    if os.path.exists(fileName):
      f = CCP4File.CI2XmlDataFile(fullPath=fileName)
      root = f.getBodyEtree()
      if root.find('folderList') is not None:
        for fEle in root.find('folderList').findall('folder'):
          fItem = CI1TreeItemFolder(self.rootItem,tree=fEle)
          self.rootItem.appendChildFolder(fItem)
      if root.find('projectList') is not None:
        #print 'loadSupplement',root.find('projectList').findall('project')
        for pEle in root.find('projectList').findall('project'):
          pItem = CI1TreeItemProject(self.rootItem,tree=pEle)
          self.rootItem.appendChildProject(pItem)
        
  def saveStatus(self):
    self.saveSupplement()
    CI1PREFERENCES().save()
    
  def saveSupplement(self,fileName=None):
    if fileName is None:
      if self.sourceFile is not None and os.access(os.path.split(self.sourceFile)[0], os.W_OK|os.X_OK):
        fileName = os.path.join(os.path.split(self.sourceFile)[0],'i2supplement.xml')
      else:
        fileName = os.path.join(CCP4Utils.getDotDirectory(),'i1supplement','temporary.xml')
    #print 'CI1ProjectModel.saveSupplement',fileName
    f = CCP4File.CI2XmlDataFile(fullPath=fileName)
    if not os.path.exists(fileName):
      f.header.setCurrent()
      f.header.function.set('I1SUPPLEMENT')
      f.header.comment = self.sourceFile
    body = ET.Element('body')
    fEleList = ET.Element('folderList')
    body.append(fEleList)
    if hasattr(self,"root") and hasattr(self.root,"childFolders"):
        for fItem in self.rootItem.childFolders:
            fEleList.append(fItem.getEtree())
    pEleList = ET.Element('projectList')
    for pItem in self.rootItem.childProjects:
      if pItem.hasAnnotation(): pEleList.append( pItem.getEtree())
    if len(pEleList)>0: body.append(pEleList)
    f.saveFile(bodyEtree=body)
    
  def loadDirectoriesDefFile(self,fileName=None):
    err = CErrorReport()
    metaData,params,err = readI1DefFile(fileName)
    if err.maxSeverity()>SEVERITY_WARNING or 'N_PROJECTS' not in params:
      return CErrorReport(self.__class__,101,details = fileName,stack=False)
    nProjects = CCP4Utils.safeInt(params.get('N_PROJECTS')[1])
    for nP in range(1,nProjects+1):
      pObj = self.getProject(params['PROJECT_ALIAS'][1][nP])
      #print 'loadDirectoriesDefFile',params['PROJECT_ALIAS'][1][nP],pObj
      if pObj is None:
        infoList= [params['PROJECT_ALIAS'][1][nP],params['PROJECT_ALIAS'][1][nP]]
        pObj = CI1TreeItemProject(self.rootItem,infoList=infoList)
        self.rootItem.appendChildProject(pObj)
      pObj.refDbIndex = nP
      pObj.setDirectory(os.path.split(params['PROJECT_DB'][1][nP])[0])
      pObj.loadDummy()
    #print 'I1 projects listed in database:',self._lastLoadedProjectIds
    self.sourceFile = fileName
    return err

  def editDirectoriesDefFile(self,projectDbIndex=None,directory = None):
    err = CErrorReport()
    
    try:
      lineList = CCP4Utils.readFile(self.sourceFile).split('\n')
    except:
      err.append(self.__class__,101,self.sourceFile)
      return err
    
    #print 'editDirectoriesDefFile',projectDbIndex,directory
    if projectDbIndex is not None and directory is not None:
      s1 = 'PROJECT_DB,'+str(projectDbIndex)
      s2 = 'PROJECT_PATH,'+str(projectDbIndex)
      done = False
      for ii in range(len(lineList)):
        key,keyType,value = splitDefLine(lineList[ii])
        #print 'editDirectoriesDefFile',ii,key,'*',keyType,'*',value
        if key is None:
          pass
        elif key == s1:
          lineList[ii] = key + ' ' + keyType + ' ' + os.path.join( directory , 'CCP4_DATABASE' )
          done = True
        elif key == s2:
          lineList[ii] = key + ' ' + keyType + ' ' + directory
          done = True
      if done:
         CCP4Utils.backupFile(self.sourceFile,delete=True)
         text = ''
         for line in lineList: text = text + line +'\n'
         #print 'editDirectories.def',text
         try:
           CCP4Utils.saveFile(self.sourceFile,text=text)
         except Exception as e:
           err.append(self.__class__,103,self.sourceFile)
      else:
        err.append(self.__class__,102,self.sourceFile)
    return err
      

  def addFolder(self,name,parent=None):
    if parent is None: parent = self.rootItem
    fItem = CI1TreeItemFolder(parent=parent, info = { 'name' : name} )
    self.beginInsertRows(QtCore.QModelIndex(),len(self.rootItem.childFolders),len(self.rootItem.childFolders))
    parent.appendChildFolder(fItem)
    self.endInsertRows()
    return fItem
    
  @QtCore.Slot(str)
  def deleteFolder(self,folderId):
    fItem = self.getFolder(folderId)
    if fItem is None:
      print('Failed to find and delete folder',folderId)   
    self.beginResetModel()
    fItem.parent().removeChildFolder(folderId)
    self.endResetModel()
    
  def canFetchMore(self,parent):
    #print 'CI1ProjectModel.canFetchMore',parent,parent.isValid(),parent.internalPointer()
    if not parent.isValid():
      return False
    else:
      return  parent.internalPointer().canFetchMore()

  def fetchMore(self,parent):
    # Load the job/files for any opened project in the tree model
    if not self.canFetchMore(parent): return
    parentNode = parent.internalPointer()
    parentNode.fetchMore()

  def index(self, row, column, parent):
    if row < 0 or column < 0 or row >= self.rowCount(parent) or column >= self.columnCount(parent):
      return QtCore.QModelIndex()
    if not parent.isValid():
      parentItem = self.rootItem
    else:
      parentItem = parent.internalPointer()
    childItem = parentItem.child(row)
    #print 'CProjectModel.index',row, column, parent,childItem.getName()
    if childItem:
      return self.createIndex(row, column, childItem)
    else:
      return QtCore.QModelIndex()
    
  def parent(self, index):
    if not index.isValid():
      return QtCore.QModelIndex()
    childItem = index.internalPointer()
    parentItem = childItem.parent()
    if parentItem == self.rootItem:
      return QtCore.QModelIndex()
    return self.createIndex(parentItem.row(), 0, parentItem)
  
  def rowCount(self, parent):
    #if parent.column() > 0:
    #  return 0
    if not parent.isValid():
      parentItem = self.rootItem
    else:
      parentItem = parent.internalPointer()
    return parentItem.childCount()

  def columnCount(self,parent):
    return 2

  def data(self,index,role):
    if not index.isValid():
      return None    
    item = index.internalPointer()
    return item.data(index.column(),role)

  def headerData(self,section,orientation,role=QtCore.Qt.DisplayRole):
    if orientation == QtCore.Qt.Horizontal:
      if role==QtCore.Qt.DisplayRole:
        if section == 0:
            return 'Project/Job/File'
        elif section == 1:
            return 'Date'
#FIXME PYQT - or maybe None? This used to return QVariant.
    return None

  def flags(self,modelIndex):
    if not modelIndex.isValid(): return QtCore.Qt.NoItemFlags
    ic = modelIndex.column()
    if ic == 0:
      return QtCore.Qt.ItemIsEnabled|QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsDragEnabled|QtCore.Qt.ItemIsDropEnabled
    else:
      return  QtCore.Qt.ItemIsSelectable|QtCore.Qt.ItemIsEnabled


  def nodeFromIndex(self, index):
    if index is not None and index.isValid():
      #print 'nodeFromIndex',index,index.internalPointer()
      return index.internalPointer()
    else:
      return self.rootItem

  def getFolder(self,fId):
    return self.getProject(fId)

  def getProject(self,pId):
    indexList = self.match(self.index(0,0,QtCore.QModelIndex()),QtCore.Qt.DisplayRole,pId,1,
                           QtCore.Qt.MatchExactly|QtCore.Qt.MatchWrap|QtCore.Qt.MatchRecursive)
    #print 'getProjectTreeItem',pId,indexList
    if len(indexList)>0:
      return indexList[0].internalPointer()
    else:
      return None
  
  def modelIndexFromProject(self,pId):
    indexList = self.match(self.index(0,0,QtCore.QModelIndex()),QtCore.Qt.DisplayRole,pId,1,
                           QtCore.Qt.MatchExactly|QtCore.Qt.MatchWrap|QtCore.Qt.MatchRecursive)
    if len(indexList)>0:
      return indexList[0]
    else:
      return None
  
  def getJob(self,pId,jId):
    indexList = self.match(self.index(0,0,QtCore.QModelIndex()),QtCore.Qt.DisplayRole,pId,1,
                           QtCore.Qt.MatchExactly|QtCore.Qt.MatchWrap|QtCore.Qt.MatchRecursive)
    #print 'getJob',pId,jId,indexList
    if len(indexList)==0: return None
    jList = self.match(self.index(0,0,indexList[0]),QtCore.Qt.UserRole,jId,1, QtCore.Qt.MatchExactly|QtCore.Qt.MatchWrap )
    if len(jList)>0:
      return jList[0].internalPointer()
    else:
      return None

  def getJob1(self,jId):
    jList = self.match(self.index(0,0,QtCore.QModelIndex()),QtCore.Qt.UserRole,jId,1, QtCore.Qt.MatchExactly|QtCore.Qt.MatchWrap )
    if len(jList)>0:
      return jList[0].internalPointer()
    else:
      return None

class CI1ProjectView(QtWidgets.QTreeView):

  rightMousePress = QtCore.Signal('QMouseEvent')
  jobClicked = QtCore.Signal('QModelIndex')
  fileClicked = QtCore.Signal(str)

  def __init__(self,parent=None):
    QtWidgets.QTreeView.__init__(self,parent)
    self.setObjectName('i1ProjectWidget')
    self.setSelectionMode(QtWidgets.QAbstractItemView.SingleSelection)
    self.setDragEnabled(True)
    self.setAcceptDrops(True)
    self.setDragDropMode(QtWidgets.QAbstractItemView.InternalMove)
    self.setExpandsOnDoubleClick(False)
    self.setRootIsDecorated(True)
    self.setIconSize(QtCore.QSize(16,16))
    self.setEditTriggers(QtWidgets.QAbstractItemView.EditKeyPressed)
    #self.setFocusPolicy(QtCore.Qt.NoFocus)
    self.setToolTip('Right mouse click for options to view jobs and files')
    self.setAlternatingRowColors(PREFERENCES().TABLES_ALTERNATING_COLOR)
    PREFERENCES().TABLES_ALTERNATING_COLOR.dataChanged.connect(self.resetAlternatingRowColors)
    self.forceUpdate = sys.platform.count('inux') or sys.platform.count('arwin')
    #print 'CProjectView.__init__ forceUpdate',self.forceUpdate

  @QtCore.Slot()
  def update(self):
    if not self.forceUpdate: return
    #print 'CProjectView.update'
    #QtWidgets.QTreeView.update(self)
    # Why am I doing this???  This is broken on Windows
    # Maybe this fixes the failure to update on Linux
    frameRect = self.frameRect()
    self.setDirtyRegion(QtGui.QRegion(frameRect.x(),frameRect.y(),frameRect.width(),frameRect.height()))
    #print 'from CProjectView.update'
    
  def mousePressEvent(self,event=None):
    #print 'mousePressEvent'
    if event.button() == QtCore.Qt.RightButton:
      self.rightMousePress.emit(event)
      event.accept()
      
      return
    else:     
      mousePressX = event.x()
      modelIndex = self.indexAt(event.pos())
      r = self.visualRect(modelIndex)
      if r.isValid():
        indent = 0
        #print 'mousePressEvent mousePressX',mousePressX,'left',r.left(),(r.left() + self.iconSize().width() + 4)
        
        if mousePressX  >r.left()+indent and mousePressX < (r.left() + self.iconSize().width() + 4 + indent):
          self.startDrag(modelIndex=modelIndex)
          event.accept()
          return
        
      node = self.model().nodeFromIndex(modelIndex)
      if node.isJob():
        self.jobClicked.emit(modelIndex)
      else:
        if isinstance(node,CI1TreeItemFile):
          self.fileClicked.emit(node.filePath())
          event.accept()
          return
    #print 'calling QTreeView.mousePressEvent'
    QtWidgets.QTreeView.mousePressEvent(self,event)


  def nodeFromEvent(self,event):
    modelIndex = self.indexAt(QtCore.QPoint(event.x(),event.y()))
    col = self.model().COLUMNS[modelIndex.column()]
    return modelIndex,self.model().nodeFromIndex(modelIndex),col

    return [None,None,None,indx]
  
  def selectRow(self,modelIndex=None):
    #print 'CProjectView.selectRow',modelIndex.row()
    sel = QtCore.QItemSelection( modelIndex.sibling(modelIndex.row(),0) , modelIndex.sibling(modelIndex.row(),2) )
    self.selectionModel().select(sel,QtCore.QItemSelectionModel.ClearAndSelect)
    #print 'CProjectView.selectRow DONE'

  
  def startDrag(self,dropActions=None,modelIndex=None):
    #print 'CI1ProjectView.startDrag'
    if modelIndex is None:
      modelIndex = self.currentIndex()
    if modelIndex is None: return
    node = self.model().nodeFromIndex(modelIndex)
    #print 'startDrag',node
    if isinstance(node,(CI1TreeItemFolder,CI1TreeItemProject,CI1TreeItemFile)):
      mimeData = node.mimeData()
      if mimeData is None: return 
      drag = QtGui.QDrag(self)
      drag.setMimeData(mimeData)    
      icon = node.data(0,QtCore.Qt.DecorationRole)
      if icon is not None:
        try:
          pixmap = icon.pixmap(18,18)
          drag.setHotSpot(QtCore.QPoint(9,9))
          drag.setPixmap(pixmap)
        except:
          pass    
      drag.exec_(QtCore.Qt.CopyAction)
      
  def dragEnterEvent(self,event):
    if event.mimeData().hasFormat('I1Project') or event.mimeData().hasFormat('I1Folder'):
      event.accept()
    else:
      event.ignore()

  def dragMoveEvent(self,event):
    if not self.indexAt(event.pos()).isValid():
      return  event.accept()    
    targetItem = self.indexAt(event.pos()).internalPointer()
    if not isinstance(targetItem,CI1TreeItemFolder):
      event.ignore()
      return
    #if event.mimeData().hasFormat('project') and dropItem is not None:
    if event.mimeData().hasFormat('I1Project') or event.mimeData().hasFormat('I1Folder') :
      event.setDropAction(QtCore.Qt.MoveAction)
      event.accept()
    else:
      event.ignore()

  def dropEvent(self,event):
    #print 'dropEvent',self.indexAt(event.pos())
    if not self.indexAt(event.pos()).isValid():
      # Attempting to drop into space - i.e. move object to root level
      targetItem = self.model().rootItem
    else:
      targetItem = self.indexAt(event.pos()).internalPointer()
    #print 'dropEvent',targetItem
    if not isinstance(targetItem,CI1TreeItemFolder):
      event.setDropAction(QtCore.Qt.IgnoreAction)
      event.ignore()
      return
    if event.mimeData().hasFormat('I1Project') or  event.mimeData().hasFormat('I1Folder'):
      pass
    else:
      event.setDropAction(QtCore.Qt.IgnoreAction)
      event.ignore()
      return
  
    targetFolderId = self.item2FolderId(targetItem)
    #print 'dropEvent targetFolder',targetFolderId
    failed = False
    try:
      self.model().beginResetModel()   
      if event.mimeData().hasFormat('I1Project'):
        movedItemName = self.event2MimeData(event,'I1Project')['name']
        movedItem = self.model().getProject(movedItemName)
        err = movedItem.parent().removeChildProject(movedItemName)
        #print 'dropEvent',movedItem.parent(),movedItem.parent().childProjects
        if len(err)>0: raise err
        #print 'dropEvent movedItem childJobs',movedItem.childJobs
        movedItem.setParent(targetItem)
        targetItem.appendChildProject(movedItem)
        #print 'dropEvent move project',movedItemName,'to',targetItem.name
      else:
        movedItemName = self.event2MimeData(event,'I1Folder')['name']
        movedItem = self.model().getFolder(movedItemName)    
        if targetItem ==  movedItem.parent():
          print('drag-n-drop to same place')
          failed = True       
        err = movedItem.parent().removeChildFolder(movedItemName)
        if len(err)>0:
          print('Failed removing folder from parent',movedItem)
          failed = True
        movedItem.setParent(targetItem)
        targetItem.appendChildFolder(movedItem)
        #print 'dropEvent move folder',movedItemName,'to',targetItem.name
    except CException as e:
      print(e.report())
      failed = True
    except Exception as e:
      print('Failed to move project',e)
      failed = True

    self.model().endResetModel()
    if failed:
      event.setDropAction(QtCore.Qt.IgnoreAction)
      event.ignore()
    else:
      event.setDropAction(QtCore.Qt.MoveAction)
      event.accept()

          
  def item2FolderId(self,item):
    if item is None: return None
    folderId = str(item.data(0,QtCore.Qt.UserRole))
    #print 'item2FolderId',folderId
    return  folderId

  def event2MimeData(self,event,label):
    text = str(event.mimeData().data(label).data())
    #print 'event2MimeData',text
    tree = ET.fromstring(text)
    params = { }
    if tree.tag != 'root' : params[str(tree.tag)]=str(tree.text)
    for child in tree:
      params[child.tag] = child.text
    #print 'event2MimeData',params
    return params
      
  @QtCore.Slot()
  def resetAlternatingRowColors(self):
    self.setAlternatingRowColors(PREFERENCES().TABLES_ALTERNATING_COLOR)
  
class CI1ProjectWidget(QtWidgets.QFrame):

  deleteFolder = QtCore.Signal(str)
  reloadProject = QtCore.Signal(str)
  makeI2Project = QtCore.Signal(str)
  associateI2Project = QtCore.Signal(str)
  collectProject = QtCore.Signal(str)
  findProject = QtCore.Signal(str)
  viewInQtrview = QtCore.Signal(str)
  viewDefFile = QtCore.Signal(str,str)
  viewInCoot = QtCore.Signal(str,str)
  viewInMg = QtCore.Signal(str,str)
  viewFile = QtCore.Signal(str)
  copyFile = QtCore.Signal(str)
  showLogFile = QtCore.Signal(str)

  MARGIN = 0

  def __init__(self,parent=None):
    #print 'CI1ProjectWidget.__init__'
    QtWidgets.QFrame.__init__(self,parent)
    layout = QtWidgets.QVBoxLayout()
    layout.setMargin(CI1ProjectWidget.MARGIN)
    layout.setContentsMargins(CI1ProjectWidget.MARGIN,CI1ProjectWidget.MARGIN,
                                CI1ProjectWidget.MARGIN,CI1ProjectWidget.MARGIN)
    layout.setSizeConstraint(QtWidgets.QLayout.SetMinAndMaxSize)
    self.setLayout(layout)

    #self.tab = QtWidgets.QTabWidget(self)
    #self.layout().addWidget(self.tab)

    model = CI1ProjectModel()
    self.projectView=CI1ProjectView(self)
    self.projectView.setModel(model)
#FIXME - Never emitted in this file and not a child of CProjectModel
    #model.redrawSigal.connect(self.projectView.update)
    self.projectView.model().resetAll("dummy")
    self.projectView.setColumnWidth(0,300)
    
    self.projectView.rightMousePress.connect(self.showPopup)
    #self.tab.addTab(self.projectView,'Job list')
    self.layout().addWidget(self.projectView)

    self.projectView.selectionModel().selectionChanged.connect(self.handleSelectionChanged)
    
    #self.doubleClicked = None
    self.popupMenu = QtWidgets.QMenu(self)
    #print 'DONE CI1ProjectWidget.__init__'

  def model(self):
    return self.projectView.model()

  @QtCore.Slot('QEvent')
  def showPopup(self,event):
    modelIndex,node,column = self.projectView.nodeFromEvent(event)
    #print 'CI1ProjectWidget.showPopup',node.getName(),column
    position = QtCore.QPoint(event.globalX(),event.globalY())
    self.popupMenu.setTitle(node.getName())
    self.popupMenu.clear()
    if node.isFolder():
      action = self.popupMenu.addAction('Delete folder')
      action.triggered.connect(functools.partial(self.deleteFolder.emit,node.name))
      if len(node.childFolders) + len(node.childProjects) > 0: action.setEnabled(False)
    elif node.isProject():
       action = self.popupMenu.addAction('Reload project')
       action.triggered.connect(functools.partial(self.reloadProject.emit,node.refProject))
       action = self.popupMenu.addAction('Make CCP4i2 project')
       action.triggered.connect(functools.partial(self.makeI2Project.emit,node.refProject))
       action = self.popupMenu.addAction('Associate with CCP4i2 project')
       action.triggered.connect(functools.partial(self.associateI2Project.emit,node.refProject))
       action = self.popupMenu.addAction('Collect input files')
       action.triggered.connect(functools.partial(self.collectProject.emit,node.refProject))
       if os.path.exists(os.path.join(node.directory,'IMPORTED_FILES')): action.setEnabled(False)
       action = self.popupMenu.addAction('Find project directory')
       action.triggered.connect(functools.partial(self.findProject.emit,node.refProject))
       if node.broken>=0: action.setEnabled(False)
    elif node.isJob():
      popupViewSubMenu = self.popupMenu.addMenu('View')
      action = popupViewSubMenu.addAction('Job results (new style,qtrview)')
      action.triggered.connect(functools.partial(self.viewInQtrview.emit,node.logFile()))
      action = popupViewSubMenu.addAction('Def file')
      action.triggered.connect(functools.partial(self.viewDefFile.emit,node.parent().refProject,node.jobId))
      action = popupViewSubMenu.addAction('In Coot')
      if not node.finished: action.setEnabled(False)
      action.triggered.connect(functools.partial(self.viewInCoot.emit,node.parent().refProject,node.jobId))
      action = popupViewSubMenu.addAction('In CCP4mg')
      if not node.finished: action.setEnabled(False)
      action.triggered.connect(functools.partial(self.viewInMg.emit,node.parent().refProject,node.jobId))
    elif node.isFile():
      action = self.popupMenu.addAction('View')
      action.triggered.connect(functools.partial(self.viewFile.emit,node.filePath()))
      action = self.popupMenu.addAction('Copy')
      action.triggered.connect(functools.partial(self.copyFile.emit,node.filePath()))
    self.popupMenu.popup(position)
    
  @QtCore.Slot('QItemSelection','QItemSelection')
  def handleSelectionChanged(self,selected,deselected):
    indices = selected.indexes()
    if len(indices)>0:
      try:
        node = self.projectView.model().nodeFromIndex(indices[0])
        #jobId = node.jobId
        #projectId = node.parent().getProjectId()
        logFile = node.logFile()
      except:
        return
      else:
        self.showLogFile.emit(logFile)
      
class CI1ProjectViewer(CCP4WebBrowser.CMainWindow):

  Instances = []
  MARGIN = 2

  ERROR_CODES = { 101 : { 'description' : 'Failed to find/open old CCP4 master directories.def file' },
                  102 : { 'description' : 'Failed to find/open old CCP4 def file' },
                  103 : { 'description' : 'Failed loading task titles from modules definition file' },
                  104 : { 'description' : 'Unknown failure loading project database file' },
                  105 : { 'severity' : SEVERITY_WARNING, 'description' : 'Project with this name already loaded - will overwrite' },
                  106 : { 'severity' : SEVERITY_WARNING, 'description' : 'Project directory does not exist' },
                  107 : { 'description' : 'Folder with this name already exists' },
                  110 : { 'severity' : SEVERITY_WARNING, 'description' : 'Failed loading subsiduary task def file' },
                  111 : { 'severity' : SEVERITY_WARNING, 'description' : 'Failed parsing line of def file' },
                  112 : { 'severity' : SEVERITY_WARNING, 'description' : 'Failed loading CCP4 task def file to extract type info' },
                  121 : { 'description' : 'Input file not found' },
                  122 : { 'severity' : SEVERITY_OK, 'description' : 'Copied' }
                  }
  def __init__(self,parent=None,fileName=None):
    CCP4WebBrowser.CMainWindow.__init__(self,parent)
    CI1ProjectViewer.Instances.append(self)
    self.setObjectName('projectViewer')
    self.preferences = CI1Preferences()
    CCP4ProjectWidget.CTreeItemJob.TODAY = time.strftime(CCP4ProjectWidget.CTreeItemJob.DATE_FORMAT,time.localtime())
    CCP4ProjectWidget.CTreeItemJob.THISYEAR = time.strftime('%y',time.localtime())
    self.actionDefinitions = {
                  
                  'i1_help' :
                               { 'text' : 'Help',
                                 'tip' : 'Help on viewing old CCP4i projects',
                                 'slot' : self.showHelp },
                  'show_load_errors' :
                               { 'text' : 'Show load errors',
                                 'tip' : 'List problems in loading CCP4i projects',
                                 'slot' : self.showLoadErrors },
                  'i1preferences' :
                               { 'text' : 'Preferences',
                                 'tip' : 'CCP4i project viewer preferences',
                                 'enabled' : True,
                                 'slot' : self.showPreferences },
                   'add_folder' :
                            { 'text' : 'Add folder',
                              'tip' : 'Create a folder to organise projects',
                              'slot' : self.createFolder }
                            }
    self.destroyed.connect(CI1ProjectViewer.updateInstances)

    self.setWindowTitle(self.version+'Old CCP4i Project Viewer')
    self.layout().setContentsMargins(CI1ProjectViewer.MARGIN,CI1ProjectViewer.MARGIN,CI1ProjectViewer.MARGIN,CI1ProjectViewer.MARGIN)
    self.layout().setSpacing(CI1ProjectViewer.MARGIN)

    # left side project widget and buttons
    leftFrame = QtWidgets.QFrame(self)
    leftFrame.setLayout( QtWidgets.QVBoxLayout())
    leftFrame.layout().setContentsMargins(CI1ProjectViewer.MARGIN,CI1ProjectViewer.MARGIN,CI1ProjectViewer.MARGIN,CI1ProjectViewer.MARGIN)
    leftFrame.layout().setSpacing(CI1ProjectViewer.MARGIN)

    self._projectWidget = CI1ProjectWidget(self)
    self._projectWidget.setMinimumSize(QtCore.QSize(370,300))
    leftFrame.layout().addWidget(self._projectWidget)
    self._projectWidget.showLogFile.connect(self.showLogFile)
    self._projectWidget.projectView.fileClicked.connect(self.handleFileClicked)

    self._projectWidget.viewInQtrview.connect(self.viewInQtrview)
    self._projectWidget.viewInCoot.connect(functools.partial(self.viewInGraphics,'coot'))
    self._projectWidget.viewInMg.connect(functools.partial(self.viewInGraphics,'ccp4mg'))
    self._projectWidget.viewDefFile.connect(self.viewDefFile)
    self._projectWidget.collectProject.connect(self.collectProject)
    #self._projectWidget.annotateProject.connect(self.annotateProject)
    self._projectWidget.makeI2Project.connect(self.makeI2Project)
    self._projectWidget.associateI2Project.connect(self.handleAssociateI2Project)
    self._projectWidget.reloadProject.connect(self.reloadProject)
    self._projectWidget.findProject.connect(self.findProjectDir)
    self._projectWidget.viewFile.connect(self.handleFileClicked)
    self._projectWidget.copyFile.connect(self.copyFile)
    self._projectWidget.deleteFolder.connect(self.model().deleteFolder)

    self.rightFrame = QtWidgets.QFrame(self)
    self.rightFrame.setLayout(QtWidgets.QVBoxLayout())
    self.qtrButton = QtWidgets.QPushButton('Show log file including graphs in qtrview',self)
    self.qtrButton.clicked.connect(functools.partial(self.viewInQtrview,'CURRENTLOG'))
    self.rightStack = QtWidgets.QStackedWidget(self)  
    self.webView= CCP4WebView.CWebView(self,blockLoading=True)
    self.rightStack.addWidget(self.webView)
    self.textView = CCP4TextViewer.CTextViewer(self)
    self.rightStack.addWidget(self.textView)
    self.rightFrame.layout().addWidget(self.rightStack)
    self.rightFrame.layout().addWidget(self.qtrButton)

    # left/right splitter
    from qtgui import CCP4TaskWidget
    self.splitterSizes = [400,CCP4TaskWidget.WIDTH+CCP4TaskWidget.PADDING_ALLOWANCE]
    mainWidget = QtWidgets.QSplitter(self)
    mainWidget.setOrientation(QtCore.Qt.Horizontal)
    mainWidget.addWidget(leftFrame)
    mainWidget.addWidget(self.rightFrame)
    # And set the splitter sizes after show ...

    self.mainToolBar = QtWidgets.QToolBar(self)
    self.setUnifiedTitleAndToolBarOnMac(False)
    self.show_load_errors_action = self.mainToolBar.addAction("Show load errors")
    self.show_load_errors_action.setToolTip('List problems in loading CCP4i projects')
    self.show_load_errors_action.triggered.connect(self.showLoadErrors)
    self.i1_help_action = self.mainToolBar.addAction("Help")
    self.i1_help_action.setToolTip('Help on viewing old CCP4i projects')
    self.i1_help_action.triggered.connect(self.showHelp)
    self.addToolBar(self.mainToolBar)

    self.updateLoadErrorsEnabled()
    
    self.setCentralWidget(mainWidget)
    self.setCurrentLogFile()
    #print 'to loadDirectoriesDefFile',fileName
    # Try reading a supplement file
    
    supFile = os.path.join(os.path.split(fileName)[0],'i2supplement.xml')
    if os.path.exists(supFile): self.model().loadSupplement(supFile)

    self.model().beginResetModel()
    self.model().loadDirectoriesDefFile(fileName=fileName)
    self.model().endResetModel()

  def model(self):
    return self._projectWidget.projectView.model()

  def handleProjectMenuExport(self):
      pass

  def updateLoadErrorsEnabled(self):
    enabled = (len(self._projectWidget.projectView.model().loadErrorReport)>0)
    #print 'updateLoadErrorsEnabled',enabled
    self.show_load_errors_action.setEnabled(enabled)

  def showLoadErrors(self):
    message = 'Errors loading old CCP4i projects:\n'
    for proj,dir in self._projectWidget.projectView.model().loadErrorProjects:
      message += '{:30} {:100}\n'.format(proj,dir)
    self._projectWidget.projectView.model().loadErrorReport.warningMessage(parent=self,windowTitle=self.windowTitle(),message=message)

  def showHelp(self):
    WEBBROWSER().loadWebPage(helpFileName='CCP4i1Projects.html',newTab=True)

  def handleCurrentJobChanged(self,jobId,projectId):
    #print 'handleCurrentJobChanged',jobId,projectId
    pass

  def setCurrentLogFile(self,logFile=None):
    if logFile is None:
      colour = 'grey'
      self.qtrButton.setEnabled(False)
    else:
      colour = 'black'
      self.qtrButton.setEnabled(True)
    self.currentLogFile = logFile

  @QtCore.Slot(str)
  def showLogFile(self,logFile):
    #print 'showLogFile',logFile
    if logFile is None:
      return
    elif logFile.endswith('html') or  logFile.endswith('htm'):
      self.webView.load(QtCore.QUrl.fromLocalFile(logFile))
      self.rightStack.setCurrentIndex(0)
    else:
      self.textView.loadText(CCP4Utils.readFile(logFile))
      self.rightStack.setCurrentIndex(1)
    self.setCurrentLogFile(logFile)

  @QtCore.Slot(str,str)
  def handleFileClicked(self,fileName=None,fileType=None):
    if fileName is None: return
    if fileType is None:
      ext = os.path.splitext(fileName)[1]
      if ext == '.mtz':
        fileType = "application/CCP4-mtz"
      elif ext == '.pdb':
        fileType = "chemical/x-pdb"
      else:
        fileType = "text/plain"
    if fileType == "application/CCP4-mtz":
      LAUNCHER().launch('viewhkl',[fileName])
    elif fileType == "chemical/x-pdb":
      self.textView.loadText(CCP4Utils.readFile(fileName))
      self.rightStack.setCurrentIndex(1)
      self.setCurrentLogFile(None)
    else:
      self.textView.loadText(CCP4Utils.readFile(fileName))
      self.rightStack.setCurrentIndex(1)
      self.setCurrentLogFile(None)
              
  @QtCore.Slot(str)
  def viewInQtrview(self,logFile=None):
    if logFile is None:
      return
    elif logFile == 'CURRENTLOG':
      if self.currentLogFile is None: return
      import copy
      logFile = copy.deepcopy(self.currentLogFile)
    if os.path.splitext(logFile)[1] == '.html':
        logFile = os.path.splitext(logFile)[0]
    LAUNCHER().launch('logview',[logFile])

  @QtCore.Slot(str,str,str)
  def viewInGraphics(self,mode='coot',projectId=None,jobId=None):
    #print 'viewInGraphics',mode,projectId,jobId
    jItem = self.model().getJob(projectId,jobId)
    if jItem is None:
      print('ERROR in viewInGraphics - could not find job:',projectId,jobId)
      return
    fileList = []
    for fItem in jItem.childOutFiles:
      if fItem.broken>0:
        path = fItem.filePath()
        if os.path.splitext(path)[1] in ['.mtz','.pdb']:
          fileList.append(fItem.filePath())
    if mode == 'coot' : mode = 'coot0'
    LAUNCHER().openInViewer(viewer=mode,fileName=fileList)

  @QtCore.Slot(str,str)
  def viewDefFile(self,projectId=None,jobId=None):
    #print 'viewDefFile',projectId,jobId
    p = self.model().getProject(projectId)
    if p is None:
      print('ERROR in viewDefFile - could not find project:',projectId)
      return

    j = self.model().getJob(projectId,jobId)
    if j is None:
      print('ERROR in viewDefFile - could not find job:',projectId,jobId)
      return
    
    defFile =os.path.join( p.directory,'CCP4_DATABASE',str(jobId)+'_'+j.taskName+'.def')
    #print 'viewDefFile',defFile
    self.textView.loadText(CCP4Utils.readFile(defFile),fixedWidthFont=True)
    self.rightStack.setCurrentIndex(1)

  @QtCore.Slot(str)
  def findProjectDir(self,projectId):
    from qtgui import CCP4FileBrowser
    self.findProjectDialog = CCP4FileBrowser.CFileDialog(self,fileMode=QtWidgets.QFileDialog.Directory,
      title='Find project directory for '+projectId)
    self.findProjectDialog.selectFile.connect(functools.partial(self.updateProjectDir,projectId))
    self.findProjectDialog.show()

  @QtCore.Slot(str,str)
  def updateProjectDir(self,projectId,directory):
    if os.path.exists(directory) and os.path.exists(os.path.join(directory,'CCP4_DATABASE')):
      pass
    elif os.path.split(directory)[1] == 'CCP4_DATABASE':
      directory = os.path.split(directory)[0]
    else:
      mess = QtWidgets.QMessageBox.warning(self,'Find project directory for '+projectId,'The selected directory does not contain a CCP4_DATABASE sub-directory - it is not a valid CCP4 project directory.')
      return

    self.model().getProject(projectId).setDirectory(directory)
    err = self.model().getProject(projectId).editDatabaseDef(resetDir=True)
    err.append(self.model().editDirectoriesDefFile(projectDbIndex=self.model().getProject(projectId).refDbIndex,directory=directory))
    if err.maxSeverity()>SEVERITY_WARNING:
      mess = QtWidgets.QMessageBox.warning(self,'Find project directory for '+projectId,"Failed saving changed project directory to old CCP4 files")
    return


  @QtCore.Slot(str)
  def collectProject(self,projectId):
    rv = QtWidgets.QMessageBox.question(self,'Collect input files for '+str(projectId),"Copy input files to project 'IMPORTED_FILES' directory.\nPlease close old CCP4i before running this.",QtWidgets.QMessageBox.Yes|QtWidgets.QMessageBox.No)
   
    if rv != QtWidgets.QMessageBox.Yes: return
    pObj = self.model().getProject(projectId)

    err = CErrorReport()
    import shutil
    
    importDir = os.path.join(pObj.directory,'IMPORTED_FILES')
    if not os.path.exists(importDir):
      try:
        os.mkdir(importDir)
      except:
        QtWidgets.QMessageBox.warning(self,'Failed creating IMPORTED_FILES directory','Do you have write permission for '+str(pObj.directory))
        return
    elif not os.access(importDir , os.W_OK|os.X_OK):
      QtWidgets.QMessageBox.warning(self,'Can not write to IMPORTED_FILES directory','You do not have write permission for '+str(pObj.directory))
      return
    
    dbFile = os.path.normpath(os.path.join(pObj.directory,'CCP4_DATABASE','database.def'))
    metaData,params,err = readI1DefFile(dbFile)
    if err.maxSeverity()>SEVERITY_WARNING or params.get('NJOBS',None) is None:
      return

    sourceLookup = {}
    movedFiles = {}
    
    nJobs =int( params['NJOBS'][1] )
    for nJ in range(1,nJobs+1):
      fList = params['INPUT_FILES'][1][nJ]
      dList = params['INPUT_FILES_DIR'][1][nJ]
      #print 'collectProject',nJ,dList,fList
      if fList is not None and dList is not None:
        fList = fList.split()
        dList = dList.split()
        for nF in range(min(len(fList),len(dList))):
          f = fList[nF].strip('"')
          d = dList[nF].strip('"')
          #print 'collectProject',nF,d,f
          if d == 'FULL_PATH':
            #print 'collectProject to copy',nJ,nF,f
            if os.path.exists(f):
              if os.path.split(f)[1] in sourceLookup:
                #Already inported file of same name
                if sourceLookup[os.path.split(f)[1]] == f:
                  target = os.path.join(importDir,os.path.split(f)[1])
                else:
                  # Same file name but from different source??
                  base,ext0 = os.path.split(f)[1].splitext()
                  base,ext1 = os.path.split(base)[1].splitext()
                  idx = 1
                  while base+'_' + ('00'+str(idx))[-2:] + ext1 + ext0 in sourceLookup:
                    idx += 1
                  filn=  base+'_' + ('00'+str(idx))[-2:] + ext1 + ext0 
                  #print 'collectProject incr target file',idx,filn
                  target = os.path.join(importDir,filn)
                  shutil.copyfile(f,target)
                  sourceLookup[filn] = f
              else:
                target = os.path.join(importDir,os.path.split(f)[1])
                #print 'copying to',target
                shutil.copyfile(f,target)
                sourceLookup[os.path.split(f)[1]] = f
              fList[nF] = target
              movedFiles[str(nJ)] = fList
              err.append(self.__class__,122,f+' input to job '+str(nJ))
            else:
              err.append(self.__class__,121,f+' input to job '+str(nJ))

    #print 'collectProject movedFiles'
    #for key,value in movedFiles.items(): print key,value

    if len(movedFiles)>0:
      err0 = pObj.editDatabaseDef(movedFiles=movedFiles)
      if err0.maxSeverity()>SEVERITY_WARNING:
        err0.warningMessage('Collect input files','There was an error saving the new filenames to the old CCP4 database',parent=self)

    self.reloadProject(projectId)
    
  @QtCore.Slot(str)
  def reloadProject(self,projectId):
    pObj = self.model().getProject(projectId)
    self.model().beginResetModel()
    pObj.loadDatabase()
    self.model().endResetModel()

  @QtCore.Slot(str)
  def annotateProject(self,projectId):
    #print 'annotateProject',projectId
    from core import CCP4Annotation,CCP4DataManager
    pObj = self.model().getProject(projectId)
    if pObj is None: return
    if getattr(self,'annotateWindow',None) is None:
      self.annotateWindow = QtWidgets.QDialog(self)
      self.annotateWindow.setLayout(QtWidgets.QVBoxLayout())
      self.tagList = CCP4Annotation.CMetaDataTagList(parent=self,
                       subItem_enumeratorsFunction = CI1PREFERENCES().getTagList,
                       subItem_addEnumeratorFunction = CI1PREFERENCES().addTag  )
      # Beware CList.set() rebuilds list so must do this before creating the widget
      #print 'annotateProject tagList',self.tagList.__dict__['_value']
      self.annotation = CCP4Annotation.CAnnotation(parent=self)
      self.tagListWidget = CCP4DataManager.DATAMANAGER().widget(model=self.tagList,parentWidget=self.annotateWindow,
                                   qualifiers= { 'listVisible':True,
                                              'title' : 'Choose keyword tags for project - new tags can be entered' })
      self.annotateWindow.layout().addWidget(self.tagListWidget)
      self.annotationWidget = CCP4DataManager.DATAMANAGER().widget(model=self.annotation,parentWidget=self.annotateWindow,
                                                    qualifiers= { 'multiLine' : True,
                                                                  'title' : 'Enter description of project' })
      self.annotateWindow.layout().addWidget(self.annotationWidget)
      bb = QtWidgets.QDialogButtonBox(self.annotateWindow)
      but = bb.addButton(QtWidgets.QDialogButtonBox.Save)
      but.clicked.connect(self.saveAnnotation)
      but = bb.addButton(QtWidgets.QDialogButtonBox.Cancel)
      but.clicked.connect(self.annotateWindow.close)
      self.annotateWindow.layout().addWidget(bb)
    
    self.annotateWindow.setWindowTitle('Annotate CCP4i project: '+projectId)
    self.annotationProject = projectId

    if pObj.annotation is not None:
      self.annotation.set(str(pObj.annotation))
    else:
      self.annotation.unSet()
    self.tagList.unSet()
    self.tagList.set(pObj.tagList)
    #print 'pObj.tagList',repr(pObj),pObj.refProject,pObj.tagList
    #print 'self.tagList set to',self.tagList
    self.annotationWidget.updateViewFromModel()
    self.tagListWidget.updateViewFromModel()
    self.tagListWidget.handleRowChange(row=0,force=True)
    
    self.annotateWindow.show()
    self.annotateWindow.raise_()

  @QtCore.Slot(str)
  def makeI2Project(self,projectId):
    try:
      pid = CCP4Modules.PROJECTSMANAGER().db().getProjectId(projectName=projectId)
    except:
      reqNewName = False
    else:
      reqNewName = True

    from qtgui import CCP4ProjectManagerGui
    if CCP4ProjectManagerGui.CNewProjectGui.insts is None:
      CCP4ProjectManagerGui.CNewProjectGui.insts = CCP4ProjectManagerGui.CNewProjectGui()
      CCP4ProjectManagerGui.CNewProjectGui.insts.projectCreated.connect(functools.partial(self.associateI2Project,projectId))
    else:
      CCP4ProjectManagerGui.CNewProjectGui.insts.clear()
    CCP4ProjectManagerGui.CNewProjectGui.insts.name.set(projectId)
    CCP4ProjectManagerGui.CNewProjectGui.insts.show()

  @QtCore.Slot(str)
  def handleAssociateI2Project(self,i1ProjectName):
    print('handleAssociateI2Project',i1ProjectName)
    from qtgui import CCP4ProjectManagerGui
    self.associateDialog = QtWidgets.QDialog(self)
    self.associateDialog.setLayout(QtWidgets.QVBoxLayout())
    self.associateDialog.layout().addWidget(QtWidgets.QLabel('Choose CCP4i2 project to associate with project '+i1ProjectName,self.associateDialog))

    #FIXME - Too complicated
    self.i2Selector =  CCP4ProjectManagerGui.CProjectsTreeView(self)
    projectsViewProxyModel = CCP4ProjectManagerGui.CProjectsViewProxyModel()
    CCP4ProjectManagerGui.repopulateTreeViewNew(projectsViewProxyModel,self.i2Selector)
    pvd = CCP4ProjectManagerGui.ProjectsViewDelegate()
    self.i2Selector.setItemDelegate(pvd)
    self.i2Selector.setDragEnabled(False)
    self.i2Selector.setAcceptDrops(False)
    #self.i2Selector.populate()

    self.associateDialog.layout().addWidget(self.i2Selector)
    butBox = QtWidgets.QDialogButtonBox(self.associateDialog)
    b = butBox.addButton(QtWidgets.QDialogButtonBox.Ok)
    b.clicked.connect(functools.partial(self.doAssociateDialog,i1ProjectName))
    b = butBox.addButton(QtWidgets.QDialogButtonBox.Cancel)
    b.clicked.connect(self.closeAssociateDialog)
    self.associateDialog.layout().addWidget(butBox)
    self.associateDialog.show()
    self.associateDialog.raise_()

    
  @QtCore.Slot()
  def closeAssociateDialog(self):
    self.associateDialog.close()
    self.associateDialog.deleteLater()
    del self.associateDialog

  @QtCore.Slot(str)
  def doAssociateDialog(self,i1ProjectName):
    projectId = self.i2Selector.selectedProjectId()
    if projectId is None: return
    self.associateI2Project(i1ProjectName,projectId)
    self.closeAssociateDialog()
  
  @QtCore.Slot(str,str)
  def associateI2Project(self,i1ProjectName,projectId=None):
    #print 'associateI2Project',i1ProjectName,projectId
    pObj = self.model().getProject(i1ProjectName)
    PROJECTSMANAGER().db().updateProject(projectId,key='I1ProjectName',value=i1ProjectName)
    PROJECTSMANAGER().db().updateProject(projectId,key='I1ProjectDirectory',value=pObj.directory)
    

  @QtCore.Slot()
  def saveAnnotation(self):
    #print 'saveAnnotation',self.tagList,self.annotation
    self.annotateWindow.close()
    pObj = self.model().getProject(self.annotationProject)  
    pObj.tagList = []
    for item in self.tagList:
      if item.tag.isSet(): pObj.tagList.append(str(item.tag))
    pObj.annotation = str(self.annotation.text)
    # Remove any tags that were created but not ultimately used for this pObj
    CI1PREFERENCES().pruneTags( pObj.tagList )
    #print 'saveAnnotation object',repr(pObj),pObj.refProject,pObj.tagList,pObj.annotation

  def getActionDef(self,name,**info):
    return self.actionDefinitions.get(name,dict(text=name))

  def loadDb(self,overwriteProject=None):
    from qtgui import CCP4FileBrowser
    dialog = CCP4FileBrowser.CFileDialog(self,
           title='Load old CCP4 project(s) from directories.def or database.def',
           filters= [ 'Old CCP4 directory.def or database.def (*.def)' ],
           defaultSuffix='.def',
           fileMode=QtWidgets.QFileDialog.ExistingFile  )
#Not possible with native browser as far as I know. SJM 22/11/2018.
    if not PREFERENCES().NATIVEFILEBROWSER:
        dialog.widget.fileDialog.setFilter(QtCore.QDir.AllEntries | QtCore.QDir.Hidden | QtCore.QDir.NoDotAndDotDot )
    dialog.show()

  def showPreferences(self):
    pass

  def createFolder(self,text=''):
    if getattr(self,'folderWidget',None) is None:
      d = QtWidgets.QDialog(self)
      d.setWindowTitle('I1 Projects Folder')
      d.setModal(True)
      d.setLayout(QtWidgets.QVBoxLayout())
      line = QtWidgets.QHBoxLayout()
      line.addWidget(QtWidgets.QLabel('Enter unique name for folder',self))
      self.folderWidget = QtWidgets.QLineEdit(self)
      line.addWidget(self.folderWidget)
      d.layout().addLayout(line)
      bb = QtWidgets.QDialogButtonBox(d)
      b = bb.addButton(QtWidgets.QDialogButtonBox.Ok)
      b.setFocusPolicy(QtCore.Qt.NoFocus)
      b.clicked.connect(self.handleCreateFolder)
      b = bb.addButton(QtWidgets.QDialogButtonBox.Cancel)
      b.setFocusPolicy(QtCore.Qt.NoFocus)
      b.clicked.connect(d.hide)
      bb.setCenterButtons(True)
      d.layout().addWidget(bb)
    self.folderWidget.setText(text)
    self.folderWidget.window().show()

  @QtCore.Slot()
  def handleCreateFolder(self):
    #print 'handleCreateFolder',self.folderWidget.text()
    name = str(self.folderWidget.text())
    indx = self.model().getProject(name)
    if indx is not None:
      QtWidgets.QMessageBox.warning(self,'Create folder','The folder name is not unique')
      return
    self.folderWidget.window().hide()
    try:
      self.model().addFolder(name)
    except CException as e:
      e.warningMessage(parent=self,windowTitle=self.windowTitle(),message='Failed to create new folder')

  @QtCore.Slot(str)
  def copyFile(self,fileName):
    #print 'CI1ProjectViewer.copyFile',fileName
    if os.path.exists(fileName):
      urlList = [QtCore.QUrl()]
      urlList[0].setPath(fileName )
      mimeData = QtCore.QMimeData()
      mimeData.setUrls(urlList)
      QTAPPLICATION().clipboard().setMimeData(mimeData)
    
  def close(self):
    self.Exit()
    CCP4WebBrowser.CMainWindow.close(self)
    
  def Exit(self):
    self.model().saveStatus()
    

  def handleSave(self):
    pass
        
  def isProjectViewer(self):
    return False

  def widgetIsSaveable(self):
    return False
  def widgetIsPrintable(self):
    return False
  def widgetIsRunable(self):
    return False
  def widgetIsSearchable(self):
    return False
  def handlePrint(self):
    pass
  def handleRun(self):
    pass
  def openFind(self):
    pass
  def isFindFrameOpen(self):
    pass
  def deleteTab(self):
    pass
  def historyBack(self):
    pass
  def historyForward(self):
    pass
  def reloadPage(self):
    pass
  def openManageImportFiles(self):
    pass
  def openApplication(self):
    pass
  def openSendReport(self):
    pass
  def openUpdate(self):
    pass
  @staticmethod
  @QtCore.Slot()
  def updateInstances(qobj):
    l = []
    for w in CI1ProjectViewer.Instances:
      if isAlive(w): l.append(w)   
    CI1ProjectViewer.Instances = l

# Put preferences in a class - could be subclassed from CContainer if
# requirements expand to warrant that
class CI1Preferences:
  insts = None
  def __init__(self):
    CI1Preferences.insts = self
    self.tagList =  [] 
    self.sources = []
    self.load()
    self.appendedTags = []

  def save(self,fileName=None):
    if fileName is None:
      fileName = os.path.join(CCP4Utils.getDotDirectory(),'i1supplement','preferences.xml')
    #print 'CI1Preferences.save',fileName
    f = CCP4File.CI2XmlDataFile(fullPath=fileName)
    if not os.path.exists(fileName):
      f.header.setCurrent()
      f.header.function.set('UNKNOWN')
      f.header.comment = 'Preferences for I1 Project Viewer'
    body = ET.Element('body')
    fEleList = ET.Element('tagList')
    body.append(fEleList)
    #print 'CI1Preferences.save tagList',self.tagList
    for tag in self.tagList:
      if tag is not None:
        fEleList.append(ET.Element('tag'))
        fEleList[-1].text = tag
    f.saveFile(bodyEtree=body)

  def load(self,fileName=None):
    if fileName is None:
      fileName = os.path.join(CCP4Utils.getDotDirectory(),'i1supplement','preferences.xml')    
    if os.path.exists(fileName):
      f = CCP4File.CI2XmlDataFile(fullPath=fileName)
      root = f.getBodyEtree()
      if root.find('tagList') is not None:
        for fEle in root.find('tagList').findall('tag'):
          self.tagList.append(str(fEle.text))
          
  def getTagList(self):
    return self.tagList

  def addTag(self,tag):
    if tag is None or len(tag)==0: return
    if self.tagList.count(tag): return self.tagList.index(tag)
    self.tagList.append(tag)
    self.appendedTags.append(tag)
    self.tagList.sort()
    return self.tagList.index(tag)

  def pruneTags(self,usedTags):
    # Remove any tags not eventually used in the annotation edit window
    for item in self.appendedTags:
      if item not in usedTags:
        if item in self.tagList: self.tagList.remove(item)
    self.appendedTags = []
        
