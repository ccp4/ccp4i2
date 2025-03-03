from __future__ import print_function

import os,shutil
from dbapi import CCP4DbApi,CCP4DbUtils
from core import CCP4File
from core.CCP4ErrorHandling import *
from core.CCP4Modules import PROJECTSMANAGER,TASKMANAGER,QTAPPLICATION,JOBCONTROLLER
from core.CCP4QtObject import CObject
from PySide2 import QtCore

# Kick of from pyi2 prompt with something like..
# import CCP4TestDb; t = CCP4TestDb.CTestDb(); t.createProject('simple1')
TEST_MASTER = '/Users/stuart/Desktop/TEST_MASTER'

class CTestDb(CObject):

  createProjectFinishedSignal = QCore.Signal(str)
  exportProjectFinishedSignal = QCore.Signal(str,str)

  def __init__(self,masterDir=None,parent=None):
    if parent is None: parent =  QTAPPLICATION(graphical=False)
    CObject.__init__(self,parent)
    print(PROJECTSMANAGER())
    self.dbApi = PROJECTSMANAGER().db()
    self.dbApi.jobFinished.connect(self.handleJobFinished)
    self.masterDir = None
    self.currentDb = None
    self.projectIds = {}
    self.createProjectName = None

    from utils import startup
    startup.startJobController()
    PROJECTSMANAGER().startCheckForFinishedJobs()    
    
    if masterDir is None: masterDir = TEST_MASTER
      
    db1 = os.path.join(masterDir,'db1.sqlite')
    if os.path.exists(masterDir) and os.path.exists(db1):
      self.masterDir = masterDir
    else:
      self.setup(masterDir)

    self.setDatabase('1')
          
  def setup(self,masterDir,nDabases=2):
    #Set up two 'databases' and two master project directories
    if not os.path.exists(masterDir):
        os.makedirs(masterDir)
        self.masterDir = masterDir

    from core import CCP4Config
    for mode in ['1','2','3','4','5'][0:nDabases]:
      os.mkdir(self.projectsDir(mode))
      self.dbApi.openDb(self.sqliteFile(mode),createDb=True)      
      config = CCP4Config.CConfig(dbFile=self.sqliteFile(mode))
      config.saveDataToXml(self.configFile(mode))

  def sqliteFile(self,mode):
    return os.path.join(self.masterDir,'db'+mode+'.sqlite')

  def projectsDir(self,mode,projectName=None):
    return os.path.join(self.masterDir,'DB'+mode)

  def configFile(self,mode):
    return os.path.join(self.masterDir,'db'+mode+'_config.xml')

  def setDatabase(self,mode='1'):
    if self.currentDb == mode: return
    self.dbApi.openDb(self.sqliteFile(mode),createDb=True)
    JOBCONTROLLER().setConfigFile(self.configFile(mode))
    self.currentDb = mode
    print('Setting database to',mode)
    
  def removeDatabase(self,mode):
     if self.currentDb == mode:
       self.dbApi.openDb(fileName=None) # Close db
       self.currentDb = None

     try:
       os.remove(self.sqliteFile(mode))
     except:
       pass
     try:
       shutil.rmtree(self.projectsDir(mode))
     except:
       pass

  def createProject(self,name=None):
    # Attempt to delete any existing project of same name - should fail quietly if there is no existing project
    err = PROJECTSMANAGER().deleteProject(projectName=name,deleteDirectory=True)
    self.projectIds[name] = PROJECTSMANAGER().createProject(name,os.path.join(self.masterDir,'DB'+self.currentDb,name))
    if name is not None: self.runJobs(name)
 

  def runJobs(self,projectName=None,scriptName=None):
    import glob,CCP4Utils
    self.jobObjList = []
    self.paramsFileList = glob.glob(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','test_projects',scriptName,'*_input_params.xml'))
    if len(self.paramsFileList)==0:
      print('No params files found in',os.path.join(CCP4Utils.getCCP4I2Dir(),'test','test_projects',scriptName))
      return
    self.paramsFileList.sort()
    self.createProjectName = projectName
    if self.createProjectName not in self.projectIds:
      self.projectIds[self.createProjectName] = self.dbApi.getProjectId(self.createProjectName)

    self.startJob(self.paramsFileList[0])


  def startJob(self,paramsFile):
    fileObj  = CCP4File.CI2XmlDataFile(fullPath=paramsFile)
    fileObj.loadFile()
    taskName = fileObj.header.pluginName.__str__()
    self.jobObjList.append(CCP4DbUtils.COpenJob(projectId = self.projectIds[ self.createProjectName]) )    
    self.jobObjList[-1].createJob(taskName=taskName)
    print('CTestDb.startJob',taskName,self.jobObjList[-1].jobId)
    paramsFile = self.editParamsFile(fileObj,self.jobObjList[-1].jobId)
    self.jobObjList[-1].loadParams(paramsFile)
    self.jobObjList[-1].runJob()

  def editParamsFile(self,fileObj,jobId):
    #If the params file has a project set to $THISPROJECTID then substitute in
    # the current projectId and save to a temporary file
    root = fileObj.getEtreeRoot()
    body = root.xpath('//ccp4i2_body')[0]
    projectNodeList = body.xpath('//project')
    ifChanged = False
    for projectNode in projectNodeList:
      if projectNode.text == '$THISPROJECTID':
        projectNode.text = self.projectIds[self.createProjectName]
        ifChanged = True
    if ifChanged:
      fileObj.setFullPath(os.path.join(PROJECTSMANAGER().jobDirectory(jobId=jobId),'tmp_input_params.xml'))
      fileObj.saveFile(bodyEtree=body)
    return fileObj.__str__()
   

  @QtCore.Slot(dict)
  def handleJobFinished(self,args):
    #print 'CTestDb.handleJobFinished',args
    if args['jobId'] == self.jobObjList[-1].jobId:
      if len(self.jobObjList)<len(self.paramsFileList):
        self.startJob(self.paramsFileList[len(self.jobObjList)])
      else:
        print('runJobs finished for project',self.createProjectName)
        self.createProjectName = None
        self.createProjectFinishedSignal.emit(args['projectId'])        

  def exportProject(self,projectId=None,jobIdList=None,compressedFile=None,projectName=None):
    if projectId is None and projectName is not None:
      projectId = self.dbApi.getProjectId(projectName)
    projectInfo = PROJECTSMANAGER().db().getProjectInfo(projectId=projectId)
    # If there is a limited set of jobs then find the input jobs that are not output by jobs on that list
    inputFilesList,inputFileIdList,fromJobIdList,errReport =  PROJECTSMANAGER().getJobInputFiles(projectDir=projectInfo['projectdirectory'],jobIdList=jobIdList,useDb=True)
    print('CTestDb.exportProject inputFilesList,fromJobIdList',inputFilesList,fromJobIdList)
    
    import time,os
    dbxml = os.path.join( projectInfo['projectdirectory'],'CCP4_TMP','DATABASE'+str(int(time.time()))+'.db.xml')
    
    # exportProjectXml returns list of TOP-LEVEL jobNumbers for the export
    jobNumberList,errReport = PROJECTSMANAGER().db().exportProjectXml(projectId,fileName=dbxml,recordExport=True,status='exportable',jobList=jobIdList,inputFileList=inputFileIdList,inputFileFromJobList=fromJobIdList)
    #print 'CTestDb.exportProject jobNumberList',jobNumberList
    if errReport.maxSeverity()>SEVERITY_WARNING:
      print('ERROR in exportProjectXml\n',errReport.report())

    if jobIdList is not None:
        directoriesList = []
    else:
        directoriesList = ['CCP4_IMPORTED_FILES','CCP4_PROJECT_FILES']
    from qtcore import CCP4Export
    self.exportThread = CCP4Export.ExportProjectThread(self,projectDir=projectInfo['projectdirectory'],dbxml=dbxml,target=compressedFile,jobList=jobNumberList,inputFilesList=inputFilesList,directoriesList=directoriesList)
    import functools
    self.exportThread.finished.connect(functools.partial(self.exportProjectFinished,projectId,compressedFile))
    self.exportThread.start()

  @QtCore.Slot(str,str)
  def exportProjectFinished(self,projectId,compressedFile):
    print('CTestDb.exportProjectFinished',projectId,compressedFile)
    self.exportProjectFinishedSignal.emit(projectId,compressedFile)

  def importProject(self,compressedFile,newProjectDir=None):
    from dbapi import CCP4DbApi
    xmlFile = PROJECTSMANAGER().extractDatabaseXml(compressedFile)
    dbImport = CCP4DbApi.CDbXml(db=self.dbApi,xmlFile=xmlFile,diagnostic=False)
    importProjectInfo = dbImport.loadProjectInfo()
    #print 'importProject importProjectInfo',importProjectInfo
    try:
      projectInfo = self.dbApi.getProjectInfo(projectId=dbImport.projectId)
    except:
      projectInfo = None
    #print 'importProject existingProjectInfo',projectInfo
    if projectInfo is None:
      # its a new project
      existingProject = None
      if newProjectDir is not None:
        dirName = newProjectDir
      else:
        dirName = os.path.join(self.masterDir,'DB'+self.currentDb,dbImport.projectName)
      if not os.path.exists(dirName): os.mkdir(dirName)
    else:
      dirName = projectInfo['projectdirectory']
      existingProject = dbImport.projectId
      
    dbImport.projectDirectory = dirName
    dbImport.createProject()
    dbImport.createTempTables()
    dbImport.loadTempTable()
    # If loading jobs to an existing project flag up jobs in temp tables that
    # are already in db
    if existingProject is not None:
      dbImport.setExclInTempTables()
    # Flag imported files to be imported (there is no checking yet that they exist)
    dbImport.setExclImportedFiles()

    from qtcore import CCP4Export   
    # Unpack project files from the tar file (possibly in separate thread) 
    # Pass import thread dbImport to enable query database and flagging loaded jobs/files
    self.importThread = CCP4Export.ImportProjectThread(self,projectDir=dirName,compressedFile=compressedFile,
                                                       dbImport=dbImport)
    errReport = self.importThread.run()
    
    dbImport.cleanupTempTables()
    #dbImport.listTempJobs('TempJobs after cleanup')
    #dbImport.listTempFiles('TempFiles after cleanup')
    #dbImport.listTempFileUses('TempFileUses after cleanup')
    stats = dbImport.importStats()
    for key,value in list(stats.items()):
      print('importProject stats', key,value)
    if errReport.maxSeverity()>SEVERITY_WARNING:
      dbImport.removeTempTables()
      text = 'ERRORS UNPACKING DATA FILES\n'
      for err in errReport: text = text + err['details'] + '\n'
      print(text)
      return

    dbImport.importTempTables()
    dbImport.removeTempTables()
    dbImport.db.projectReset.emit({'projectId':dbImport.projectId})
    dbImport = None

    

  def compareDb(self,projectId=None,diffFile=None,projectName=None):
    if projectId is None and projectName is not None:
      projectId = self.dbApi.getProjectId(projectName)
    self.setDatabase('1')  
    node1,jobNumberList1,errReport1  = self.dbApi.getTablesEtree(projectId=projectId)
    self.setDatabase('2')  
    node2,jobNumberList2,errReport2  = self.dbApi.getTablesEtree(projectId=projectId)
    if errReport1.maxSeverity()>SEVERITY_WARNING or errReport2.maxSeverity()>SEVERITY_WARNING:
      print('ERRORS convering Db to Etree\ndb1\n'+errReport1.report()+'\ndb2\n'+errReport2.report())
      return
    diffs = CCP4File.compareEtreeNodes(node1,node2)
    if diffFile is not None:
      from core import CCP4Utils
      CCP4Utils.saveFile(fileName=diffFile,text=diffs.report(),overwrite=1)
    for d in diffs:
      print(d)

  @QtCore.Slot(str,list)
  def test1Export(self,projectId,jobIdList=None):
    print('Exporting from ',self.currentDb,'projectId:',projectId)
    from qtcore import CCP4Export
    self.exportProjectFinished.connect(self.test1Import)
    compressedFile = os.path.join(self.masterDir,'db1_export.'+CCP4Export.COMPRESSED_SUFFIX)
    try:
      os.remove(compressedFile)
    except:
      pass
    self.exportProject(projectId,compressedFile=compressedFile,jobIdList=jobIdList)

  @QtCore.Slot(str,str)
  def test1Import(self,projectId,compressedFile):
    print('Importing to',self.currentDb,'compressedFile',compressedFile)
    self.setDatabase('2')
    self.importProject(compressedFile)

    #self.compareDb(projectId)

  def test1(self,source='simple1'):
    self.createProjectFinishedSignal.connect(self.test1Export)
    self.createProject(source)

