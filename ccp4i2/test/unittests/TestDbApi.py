def TESTSUITE():
  suite = unittest.defaultTestLoader.loadTestsFromTestCase(testSqliteDb)
  #unittest.TestLoader.testMethodPrefix = 'testx'
  #suite = unittest.defaultTestLoader.loadTestsFromTestCase(testQtDb)
  return suite

def testModule():

  CCP4File.CDataFile.BLOCK_LOAD_FILES = True

  if os.path.exists(testQtDb.TESTDBFILE): os.remove(testQtDb.TESTDBFILE)
  if os.path.exists(testSqliteDb.TESTDBFILE): os.remove(testSqliteDb.TESTDBFILE)
  for testDir in [testDb.TESTDBDIR]:
    if not os.path.exists(testDir):
      os.mkdir(testDir)
      jobsDir = os.path.join(testDir,'CCP4_JOBS')
      os.mkdir(jobsDir)

      os.mkdir(os.path.join(jobsDir,'job_1'))
      os.mkdir(os.path.join(jobsDir,'job_2'))
      os.mkdir(os.path.join(jobsDir,'job_3'))
      os.mkdir(os.path.join(jobsDir,'job_3','job_1'))
      os.mkdir(os.path.join(jobsDir,'job_3','job_2'))
      os.mkdir(os.path.join(jobsDir,'job_3','job_3'))
      os.mkdir(os.path.join(jobsDir,'job_5'))
      # Emulating imported data
      for name,job in [['starting_data.mtz','job_2'],['starting_data.seq','job_1']]:
        fileName = os.path.join(testDir,name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
        fileName = os.path.join(testDir,'CCP4_JOBS',job,name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
      # Emulating pipeline
      for name in ['model_1.pdb','model_2.pdb']:
        fileName = os.path.join(jobsDir,'job_3','job_1',name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
      for name in ['built.pdb','built.mtz']:
        fileName = os.path.join(jobsDir,'job_3','job_2',name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
      for name in ['built.pdb','built.mtz']:
        fileName = os.path.join(jobsDir,'job_3','job_3',name)
        CCP4Utils.saveFile(fileName=fileName,text='Dummy file')
      fileName = os.path.join(jobsDir,'job_5','XYZOUT.pdb')
      CCP4Utils.saveFile(fileName=fileName,text='Dummy file')

  pm = PROJECTSMANAGER()

  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)

  CCP4File.CDataFile.BLOCK_LOAD_FILES = False

class testDb(unittest.TestCase):
  TESTDBDIR = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi.test.dir')


  def test0User(self):
    uid0 = self.db.createUser('metoo')
    uid1 = self.db.getUserId('metoo')
    self.assertEqual(uid0,uid1,'Failed to create and retrieve user id')
    try:
      self.db.createUser('metoo')
    except CException as e:
      self.assertEqual( e[0]['code'], 120,'Wrong error when attempting to create user with same name')
    except:
      self.fail('Unknown error when attempting to create user with same name')
    else:
      self.fail('No error when attempting to create user with same name')


    self.db._userRole = USER_ROLE_USER
    try:
      self.db.createUser('another')
    except CException as e:
      self.assertEqual( e[0]['code'], 106,'Wrong error when attempting to create user with wrong permission')
    except:
      self.fail('Unknown error when attempting to create user with wrong permission')
    else:
      self.fail('No error when attempting to create user with wrong permission')

    self.db._userRole = USER_ROLE_OWNER
    self.db.updateUser('metoo','another')

    self.db.commit()
    self.db.close()

    self.db = CDbApi(mode=self.mode,fileName=self.TESTDBFILE,userName='me')
    uList = self.db.listUsers(True)
    print('test0User uList',uList)
    self.assertEqual(uList[1][1],'another','Failed changing user name and listing')

    info = self.db.getUserInfo(uid0)
    #print 'after getUserInfo',info
    self.assertEqual(info['username'],'another','Error calling getUserInfo')


  def test1Project(self):

    pid0 = self.db.createProject('myproject',projectDirectory=self.TESTDBDIR)
    pid1 = self.db.getProjectId('myproject')
    self.assertEqual(pid0,pid1,'Error creating project and retrieving project ID')

    baby = self.db.createProject('babyproject',userName='me',parentProjectId=pid1,projectDirectory=os.path.join(self.TESTDBDIR,'baby'))
    self.db.commit()
    self.db.close()

    self.db = CDbApi(mode=self.mode,fileName=self.TESTDBFILE,userName='me')
    #print 'test1Project _userName',self.db._userName
    pList = self.db.listProjects(True)
    self.assertEqual(pList[1][3],pid0,'Failed creating child project and listing')

    self.db.updateProject(baby,key='projectDirectory',value=os.path.join(self.TESTDBDIR,'baby2'))
    pList = self.db.listProjects(True)
    self.assertEqual(os.path.split(pList[1][2])[-1],'baby2','Failed updating project directory')


    '''
    # Recent projects concept currently dropped
    #self.db.listContexts(True)
    recent = self.db.getRecentProjects()
    #print 'test1Project getRecentProjects',recent
    self.assertEqual(recent,[[1,'myproject']],'Failed retrieving recent projects')
    '''


  def test2Aliases(self):
     self.db.setDirectoryAlias('CCP4I2_TOP',CCP4Utils.getCCP4I2Dir())
     self.db.setDirectoryAlias('HOME',CCP4Utils.getHOME())
     self.db.commit()

     aliasList = self.db.listDirectoryAliases(toTerm=True)
     self.assertEqual(len(aliasList),2,'Alias list wrong length')
     self.assertTrue(('CCP4I2_TOP' in aliasList),'Alias list does not have CCP4I2_TOP')



  def test3Jobs(self):
    #pid = self.db.createProject('myproject',projectDirectory=self.TESTDBDIR)
    pid = self.db.getProjectId('myproject')
    print('test3Jobs pid',pid)
    self.db.listProjects(True)

    # Two jobs loading contents and exptal data
    jid1 =  self.db.createJob(pid,'contents',status=JOB_STATUS_FINISHED)
    fo1 = CCP4ModelData.CSeqDataFile(project='myproject',relPath='CCP4_JOBS/job_1',baseName='starting_data.seq')
    f1=self.db.createFile(jobId=jid1,fileObject=fo1,sourceFileName=os.path.join(CCP4Utils.getTestTmpDir(),'starting_data.seq'))

    jid2 =  self.db.createJob(pid,'import_xtal',status=JOB_STATUS_FINISHED)
    fo2 = CCP4XtalData.CMtzDataFile(project='myproject',relPath='CCP4_JOBS/job_2',baseName='starting_data.mtz')
    f2=self.db.createFile(jobId=jid2,fileObject=fo2,sourceFileName=os.path.join(CCP4Utils.getTestTmpDir(),'starting_data.mtz'))


    # An 'auto mr' pipeline
    jid3 =  self.db.createJob(pid,'auto_mr')
    self.db.createFileUse(jobId=jid3,fileId=f1)
    self.db.createFileUse(jobId=jid3,fileId=f2)
    jobNo = self.db.getJobInfo(jid3,'jobNumber')
    self.assertEqual(jobNo,'3','Error retrieving jobNumber')

    # mr run creating two pdbs files
    jid4 =  self.db.createJob(pid,'mr',parentJobId=jid3,jobTitle='testing a parentJobId and jobTitle')
    self.db.createFileUse(jobId=jid4,fileId=f1)
    self.db.createFileUse(jobId=jid4,fileId=f2)
    fo3 = CCP4ModelData.CPdbDataFile(project='myproject',relPath='CCP4_JOBS/job_3/job_1',baseName='model_1.pdb')
    fo4 = CCP4ModelData.CPdbDataFile(project='myproject',relPath='CCP4_JOBS/job_3/job_1',baseName='model_2.pdb')
    f3=self.db.createFile(jobId=jid4,fileObject=fo3)
    f4=self.db.createFile(jobId=jid4,fileObject=fo4)
    self.assertEqual(self.db.getJobInfo(jid4,'ParentJobId'),jid3,'Error setting and retrieving parent job number')
    self.db.updateJobStatus(jid4,JOB_STATUS_FINISHED)

    # Two runs of model building
    jid5 =  self.db.createJob(pid,'build',parentJobId=jid3,jobTitle='build')
    self.db.createFileUse(jobId=jid5,fileId=f2)
    self.db.createFileUse(jobId=jid5,fileId=f3)
    fo5 = CCP4ModelData.CPdbDataFile(project='myproject',baseName='built.pdb',relPath='CCP4_JOBS/job_3/job_2')
    f5=self.db.createFile(jobId=jid5,fileObject=fo5,fileTypeName='chemical/x-pdb')
    self.db.updateJobStatus(jid5,JOB_STATUS_FINISHED)

    fileIds = self.db.getJobFiles(jobId=jid5,role=FILE_ROLE_IN,mode='fileId')
    self.assertTrue(fileIds==[f2,f3] or fileIds==[f3,f2],'Failed restoring input file ids')
    fileNames = self.db.getJobFiles(jobId=jid5,role=FILE_ROLE_OUT,mode='fileName')
    self.assertEqual(fileNames,['built.pdb'],'Failed restoring input file ids')

    # Lets say second model building failed
    jid6 =  self.db.createJob(pid,'build',parentJobId=jid3,jobTitle='build')
    self.db.createFileUse(jobId=jid6,fileId=f2)
    self.db.createFileUse(jobId=jid6,fileId=f4)
    #f6=self.db.createFile(jobId=jid6,fileName='job_6_built.pdb',relPath='CCP4_JOBS/job_6',fileTypeName='text/pdb')
    self.db.updateJobStatus(jid6,JOB_STATUS_FAILED)

    # Finish off the pipeline

    self.db.updateJobStatus(status='Finished',projectName='myproject',jobNumber='3')
    self.db.createFileUse(jobId=jid3,fileId=f5,role=FILE_ROLE_OUT)
    info = self.db.getJobInfo(jid3)
    self.assertEqual(info['status'],'Finished','Error changing job status')
    self.assertEqual(info['taskname'],'auto_mr','Error retrieving taskName')
    fileTypeList = self.db.getJobFiles(jobId=jid3,role=FILE_ROLE_OUT,mode='fileType')
    # Output should be one PDB file
    #print 'test3Jobs fileTypeList',fileTypeList
    self.assertEqual(fileTypeList,['chemical/x-pdb'])

    fileList = self.db.getJobFiles(jobId=jid5,role=FILE_ROLE_OUT,mode='fullPath')
    self.assertEqual(fileList,[os.path.join(self.TESTDBDIR,'CCP4_JOBS','job_3','job_2','built.pdb')],'Failed running getJobFiles')
    fileId,jobId = self.db.matchFileName(fileName = os.path.join(self.TESTDBDIR,'CCP4_JOBS','job_3','job_2','built.pdb'))
    #print 'test3Jobs matchFileName',fileId,jobId ,'f5',f5
    self.assertEqual(fileId,f5,'Error using matchFileName')

    info = self.db.getProjectJobListInfo(projectName='myproject')
    #print 'testJob getJobListInfo',info


  def test4Task(self):
    pid = self.db.getProjectId('myproject')
    jid = self.db.createJob(pid,'pdbset')
    self.db.commit()


  def test5glean(self):
    # This will e jobNumber=5
    pid = self.db.getProjectId('myproject')
    jid7 =  self.db.createJob(projectId=pid,taskName='pdbset',jobTitle='test5glean')
    shutil.copyfile(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','test_dbapi.params.xml'),PROJECTSMANAGER().makeFileName(jobId=jid7,mode='PARAMS'))
    self.db.updateJobStatus(jid7,status='Finished')
    #self.db.setDiagnostic(True)
    #print 'test5glean jid7',jid7,PROJECTSMANAGER().makeFileName(jobId=jid7,mode='PARAMS')
    rv = self.db.gleanJobFiles(jobId=jid7)
    #print 'from gleanJobFiles',rv.report()

    fileList = self.db.getJobFiles(jobId=jid7,role=FILE_ROLE_IN,mode='fullPath')
    #print 'test5glean fileList',fileList
    self.assertEqual(len(fileList),1,'Error gleaning files from params - wrong number of input files')
    self.assertEqual(os.path.split(fileList[0])[1],'model_1.pdb','Error gleaning files from params - wrong file name')

    fileList = self.db.getJobFiles(jobId=jid7,role=FILE_ROLE_OUT)
    self.assertEqual(len(fileList),1,'Error gleaning files from params - wrong number of output files')
    filePath = self.db.getFullPath(fileId=fileList[0])
    self.assertEqual(os.path.split(filePath)[1],'XYZOUT.pdb','Error gleaning files from params - wrong output file name')

  def test6glean(self):
    pid = self.db.getProjectId('myproject')
    pf =  os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','test_dbapi_2.params.xml')
    jid8 =  self.db.createJob(pid,'newProject',jobTitle='test6glean')
    shutil.copyfile(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','test_dbapi_2.params.xml'),PROJECTSMANAGER().makeFileName(jobId=jid8,mode='PARAMS'))

    PROJECTSMANAGER().importFiles(jobId=jid8)
    self.db.updateJobStatus(jid8,status='Finished')
    rv = self.db.gleanJobFiles(jobId=jid8)

    containsSe = self.db.getXDataByJobContext(contextJobId=jid8,dataClass='CContainsSeMet',projectId=pid)
    #print 'CCP4DbApi.test6glean containsSe',containsSe

    fileName=os.path.join(CCP4Utils.getTestTmpDir(),'export_db.xml')
    if os.path.exists(fileName): os.remove(fileName)
    self.db.exportProjectXml(projectId=pid,fileName=fileName)
    self.assertTrue(os.path.exists(fileName),'Export xml file not created')

  def test7comment(self):

    jid3 =  self.db.getJobId(projectName='myproject',jobNumber='3')
    cid = self.db.createComment(jobId=jid3,comment='Well whatever really')

    commentList = self.db.getComments(jobId=jid3)
    #print 'test7comment commentList',commentList
    self.assertEqual('Well whatever really',commentList[0][3],'Error creating and retrieving comment')


class testSqliteDb(testDb):

  TESTDBFILE = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi.testsqlite.db')

  def setUp(self):

    self.db = CDbApi(mode='sqlite',userName='me',userPassword='foo',fileName=testSqliteDb.TESTDBFILE)
    self.mode = 'sqlite'
    PROJECTSMANAGER().setDatabase(self.db)

class testQtDb(testDb):

  TESTDBFILE = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi.testqt.db')

  def setUp(self):
    self.db = CDbApi(mode='qt_sqlite',fileName=testQtDb.TESTDBFILE,userName='me',userPassword='foo')
    self.mode = 'qt_sqlite'
    PROJECTSMANAGER().setDatabase(self.db)
