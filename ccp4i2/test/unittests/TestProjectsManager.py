def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testManager)
    #suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testProject))
    #suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testJob))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testJob(unittest.TestCase):

    def test1(self):
        j = CJob(jobId=99,taskName='testing')  # KJS : Looks like this may not work ...
        j.inputFiles.HKLIN = {  'baseName' : 'foo.mtz', 'project' : 'myProject' }
        j.outputFiles.set({'HKLOUT' : { 'baseName' : 'foo_99.mtz', 'project' : 'myProject' }} )
        tree = j.getEtree()
        ET.indent(tree)
        text = ET.tostring(tree)
        #print text
        k = CJob(jobId=99)
        k.setEtree(tree)
        self.assertEqual(str(k.jobId),'job_99','Wrong job id')
        self.assertEqual(k.inputFiles.HKLIN.get(),{ 'baseName' : 'foo.mtz', 'project' : 'myProject' , 'relPath' : None },'Wrong inputFile')
        self.assertEqual(k.outputFiles.HKLOUT.get(),{ 'baseName' : 'foo_99.mtz', 'project' : 'myProject' , 'relPath' : None },'Wrong outputFile')

class testProject(unittest.TestCase):

    def test1(self):
        p = CProject(name='myproject')     # KJS : More problems ?
        jobId = p.newJob(taskName='testing',title='Just testing projects')
        print('testProject.test1',jobId,list(p.jobList.keys()))
        print(p.jobList[jobId].taskName)
        self.assertEqual(p.jobList[jobId].taskName,'testing','Error in CProject.newJob')

    def test2(self):
        p = CProjectDirectories()
        p.setDirectory('/Users/me/myProject')
        self.assertEqual(p.getDirectory(),'/Users/me/myProject','Failed setting CProjectDirectories')
        user = CCP4Utils.getUserId() + '@' + platform.node()
        self.assertEqual(p[user],'/Users/me/myProject','Error getting currentUserAddress')
        p = CProject(name='myproject')
        p.directories.setDirectory('/Users/me/myProject')
        self.assertEqual(p.directories.getDirectory(),'/Users/me/myProject','Failed setting CProject.directories')
        jobId = p.newJob(taskName='testing',title='Just testing projects')
        fileName = p.makeFileName(jobId=jobId,mode='LOG')
        self.assertEqual(fileName,'/Users/me/myProject/'+jobId+'_testing.log')
        tree = p.getEtree()
        q = CProject(name='myproject')
        q.setEtree(tree)
        jobId = p.newJob(taskName='anotherTest',title='Another test job')
        self.assertEqual(jobId,'job_2','Error inn XML cycle returning jobId')


class testManager(unittest.TestCase):

    def setUp(self):
        self.testDir = os.path.join(getHOME(),'CCP4I2_TEST')
        if not os.path.exists(self.testDir):
            os.mkdir(self.testDir)
        dbFile = os.path.join(self.testDir,'test_database.sql')
        if os.path.exists(dbFile): os.remove(dbFile)
        self.pm = CProjectsManager(name='TESTPROJECTMANAGER',databaseFile=dbFile)
        self.pm.createProject('myproject',os.path.join(self.testDir,'myproject'))
        self.pm.createProject('anotherproject',os.path.join(self.testDir,'anotherproject'))

    def test1(self):
        pDir = self.pm.getProjectDirectory(projectName='myproject')
        self.assertTrue(pDir,os.path.join(self.testDir,'myproject'),'Error testing project directory')

    def test2(self):
        try:
            self.pm.createProject('myproject',os.path.join(self.testDir,'whatever'))
        except CException as e:
            self.assertEqual(e[0]['code'] , 110, 'Wrong error code (not 110) creating project with same name')
        except:
            self.fail('Failed unknown reason attempting to create project with same name')
        try:
            self.pm.createProject('newproject',os.path.join(self.testDir,'myproject'))
        except CException as e:
            self.assertEqual(e[0]['code'] , 117, 'Wrong error code (not 117) creating project with same dir')
        except:
            self.fail('Failed with unknown reason attempting to create project with same dir')
