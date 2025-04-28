def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testMakeProjectDbXml)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testJobDbBackup))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)


class testMakeProjectDbXml(unittest.TestCase):

    def setUp(self):
        self.TESTDBFILE = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi.testsqlite.db')
        if os.path.exists(self.TESTDBFILE): os.remove(self.TESTDBFILE)
        # Create a temporary database (only used for getFileTypeId)
        self.db = CCP4DbApi.CDbApi(mode='sqlite',userName='me',userPassword='foo',fileName=self.TESTDBFILE)
        self.mode = 'sqlite'
        PROJECTSMANAGER().setDatabase(self.db)
        self.makeDb = CMakeProjectDbXml(db=self.db)

    def test1(self):
        self.makeDb.loadProject(projectDir = '/Users/lizp/Desktop/test_projects/tt3')
        self.makeDb.saveXmlFile(xmlFile= '/Users/lizp/Desktop/make_tt3.ccp4db.xml')


class testJobDbBackup(unittest.TestCase):

    def setUp(self):
        self.TESTDBFILE = os.path.join(CCP4Utils.getTestTmpDir(),'CCP4DbApi.testsqlite.db')
        self.db = CCP4DbApi.CDbApi(mode='sqlite',userName='me',userPassword='foo',fileName=self.TESTDBFILE)
        print('testJobDbBackup.setUp..')
        self.db.listProjects()

    def test1(self):
        self.backup = CJobDbBackup(jobId=1,jobNumber='1',fileName='/Users/lizp/Desktop/testJobDbBackup.ccp4db.xml')
        self.backup.updateJob('preceecingjobid',999)
        self.backup.editImportFile(fileId=111,importId=1,sourceFileName='/whereever/whatever.pdb',fileName='HKLIN.pdb',annotation='testing testing',creationTime=9999999.9)
        self.backup.editExportFile(fileId=111,exportId=1,exportFileName='/whereever/whatever.pdb',creationTime=9999999.9)
        self.backup.editComment(commentId=222,userName='me',timeOfComment=9999999.9,comment='Try this')
        self.backup.save()
