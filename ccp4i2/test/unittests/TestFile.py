def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testProject)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testFilePath))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testDataFile))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testI2XmlDataFile))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testXmlDataFile))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testProject(unittest.TestCase):

    def testCProjectName1(self):
        p = CProjectName()
        p.set('CCP4I2_TEST')
        self.assertEqual(p,'CCP4I2_TEST','Failed to set good project name')
    """
    def testCProjectName2(self):
        p = CProjectName()
        try:
          p.set('rubbish')
        except CException as e:
          self.assertEqual(len(e),1,'Unexpected exception length in setting bad CProjectName')
          self.assertEqual(e[0]['code'],104,'Unexpected exception in setting bad CProjectName')
        except:
          self.fail('Unexpected exception in setting bad CProjectName')
        else:
          self.fail('No exception in setting bad CProjectName')

    def testCProjectName3(self):
        PROJECTSMANAGER().createProject('DUMMY_CCP4I2_TEST',os.path.join(getHOME(),'DUMMY_CCP4I2_TEST'))
        shutil.rmtree(os.path.join(getHOME(),'DUMMY_CCP4I2_TEST'))

        p = CProjectName()
        try:
          p.set('DUMMY_CCP4I2_TEST')
        except CException as e:
          self.assertEqual(len(e),1,'Unexpected exception length in setting project with nonexistant directory')
          self.assertEqual(e[0]['code'],106,'Unexpected exception in setting project with nonexistant directory')
        except:
          self.fail('Unexpected exception in setting bad project with nonexistant directory')
        else:
          self.fail('No exception in setting project with nonexistant directory')
        PROJECTSMANAGER().deleteProject('DUMMY_CCP4I2_TEST')
    """

class testFilePath(unittest.TestCase):

    def testFilePath1(self):
        f = CFilePath('foo/bar/fo_ehue.tmp')
        self.assertEqual(f, 'foo/bar/fo_ehue.tmp', 'Error setting CFilePath')

    def testFilePath2(self):
        try:
            f = CFilePath('foo/bar/fo_e*ue.tmp', allowedCharactersMode=CFilePath.ALLOWED_CHARACTERS_FAIL)
        except CException as e:
            self.assertEqual(len(e), 1, 'Unexpected exception length in setting bad CFilePath')
            self.assertEqual(e[0]['code'], 101, 'Unexpected exception in setting project bad CFilePath')
        except:
            self.fail('Unexpected exception in setting bad CFilePath')
        else:
            self.fail('No exception in setting bad CFilePath')

    def testFilePath3(self):
        try:
            f = CFilePath()
            f.set(f.fix('foo/bar/fo_e*ue.tmp'))
        except:
            self.fail('Unexpected exception in fixing bad CFilePath')
        self.assertEqual(f, 'foo/bar/fo_e_ue.tmp', 'Error setting CFilePath')

class testDataFile(unittest.TestCase):

    def setUp(self):
        PROJECTSMANAGER().makeDefaultProject('CCP4I2_TEST')

    def testDataFile0(self):
        f = CDataFile()
        self.assertEqual(f.baseName.get(), None,'Initialising CDataFile does not have baseName as None')
        self.assertEqual(f.fullPath.get(), None,'Initialising CDataFile does not have fullPath as None')

    def testDataFile1(self):
        f = CDataFile(projectName='CCP4I2_TEST', relPath='foo/bar', baseName='myfile.txt')
        #print 'testDataFile1',f.fullPath,f.relPath,f.project.directory()
        #print 'testDataFile1',PROJECTSMANAGER().getProjectDirectory(projectName='CCP4I2_TEST')
        #print 'testDataFile1',PROJECTSMANAGER().db().listProjects(toTerm=True)
        self.assertEqual(f.fullPath,os.path.join(PROJECTSMANAGER().getProjectDirectory(projectName='CCP4I2_TEST'),'foo/bar/myfile.txt'))

    def testDataFile2(self):
        f = CDataFile(mustExist=True)
        projectId=PROJECTSMANAGER().db().getProjectId(projectName = 'CCP4I2_TEST')
        try:
            f.set(projectId=projectId,relPath = 'foo/bar',baseName = 'myfile.txt')
        except CException as e:
            self.assertEqual(len(e),1,'Unexpected exception length in setting non-existant CDataFile')
            self.assertEqual(e[0]['code'],101,'Unexpected exception in setting non-existant CDataFile')
        except:
            self.fail('Unexpected exception in setting non-existant CDataFile')
        else:
            self.fail('No exception in setting non-existant CDataFile')

    def testDataFile3(self):
        f = CDataFile()
        f.fullPath = os.path.join(PROJECTSMANAGER().getProjectDirectory(projectName='CCP4I2_TEST'),'foo/bar/myfile.txt')
        projectName = PROJECTSMANAGER().db().getProjectInfo(str(f.project),'projectName')
        self.assertEqual(projectName,'CCP4I2_TEST','CDataFile setting full path gives wrong project')
        self.assertEqual(f.relPath,'foo/bar','CDataFile setting full path gives wrong relPath')
        self.assertEqual(f.baseName,'myfile.txt','CDataFile setting full path gives wrong baseName')

    def testDataFile4(self):
        f = CDataFile(fullPath = os.path.join(PROJECTSMANAGER().getProjectDirectory(projectName='CCP4I2_TEST'),'foo/bar/myfile.txt'))
        projectName = PROJECTSMANAGER().db().getProjectInfo(str(f.project),'projectName')
        self.assertEqual(projectName,'CCP4I2_TEST','CDataFile setting full path gives wrong project')
        self.assertEqual(f.relPath,'foo/bar','CDataFile setting full path gives wrong relPath')
        self.assertEqual(f.baseName,'myfile.txt','CDataFile setting full path gives wrong baseName')

    def testDataFile5(self):
        bindir = os.path.join(os.environ['CCP4I2_TOP'],'bin')
        f = CDataFile(relPath=bindir,baseName='browser')
        #print 'testDataFile5',bindir,f.fullPath
        h = CFilePath(os.path.join(bindir,'browser'))
        self.assertTrue(f.samefile(os.path.join(bindir,'browser')),'Failed CDataFile.samefile(Python string)')
        self.assertTrue(f.samefile(h),'Failed CDataFile.samefile(CFilePath)')


class testI2XmlDataFile(unittest.TestCase):

    def test1(self):
        c = CI2XmlDataFile( projectName='CCP4I2_TOP',relPath='test/data',baseName='pdbset.def.xml')
        self.assertEqual(c.header.pluginTitle,'PDBSet','Error loading header of I2XmlDataFile')

    def test2(self):
        c = CI2XmlDataFile(projectName='CCP4I2_TOP',relPath='test/data',baseName='pdbset.def.xml')
        e = c.getBodyEtree()
        self.assertEqual(str(e.find('container').get('id')),'inputData','Error in getBodyEtree')

    def test3(self):
        c = CI2XmlDataFile(projectName='CCP4I2_TEST',baseName='testXmlDataFile.xml')
        if c.fullPath.exists() : os.remove(c.fullPath.get())
        ele = etree.Element(CI2XmlDataFile.BODY_TAG)
        c.saveFile(bodyEtree=ele)
        self.assertTrue(os.path.exists(c.fullPath.get()),'No file written by saveFile')


class testXmlDataFile(unittest.TestCase):

    def test1(self):
        c = CXmlDataFile(projectName='CCP4I2_TOP',relPath='test/data',baseName='refmac_as2m1_n.xml')
        e = c.loadFile()
        self.assertEqual(e.tag,'REFMAC','Failed to load refmac xml to CXmlDataFile')

        f = CDataFile(projectName='CCP4I2_TEST',baseName='refmac_as2m1_n.refmac.xml')
        if f.fullPath.exists():  os.remove(str(f.fullPath))
        d = c.makeI2XmlDataFile(fileName=f,function='REFMAC',pluginVersion='5.1.1')
        self.assertTrue( f.fullPath.exists(),'Failed to create I2 version of refmac xml file')


def testTiming():
    f = CDataFile()
    start = time.perf_counter()
    for n in range(0,1000):
        v = f.validity('/foo')
        print('testTiming',time.perf_counter()-start)
