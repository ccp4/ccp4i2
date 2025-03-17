def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testEtreeTools)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testFileOpen))
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testBackup))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testEtreeTools(unittest.TestCase):

    def test1(self):
        fileName = os.path.join(getCCP4I2Dir(),'test','data','test_job_104.def.xml')
        root = openFileToEtree(fileName=fileName)
        self.assertEqual(root.tag,'ccp4i2','Failed to read first item in eTree')
 
    def test2(self):
        fileName = os.path.join(getTestTmpDir(),'CCP4Utils_test.def.xml')
        x = CCP4Data.CFloat(42,default=0)
        ele = x.getEtree()
        saveEtreeToFile(tree=ele,fileName=fileName)
        if not os.path.exists(fileName):
            self.fail('Failed to write eTree file')

class testFileOpen(unittest.TestCase):

    def test1(self):
        fileName = os.path.join(getCCP4I2Dir(),'test','data','test_job_104.def.xml')
        text = readFile(fileName)
        textLines = text.split('\n')
        self.assertEqual(textLines[0],"<?xml version='1.0'?>",'Error reading file with readFile')

    def test2(self):
        text = 'The quick brown fox jumped over the lazy dog'
        fileName = os.path.join(getTMP(),'fox.dog')
        saveFile(fileName,text)
        self.assertTrue(os.path.exists(fileName),'Failed to save file with saveFile')

    def test3(self):
        text = 'The quick brown fox jumped over the lazy dog'
        fileName = '/foo/bar/fox.dog'
        try:
            saveFile(fileName,text)
        except CException as e:
            #print 'testFileOpen.test3',e.report()
            self.assertEqual(e[0]['code'],110,'saveFile gives unexpected CException code')
        except:
            self.fail('saveFile gives unexpected Python exception')

class testBackup(unittest.TestCase):
    def test1(self):
        tempdir = getTMP()
        fileName = os.path.join(tempdir,'fox.PARAMS.xml')
        try:
            os.remove(os.path.join(tempdir,'fox.backup_1.PARAMS.xml'))
        except:
            pass
        saveFile(fileName,'whatever')
        newFile = backupFile(fileName=fileName,delete=True)
        self.assertEqual(newFile,os.path.join(tempdir,'fox.backup_1.PARAMS.xml'),'backupFile gives wrong file name')   
        if not os.path.exists(newFile):
            self.fail('backupFile does not create new file')
        if os.path.exists(fileName):
            self.fail('backupFile old file still exists')
