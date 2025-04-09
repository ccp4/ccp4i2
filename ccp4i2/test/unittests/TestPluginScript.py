# NB Need to run wrapper/mtzdump pipelines/demo_copycell tests to test running process


class testCPluginScript(unittest.TestCase):

    def setUp(self):
        self.app = QTAPPLICATION()
        self.script = CPluginScript(name='test_CPluginScript', parent=self.app, workDirectory='/foo/bar')
        self.hklin = CCP4File.CDataFile(project='CCP4I2_TOP', relPath='test/data', baseName='1df7.pdb')

    def test1(self):
        self.script.appendCommandLine(['HKLIN', str(self.hklin)])
        self.script.appendCommandLine(['test', 1])
        filePath = os.path.join(CCP4Utils.getCCP4I2Dir(), 'test', 'data', '1df7.pdb')
        #print 'commandLine',self.script.commandLine
        self.assertEqual(self.script.commandLine,['HKLIN', filePath, 'test', '1'], 'Error in appendCommandLine')

    def test2(self):
        filePath = os.path.join(CCP4Utils.getCCP4I2Dir(), 'test', 'data', '1df7.pdb')
        comFile = self.script.makeFileName('COM')
        comText = '''TEST 1
TEST 2
TEST 3
HKLIN ''' + filePath + '\n'
        if os.path.exists(comFile):
            os.remove(comFile)
        self.script.appendCommandScript(['TEST 1', 'TEST 2'])
        self.script.appendCommandScript(['TEST', '3'], oneLine=True)
        self.script.appendCommandScript(['HKLIN', self.hklin], oneLine=True)
        self.script.writeCommandFile()
        self.assertTrue(os.path.exists(comFile), 'Failed to write command file')
        text = CCP4Utils.readFile(comFile)
        #print 'comFile',text
        self.assertEqual(text, comText, 'Wrong text in com file')

    def test3(self):
        s = self.script.makePluginObject('buccaneer')
        self.assertTrue(isinstance(s, CPluginScript), 'Failed to run makePluginObject')
        self.assertEqual(str(s.objectName()), 'test_CPluginScript_1', 'makePluginObject created objected has wrong name')
        self.assertEqual(s.workDirectory, os.path.join(self.script.workDirectory, 'job_1'), 'makePluginObject created objected has wrong work directory')

def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testCPluginScript)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
