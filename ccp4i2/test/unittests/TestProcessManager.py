# Python shell test:
# p = CCP4ProcessManager.CProcessManager(); p.startProcess(command='mtzdump',argList=['HKLIN','/Users/lizp/Desktop/test_data/rnase25_phases.mtz'],inputText='HEADER\nGO\n')

def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testProcessManager)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
  

# This test class can run without unittest framework - so the PROCESSMANAGER
# can work 'naturally' without the waitForFinish


#*  unittest
class testProcessManager(unittest.TestCase):
#* hand testing
#class testProcessManager():
#  def __init__(self):
#    self.isUnitTest=False
#    self.setUp()
# *
    def setUp(self):
        if not hasattr(self, 'isUnitTest'):
            self.isUnitTest = True
        self.pdbFile = os.path.join(CCP4Utils.getCCP4I2Dir(), 'test', 'data', '1df7.pdb')
        if self.isUnitTest:
            PROCESSMANAGER().setWaitForFinished(1000)

    def tearDown(self):
        PROCESSMANAGER().setWaitForFinished(-1)

    def evalCRYST(self, pdbFile=None):
        text = CCP4Utils.readFile(pdbFile)
        for line in text.split('\n'):
            if line[0:6] == 'CRYST1':
                words = line.split()
                cryst = []
                for i in range(1, 6):
                    cryst.append(float(words[i]))
                #print 'cryst',cryst
                if 49.999 < cryst[0] < 50.001:
                    return True
        return False

    def printReview(self, processID):
        print(' ')
        print('processID exitStatus', processID, PROCESSMANAGER().getJobData(processID, 'exitStatus'))
        print('Error:', PROCESSMANAGER().getJobData(processID, 'processError'))
        print('Is PDB file created:', self.pdbOut, str(os.path.exists(self.pdbOut)))
        if os.path.exists(self.pdbOut):
            print('Output PDB contains correct data',str(self.evalCRYST(self.pdbOut)))
        print(' ')

    def test1(self):
        # This ine should run OK
        self.pdbOut = CCP4Utils.makeTmpFile(name='testProcessManager_test1')
        print('test1 pdbOut', self.pdbOut)
        if self.isUnitTest:
            handler = None
        else:
            handler = [self.review1, {}]
        processID = PROCESSMANAGER().startProcess(command='pdbset', args=['XYZIN', self.pdbFile, 'XYZOUT', self.pdbOut],
                                                              inputText='CELL 50.0 50.0 70.0 90.0 90.0 90.0\nEND\n', handler=handler)
        if self.isUnitTest:
            self.review1(processID)

    def review1(self, processID=None):
        if self.isUnitTest:
            self.assertTrue(os.path.exists(self.pdbOut), 'No output PDB file created')
            self.assertTrue(self.evalCRYST(self.pdbOut), 'Output PDB file does not contain correct data')
        else:
            self.printReview(processID)

    def test2(self):
        # Input file does not exist
        self.pdbOut = CCP4Utils.makeTmpFile(name='testProcessManager_test1')
        if self.isUnitTest:
            handler = None
        else:
            handler = [self.review2, {}]
        processID = PROCESSMANAGER().startProcess(command='pdbset', args=['XYZIN', 'foobar', 'XYZOUT', self.pdbOut],
                                  inputText='CELL 50.0 50.0 70.0 90.0 90.0 90.0\nEND\n', handler=handler)
        if self.isUnitTest:
            self.review2(processID)

    def review2(self,processID=None):
        if self.isUnitTest:
            exitCode = PROCESSMANAGER().getJobData(processID, 'exitStatus')
            error = PROCESSMANAGER().getJobData(processID, 'processError')
            self.assertEqual(exitCode, 1, 'Wrong exit code when bad input filename')
            self.assertEqual(error.count('No such file or directory'), 1, 'Wrong error message when bad input filename')
        else:
            self.printReview(processID)

    def test3(self):
        # Calling a non-existant executable
        self.pdbOut =  CCP4Utils.makeTmpFile(name='testProcessManager_test3')
        if self.isUnitTest:
            handler = None
        else:
            handler = [self.review2,{}]
        processID = PROCESSMANAGER().startProcess(command='pdbset-no', args=['XYZIN', 'foobar', 'XYZOUT', self.pdbOut],
                                                              inputText='CELL 50.0 50.0 70.0 90.0 90.0 90.0\nEND\n', handler=handler)
        if self.isUnitTest:
            self.review3(processID)

    def test4(self):
        # Bad input (CELL incomplete) - expect process to return with error code
        self.pdbOut = CCP4Utils.makeTmpFile(name='testProcessManager_test3')
        if self.isUnitTest:
            handler = None
        else:
            handler = [self.review2, {}]
        processID = PROCESSMANAGER().startProcess(command='pdbset', args=['XYZIN','foobar','XYZOUT',self.pdbOut],
                                                              inputText='CELL 50.0 \nEND\n', handler=handler)
        if self.isUnitTest:
            self.review4(processID)
