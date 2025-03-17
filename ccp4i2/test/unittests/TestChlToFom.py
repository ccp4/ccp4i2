
class testchltofom(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    self.app = QTAPPLICATION()

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     workDirectory = CCP4Utils.getTestTmpDir()
     logFile = os.path.join(workDirectory,'chltofom_test1.log')
     # Delete any existing log file
     if os.path.exists(logFile): os.remove(logFile)
     outFile = os.path.join(workDirectory,'chltofom_test1.mtz')
     print('testchltofom outFile',outFile)
     if os.path.exists(outFile): os.remove(outFile)
     wrapper = chltofom(parent=QTAPPLICATION(),name='chltofom_test1')
     wrapper.container.inputData.HKLIN.setFullPath(os.path.join(CCP4Utils.getCCP4I2Dir(),'test','data','rnase25_mini_HL.mtz'))
     wrapper.container.outputData.HKLOUT.setFullPath(outFile)
     wrapper.container.controlParameters.OUTPUTMINIMTZ.set(True)
     pid = wrapper.process()
     self.assertTrue(os.path.exists(  outFile),'No output file from chltofom_test1')                                     
     wrapper.container.outputData.HKLOUT.loadFile()
     columns = wrapper.container.outputData.HKLOUT.fileContent.getListOfColumns()
     print('chltofom.processOutputFiles',columns)
     self.assertEqual(len(columns),2,'Output from chltofom_test1 has wrong number of columns')
     self.assertTrue(columns[0].columnLabel.__str__() in ['PHI','FOM'] and columns[1].columnLabel.__str__() in ['PHI','FOM'],'Output from chltofom_test1 has wrong column labels')


def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testchltofom)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
