
class testcad_copy_column(unittest.TestCase):

   def setUp(self):
    self.app = QTAPPLICATION()
    # make all background jobs wait for completion
    # this is essential for unittest to work
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     workDirectory = CCP4Utils.getTestTmpDir()
     logFile = os.path.join(workDirectory,'cad_copy_column_test1.log')
     # Delete any existing log file
     if os.path.exists(logFile): os.remove(logFile)

     self.wrapper = cad_copy_column(parent=QTAPPLICATION(),name='cad_copy_column_test1',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','cad_copy_column','test_data','test1.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())
     #self.assertTrue(os.path.exists(logFile),'No log file found')
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testcad_copy_column)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
