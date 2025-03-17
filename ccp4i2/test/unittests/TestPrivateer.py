class testPrivateer(unittest.TestCase):

  def test1(self):
    # Test creation of log file using ../test_data/test1.params.xml input
    workDirectory = CCP4Utils.getTestTmpDir()
    logFile = os.path.join(workDirectory,'privateer_test1.log')
    # Delete any existing log file
    if os.path.exists(logFile): os.remove(logFile)
    self.wrapper = privateer(name='privateer_test1',workDirectory=workDirectory)
    self.wrapper.container.loadDataFromXml(os.path.join(getCCP4I2Dir(),'wrappers','privateer','test_data','test1.params.xml'))
    self.wrapper.setWaitForFinished(1000000)
    pid = self.wrapper.process()
    self.wrapper.setWaitForFinished(-1)
    if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())
    #self.assertTrue(os.path.exists(logFile),'No log file found')


def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testPrivateer)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
