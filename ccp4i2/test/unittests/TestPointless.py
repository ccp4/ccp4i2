
class testpointless(unittest.TestCase):

   def setUp(self):
    self.app = QTAPPLICATION()
    # make all background jobs wait for completion
    # this is essential for unittest to work
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     workDirectory = os.path.join(CCP4Utils.getTestTmpDir(),'test1')
     if not os.path.exists(workDirectory): os.mkdir(workDirectory)

     self.wrapper = pointless(parent=QTAPPLICATION(),name='test1',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','pointless','test_data','test1.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())

   def test_2(self):
     workDirectory = os.path.join(CCP4Utils.getTestTmpDir(),'test2')
     if not os.path.exists(workDirectory): os.mkdir(workDirectory)

     self.wrapper = pointless(parent=QTAPPLICATION(),name='test2',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','pointless','test_data','test2.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())

   def test_3(self):
     workDirectory = os.path.join(CCP4Utils.getTestTmpDir(),'test3')
     if not os.path.exists(workDirectory): os.mkdir(workDirectory)

     self.wrapper = pointless(parent=QTAPPLICATION(),name='test3',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','pointless','test_data','test3.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())

   def test_4(self):
     workDirectory = os.path.join(CCP4Utils.getTestTmpDir(),'test4')
     if not os.path.exists(workDirectory): os.mkdir(workDirectory)

     self.wrapper = pointless(parent=QTAPPLICATION(),name='test4',workDirectory=workDirectory)
     self.wrapper.container.loadDataFromXml(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','pointless','test_data','test4.data.xml'))

     self.wrapper.setWaitForFinished(1000000)
     pid = self.wrapper.process()
     self.wrapper.setWaitForFinished(-1)
     if len(self.wrapper.errorReport)>0: print(self.wrapper.errorReport.report())

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testpointless)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
