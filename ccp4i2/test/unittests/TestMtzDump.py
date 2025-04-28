
class testMtzdump(unittest.TestCase):
  
  def setUp(self):
    # make all background jobs wait for completion
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

  def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

  def testMtzdump(self):
    self.wrapper = mtzdump(parent=QTAPPLICATION(),name='test_mtzdump')
    self.wrapper.container.inputData.HKLIN.set({'project':'CCP4I2_TOP','baseName':'gere_nat.mtz','relPath':'test/data'})
    pid = self.wrapper.process()
    print(self.wrapper.container.outputData.CELL)
    if len(self.wrapper.errorReport)>0: self.wrapper.errorReport.report()
    self.assertEqual(self.wrapper.container.outputData.CELL.a,108.742,'Mtzdump output CELL wrong')


def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testMtzdump)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
