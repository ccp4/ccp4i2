
class testimport_mosflm(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     wrapper = import_mosflm(parent=QTAPPLICATION(),name='import_mosflm_test1')
     wrapper.container.loadDataFromXml()

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testimport_mosflm)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
