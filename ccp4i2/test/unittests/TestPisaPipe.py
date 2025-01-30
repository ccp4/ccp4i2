# PLUGIN TESTS
# See Python documentation on unittest module

class testpisa(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    # this is essential for unittest to work
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def test_1(self):
     wrapper = pisapipe(parent=QTAPPLICATION(),name='pisa_test1')
     wrapper.container.loadDataFromXml()
     

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testpisa)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
