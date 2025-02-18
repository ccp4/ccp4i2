
class test_demo_multi_mtzdump(unittest.TestCase):

  def setUp(self):
    # make all background jobs wait for completion
    if not CONFIG().graphical:
      PROCESSMANAGER().setWaitForFinished(10000)

  def tearDown(self):
    if not CONFIG().graphical:
      PROCESSMANAGER().setWaitForFinished(-1)

  def test_1(self):
    # Run the pipeline
    pipe = demo_multi_mtzdump(parent=QTAPPLICATION(),name='demo_multi_mtzdump')
    pipe.process()
      

def testSuite():
  suite = unittest.TestLoader().loadTestsFromTestCase(test_demo_multi_mtzdump)
  return suite

def runAllTests():
  suite = testSuite()
  unittest.TextTestRunner(verbosity=2).run(suite)
