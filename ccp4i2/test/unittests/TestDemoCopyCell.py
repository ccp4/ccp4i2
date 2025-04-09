
class test_demo_copycell(unittest.TestCase):
  
  def setUp(self):
    # make all background jobs wait for completion
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

  def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

  def test_1(self):
    # Run the pipeline
    wrapper = demo_copycell(parent=QTAPPLICATION(),name='demo_copycell')
    wrapper.container.loadDataFromXml(os.path.join(getCCP4I2Dir(),'pipelines','demo_copycell','test_data','test_1.params.xml'))
    # Ensure no output file exists already
    xyzout = wrapper.container.outputData.XYZOUT.fullPath.get()
    if xyzout is not None and os.path.exists(xyzout): os.remove(xyzout)
    pid = wrapper.process()

    #test if output file created
    self.assertEqual(os.path.exists( xyzout),1,'Failed to create copied pdb file '+xyzout)                             


def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(test_demo_copycell)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
