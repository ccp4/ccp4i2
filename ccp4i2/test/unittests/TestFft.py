 
class testFft(unittest.TestCase):
  
  def setUp(self):
    # make all background jobs wait for completion
    PROCESSMANAGER().setWaitForFinished(10000)

  def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)


  def testFft(self):
    inputData =  CScriptDataContainer(name='fft_test',containerType='inputData',initialise=fft.INPUTDATA)
    outputData =  CScriptDataContainer(name='fft_test',containerType='outputData',initialise=fft.OUTPUTDATA)
    try:
      inputData.importXML(os.path.join(getCCP4I2Dir(),'wrappers','fft','test_data','fft_test_1.def.xml'))
    except CException as e:
      self.fail(e.errorType)
    try:
      outputData.importXML(os.path.join(getCCP4I2Dir(),'wrappers','fft','test_data','fft_test_1.def.xml'))
    except CException as e:
      self.fail(e.errorType)
      
    wrapper = fft()
    pid = wrapper.process()


def testSuite():
  suite = unittest.TestLoader().loadTestsFromTestCase(testFft)
  return suite

def runAllTests():
  suite = testSuite()
  unittest.TextTestRunner(verbosity=2).run(suite)
