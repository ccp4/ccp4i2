
class testCmapcoeff(unittest.TestCase):

  def setUp(self):
    # make all background jobs wait for completion
    PROCESSMANAGER().setWaitForFinished(10000)

  def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)


  def testCmapcoeff(self):
    inputData =  CScriptDataContainer(name='cmapcoeff_test',containerType='inputData',initialise=cmapcoeff.INPUTDATA)
    outputData =  CScriptDataContainer(name='cmapcoeff_test',containerType='outputData',initialise=cmapcoeff.OUTPUTDATA)
    try:
      inputData.importXML(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','cmapcoeff','test_data','cmapcoeff_test_1.def.xml'))
    except CException as e:
      self.fail(e.errorType)
    try:
      outputData.importXML(os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','cmapcoeff','test_data','cmapcoeff_test_1.def.xml'))
    except CException as e:
      self.fail(e.errorType)

    wrapper = cmapcoeff()
    pid = wrapper.process()


def testSuite():
  suite = unittest.TestLoader().loadTestsFromTestCase(testCmapcoeff)
  return suite

def runAllTests():
  suite = testSuite()
  unittest.TextTestRunner(verbosity=2).run(suite)
