class test_arcimboldo ( unittest.TestCase ) :
    def setUp(self):
      # make all background jobs wait for completion
      PROCESSMANAGER().setWaitForFinished(10000)

    def tearDown(self):
      PROCESSMANAGER().setWaitForFinished(-1)

    def test_arcimboldo(self):
      inputData =  CScriptDataContainer(name='test_arcimboldo_test',containerType='inputData',initialise=test_arcimboldo.INPUTDATA)
      outputData =  CScriptDataContainer(name='test_arcimboldo_test',containerType='outputData',initialise=test_arcimboldo.OUTPUTDATA)
      try:
         inputData.importXML(os.path.join(getCCP4I2Dir(),'wrappers','test_arcimboldo','test_data','test_arcimboldo_test_1.def.xml'))
      except CException as e:
         self.fail(e.errorType)
      try:
         outputData.importXML(os.path.join(getCCP4I2Dir(),'wrappers','test_arcimboldo','test_data','test_arcimboldo_test_1.def.xml'))
      except CException as e:
         self.fail(e.errorType)
      
      wrapper = test_arcimboldo()
      pid = wrapper.process()

    def testSuite():
      suite = unittest.TestLoader().loadTestsFromTestCase(test_test_arcimboldo)
      return suite

    def runAllTests():
      suite = testSuite()
      unittest.TextTestRunner(verbosity=2).run(suite)
