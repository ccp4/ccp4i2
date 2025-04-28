
class testController(unittest.TestCase):

    def setUp(self):
        self.controller = CJobController()

    def testRunFreerflag(self):
        controlFile = os.path.join(CCP4Utils.getCCP4I2Dir(),'wrappers','freerflag','test_data','test1.data.xml')
        stdoutFile = os.path.join(CCP4Utils.getTestTmpDir(),'CJobController_test.stdout')
        c = CCP4Container.CContainer()
        c.loadDataFromXml(controlFile)
        output = str(c.outputData.HKLOUT)
        print('Output file:',output,'Stdout:',stdoutFile)
        if os.path.exists(output): os.remove(output)
        self.controller.runTask(fileName=controlFile)
        self.assertTrue(os.path.exists(output),'No output file created')


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(testController)
    return suite
    
def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)
