
class testcsymmatch( unittest.TestCase ) :

   def setUp(self):
    self.app = QTAPPLICATION()
    # make all background jobs wait for completion
    # this is essential for unittest to work
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def test1( self ) :
      workDirectory = CCP4Utils.getTestTmpDir()
      xmlInput = os.path.join( CCP4Utils.getCCP4I2Dir(), 'wrappers', 'csymmatch', 'test_data', 'test1'+'.params.xml' )
      self.wrapper = csymmatch(parent=QTAPPLICATION(), name='csymmatch_test1',workDirectory=workDirectory)
      self.wrapper.container.loadDataFromXml( xmlInput )

      self.wrapper.setWaitForFinished( 1000000 )
      pid = self.wrapper.process()
      self.wrapper.setWaitForFinished( -1 )
      if len(self.wrapper.errorReport)>0:
         print(self.wrapper.errorReport.report())

def TESTSUITE() :

   suite = unittest.TestLoader().loadTestsFromTestCase( testcsymmatch )
   return suite

def testModule() :

   suite = TESTSUITE()
   unittest.TextTestRunner( verbosity=2 ).run( suite )

