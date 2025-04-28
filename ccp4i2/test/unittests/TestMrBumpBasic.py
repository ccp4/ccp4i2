
class testmrbump_basic( unittest.TestCase ) :

#- def setUp( self ) :
#- def tearDown( self ) :
#- def test2( self ) :

   def test1( self ) :
      xmlInput = os.path.join( getCCP4I2Dir(), 'wrappers', 'mrbump_basic', 'test_data', 'test1'+'.params.xml' )
      self.wrapper = mrbump_basic( name='job' )
      self.wrapper.container.loadDataFromXml( xmlInput )
      self.wrapper.setWaitForFinished( 1000000 )

      pid = self.wrapper.process()
      self.wrapper.setWaitForFinished( -1 )
      if len(self.wrapper.errorReport)>0:
         print(self.wrapper.errorReport.report())

def TESTSUITE() :
   suite = unittest.TestLoader().loadTestsFromTestCase( testmrbump_basic )
   return suite

def testModule() :
   suite = TESTSUITE()
   unittest.TextTestRunner( verbosity=2 ).run( suite )
