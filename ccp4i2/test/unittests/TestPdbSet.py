
class testPdbset(unittest.TestCase):

   def setUp(self):
    # make all background jobs wait for completion
    self.app = QTAPPLICATION()
    PROCESSMANAGER().setWaitForFinished(10000)

   def tearDown(self):
    PROCESSMANAGER().setWaitForFinished(-1)

   def testPdbset(self):
     wrapper = pdbset(parent=QTAPPLICATION(),name='pdbset')
     
     wrapper.container.inputData.XYZIN.set(project='CCP4I2_TOP',relPath='wrappers/pdbset/test_data',baseName='1df7.pdb')
     wrapper.container.inputData.CELL.set(a=100.0,b=120.0,c=30.0,alpha=90.0,beta=90.0,gamma=89.1)
     wrapper.container.outputData.XYZOUT.set(project='CCP4I2_TEST',baseName='1df7_mangled.pdb')
     print('XYZOUT',wrapper.container.outputData.XYZOUT.fullPath)

     # Ensure no output file exists
     if wrapper.container.outputData.XYZOUT.exists():  os.remove(wrapper.container.outputData.XYZOUT.fullPath.get())
     wrapper.process()

     #test if output file created
     self.assertEqual(wrapper.container.outputData.XYZOUT.exists(),True,'Failed to create copied pdb file')                             

def TESTSUITE():
  suite = unittest.TestLoader().loadTestsFromTestCase(testPdbset)
  return suite

def testModule():
  suite = TESTSUITE()
  unittest.TextTestRunner(verbosity=2).run(suite)
