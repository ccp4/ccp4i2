
class test_coot_fit_residues(unittest.TestCase):
    
    def setUp(self):
        # make all background jobs wait for completion
        # this is essential for unittest to work
        self.app = QTAPPLICATION()
        PROCESSMANAGER().setWaitForFinished(10000)
    
    def tearDown(self):
        PROCESSMANAGER().setWaitForFinished(-1)
    
    def test_1(self):
        wrapper = coot_fit_residues(parent=QTAPPLICATION(),name='coot_fit_residues_test1')
        wrapper.container.loadDataFromXml()


def TESTSUITE():
    suite = unittest.TestLoader().loadTestsFromTestCase(test_coot_fit_residues)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

