def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testLaunch)
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testLaunch(unittest.TestCase):

    def test1(self):
        l = CLauncher()
        com = l.makeCommand(viewer='CCP4mg',command='openFile',data='foo/bar')
        print('testLaunch.test1',com)
        self.assertEqual(com[0:7],'<begin>','makeCommand failed for openFile')
