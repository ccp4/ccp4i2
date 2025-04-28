def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testBuildLookups)
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(testLoadLookups))
    return suite

def testModule():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

class testBuildLookups(unittest.TestCase):
    def test1(self):
        taskManager = CTaskManager()
        taskManager.buildLookupFromScratch()
        self.assertTrue(len(taskManager.taskLookup) >= 97,msg="Built less than 97 tasks...possibly something fishy")

class testLoadLookups(unittest.TestCase):
    def test1(self):
        taskManager = CTaskManager()
        taskManager.loadCachedClassLookups()
        self.assertTrue(len(taskManager.taskLookup) >= 97,msg="Loaded less than 97 tasks...possibly something fishy")
