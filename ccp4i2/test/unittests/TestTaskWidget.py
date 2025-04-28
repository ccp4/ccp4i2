def TESTSUITE():
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(testTaskWidget)
    return suite

def runAllTests():
    suite = TESTSUITE()
    unittest.TextTestRunner(verbosity=2).run(suite)

#-------------------------------------------------------------------
class CSummatTask(CFolderTaskWidget):
#-------------------------------------------------------------------
# Subclass CTaskWidget to give specific task window

    def __init__(self, parent):
        CTaskWidget.__init__(self, parent=None)

    def drawContents(self):
        #self.setProgramHelpFile('fft')
        #self.openFolder(folderFunction='protocol')
        self.openFolder(title='Test folder')
        #self.createTitleLine()
        self.createLine(['widget', 'nCycles', 'label', 'cycles with cutoff', 'widget', 'cutoff'])
        self.createLine(['label', 'Range of gubbins', 'widget', 'range' ] )


#class testTaskWidget(unittest.TestCase):
class testTaskWidget():

    def __init__(self):
        self.setUp()

    def setUp(self):
        self.app = QTAPPLICATION()
        self.window = QtWidgets.QMainWindow()

    def test1(self):
        summat = CContainer(parent=self.app, definitionFile='/Users/lizp/Desktop/dev/ccp4i2/sandpit/summat.contents.xml')
        self.assertEqual(len(summat.CONTENTS), 4, 'Container - Wrong content length')
        self.assertEqual(summat.gubbins.startValue, 12.0, 'Container - sub-container CFloat wrong initial value')

    def test2(self):
        self.container = CContainer(parent=self.app, definitionFile='/Users/lizp/Desktop/dev/ccp4i2/sandpit/summat.contents.xml')
        self.task = CSummatTask(self.window)
        self.task.setContainer(self.container)
        self.task.draw()
        self.window.setCentralWidget(self.task)
        self.window.show()
        sys.exit(self.app.exec_())

    def tearDown(self):
        self.app.quit()
