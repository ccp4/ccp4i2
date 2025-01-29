# Weird - this is failing with Qt error message - can not create widget when no GUI
# But works OK when do same thing manually 


class testDataManager(unittest.TestCase):

    def setUp(self):
        '''
        if GRAPHICAL():
          self.app = QTAPPLICATION()
          print 'testDataManager.setUp',GRAPHICAL(),type(self.app)
          self.dialog = QtWidgets.QDialog()
        '''
        self.testData = CFloat()
        self.manager = DATAMANAGER()

    def test1(self):
        if GRAPHICAL():
            app = QTAPPLICATION()
            widgetClass = self.manager.getWidgetClass(self.testData)
            self.assertEqual(widgetClass, CCP4Widgets.CFloatView, 'Failed to return correct widget class')
        else:
            self.fail('Can not test CCP4DataMananger in non-graphical mode')

    def test2(self):
        if GRAPHICAL():
            app = QTAPPLICATION()
            dialog = QtWidgets.QDialog()
            widget = self.manager.widget(model=self.testData,parentWidget=dialog)
            self.assertTrue(isinstance(widget,CCP4Widgets.CFloatView),'Failed to create CFloatView widget')
        else:
            self.fail('Can not test CCP4DataMananger in non-graphical mode')
