from __future__ import print_function

from ccp4i2.baselayer import QtCore,QtGui, QtWidgets
import os

MYAPPLICATION = None

class CApplication(QtCore.QCoreApplication):

    doCheckForFinishedJobs = QtCore.Signal()

    def __init__(self,args):
        QtCore.QCoreApplication.__init__(self,args)
        """
        if os.path.exists(os.path.join(os.environ["CCP4"],"lib","qt4","plugins")):
            self.addLibraryPath(os.path.join(os.environ["CCP4"],"lib","qt4","plugins"))
        elif os.path.exists(os.path.join(os.environ["CCP4"],"plugins")):
            self.addLibraryPath(os.path.join(os.environ["CCP4"],"plugins"))
        elif os.path.exists(os.path.join(os.environ["CCP4"],"QtPlugins")):
            self.addLibraryPath(os.path.join(os.environ["CCP4"],"QtPlugins"))
        elif os.path.exists(os.path.join(os.environ["CCP4MG"],"QtPlugins")):
            self.addLibraryPath(os.path.join(os.environ["CCP4MG"],"QtPlugins"))
        """

        def mainWindow(self):
            return None
      
        def objectPath(self):
            return ''


class CGuiApplication(QtWidgets.QApplication):

    doCheckForFinishedJobs = QtCore.Signal()

    prefsChanged = QtCore.Signal()
    def __init__(self,args):
        #print 'CGuiApplication.__init__',self
        QtWidgets.QApplication.__init__(self,args)
        """
        if os.path.exists(os.path.join(os.environ["CCP4"],"lib","qt4","plugins")):
            self.addLibraryPath(os.path.join(os.environ["CCP4"],"lib","qt4","plugins"))
        elif os.path.exists(os.path.join(os.environ["CCP4"],"plugins")):
            self.addLibraryPath(os.path.join(os.environ["CCP4"],"plugins"))
        elif os.path.exists(os.path.join(os.environ["CCP4"],"QtPlugins")):
            self.addLibraryPath(os.path.join(os.environ["CCP4"],"QtPlugins"))
        elif os.path.exists(os.path.join(os.environ["CCP4MG"],"QtPlugins")):
            self.addLibraryPath(os.path.join(os.environ["CCP4MG"],"QtPlugins"))
        """
        #import CCP4StyleSheet
        #CCP4StyleSheet.setStyleSheet(self)

    def objectPath(self):
        return ''

    def setNamedStyle(self,styleName):
        print('CGuiApplication.setStyle',styleName)
        factory = QtWidgets.QStyleFactory()
        self.setStyle(factory.create(str(styleName)))

    def getStyleKeys(self):
        factory = QtWidgets.QStyleFactory()
        keyList = list(factory.keys())
        print('getStyleKeys',type(keyList),keyList)
        pyKeyList = []
        for item in keyList: pyKeyList.append(str(item))
        return pyKeyList

    def screenSize(self):
        rect = self.desktop().screenGeometry()
        return rect.width(),rect.height()
