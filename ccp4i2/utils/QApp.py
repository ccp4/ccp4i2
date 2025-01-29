import sys

from PySide2 import QtCore, QtWidgets

from ..core.CCP4Config import GRAPHICAL


MYAPPLICATION = None


def QTAPPLICATION(graphical=None):
    global MYAPPLICATION
    # NB can not use QApplication.instance() as it returns the QApplication object
    # rather than CApplication. Suppose could reimplement CApplication.instance() but
    # would just do the same as this function.
    if graphical is None:
        graphical = GRAPHICAL()
    if MYAPPLICATION is None:
        #print 'Starting Qt, graphical mode:',graphical
        if graphical:
            MYAPPLICATION = CGuiApplication(sys.argv)
        else:
            MYAPPLICATION = CApplication(sys.argv)
    return MYAPPLICATION


class CApplication(QtCore.QCoreApplication):

    doCheckForFinishedJobs = QtCore.Signal()

    def __init__(self,args):
        QtCore.QCoreApplication.__init__(self,args)

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
