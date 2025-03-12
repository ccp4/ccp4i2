import sys

from PySide2 import QtCore, QtWidgets

from ..core.CCP4Config import CONFIG


_QTAPPLICATION = None


def QTAPPLICATION(graphical=None):
    global _QTAPPLICATION
    if _QTAPPLICATION is None:
        if graphical is None:
            graphical = CONFIG().graphical
        if graphical:
            _QTAPPLICATION = CGuiApplication(sys.argv)
        else:
            _QTAPPLICATION = CApplication(sys.argv)
    return _QTAPPLICATION


class CApplication(QtCore.QCoreApplication):

    doCheckForFinishedJobs = QtCore.Signal()

    def __init__(self, args):
        super().__init__(self, args)


class CGuiApplication(QtWidgets.QApplication):

    doCheckForFinishedJobs = QtCore.Signal()
    prefsChanged = QtCore.Signal()

    def __init__(self, args):
        super().__init__(self, args)

    def objectPath(self):
        return ''

    def setNamedStyle(self, styleName):
        factory = QtWidgets.QStyleFactory()
        self.setStyle(factory.create(str(styleName)))

    def getStyleKeys(self):
        factory = QtWidgets.QStyleFactory()
        return map(str, factory.keys())

    def screenSize(self):
        rect = self.desktop().screenGeometry()
        return rect.width(), rect.height()
