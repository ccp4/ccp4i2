from PySide2 import QtCore

class CObject(QtCore.QObject):

    dataChanged = QtCore.Signal()
    finished = QtCore.Signal(int)

    def __init__(self, parent=None, name=None):
        QtCore.QObject.__init__(self, parent)
        if name is not None:
            self.setObjectName(name)

    @QtCore.Slot()
    def emitDataChanged(self):
        self.dataChanged.emit()

    def objectName(self):
        name = QtCore.QObject.objectName(self)
        if name is None:
            return None
        else:
            return str(name)
    
    def className(self):
        return str(self.__class__.__name__)

    def connectSignal(self,origin=None,signal='',handler=None):
        if signal == "finished":
            origin.finished.connect(handler)
