"Liz Potterton Aug 2010 - Base class to mimic QtCore.QObject"

#FIXME - SJM 2020 - Bleh, why...

class CObject():

    def __init__(self,parent=None,name=''):
        self.setParent(parent)
        self.setObjectName(name)
        if name is not None:
            if parent is not None:
                try:
                    self.__dict__['_objectPath'] = parent.objectPath() + '.' + name
                except:
                    self.__dict__['_objectPath'] = name
            else:
                    self.__dict__['_objectPath'] = name
        else:
            self.__dict__['_objectPath'] = ''
        self.__dict__['_signals'] = {}

    def parent(self):
        return self.__dict__['_parent']

    def setParent(self,parent):
        self.__dict__['_parent'] = parent

    def setObjectName(self,name=''):
        self.__dict__['_name'] = name

    def objectName(self):
        return self.__dict__['_name']

    def className(self):
        return str(self.__class__.__name__)
      
    def objectPath(self):
        return self.__dict__['_objectPath']

    def emitSignal(self,signal='finished',status=None):
        # BIG PROBLEM : what if handler object deleted - or worse
        # is a zombie because it has an instance here??
        if signal in self.__dict__['_signals']:
            for handler in self.__dict__['_signals'][signal]:
                try:
                    handler()
                except:
                    pass

    def connectSignal(self,origin=None,signal='',handler=None):
        return

    def disConnectSignal(self,origin=None,signal='',handler=None):
        pass
        
    def emit(self,**kw):
        pass

    def connect(self,**kw):
        pass

    def emitDataChanged(self):
        self.emitSignal(signal='dataChanged')

    def blockSignals(self):
        pass
