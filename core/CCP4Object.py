from __future__ import print_function

"""
     CCP4Object.py: CCP4 GUI Project
     Copyright (C) 2010 University of York

     This library is free software: you can redistribute it and/or
     modify it under the terms of the GNU Lesser General Public License
     version 3, modified in accordance with the provisions of the 
     license to address the requirements of UK law.
 
     You should have received a copy of the modified GNU Lesser General 
     Public License along with this library.  If not, copies may be 
     downloaded from http://www.ccp4.ac.uk/ccp4license.php
 
     This program is distributed in the hope that it will be useful,
     but WITHOUT ANY WARRANTY; without even the implied warranty of
     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
     GNU Lesser General Public License for more details.
"""

"""
   Liz Potterton Aug 2010 - Base class to mimic QtCore.QObject
"""

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
        #print 'objectPath',self.objectName(),self.__dict__['_objectPath']
        return self.__dict__['_objectPath']  

    '''
    def getName(self):
      if self._parent is None:
        return None
      else:
        for key,obj in self._parent._value.items():
          if obj is not None and obj == self:
            return key
        return None
    '''

    def emitSignal(self,signal='finished',status=None):
        # BIG PROBLEM : what if handler object deleted - or worse
        # is a zombie because it has an instance here??
        '''
        if status is not None:
            self.emit(QtCore.SIGNAL(signal),status)
        else:
            self.emit(QtCore.SIGNAL(signal))
        '''
        if signal in self.__dict__['_signals']:
            for handler in self.__dict__['_signals'][signal]:
                try:
                    handler()
                except:
                    pass

    def connectSignal(self,origin=None,signal='',handler=None):
        '''
        self.connect(origin,QtCore.SIGNAL(signal),handler)
        '''
        return
        if origin != self:
            try:
                origin.connectSignal(origin=origin,signal=signal,handler=handler)
            except:
                pass
        if signal not in self._signals:
            self._signals[signal] = []
        if not self._signals[signal].count(handler): 
            self._signals[signal].append(handler)

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
