
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
   Liz Potterton Aug 2010 -Wrapper fQtCore.or QtCore.QObject complementing CCP4Object
"""

from PySide6 import QtCore

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

    '''
    def getName(self):
      if self.parent() is None:
        return None
      else:
        for key,obj in self.parent()._value.items():
          print 'CObject.getName',key,obj
          if obj is not None and obj == self:
            return key
        return None
    '''

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
