from __future__ import print_function

"""
     CCP4ProgramLog.py: CCP4 GUI Project
     Copyright (C) 2011 University of York

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
   Liz Potterton Aug 2011 - Separate ProgramLog from ErrorHandling
"""
import sys
import time
import traceback

from core.CCP4Config import QT
if QT():
    from core.CCP4QtObject import CObject
else:
    from core.CCP4Object import CObject
from core.CCP4ErrorHandling import *


class CProgramLog(CErrorReport, CObject):
    insts = None
    def __init__(self, parent=None):
        CErrorReport.__init__(self)
        CObject.__init__(self, parent=None)
      
    def append(self, cls=None, code=0, details=None, name=None):
        CErrorReport.append(self, cls=cls, code=code, details=details, name=name, recordTime=True)
        self.dataChanged.emit()
    
    def extend(self, other=None, label=None):
        CErrorReport.extend(self, other=other, label=label, recordTime=True)
        self.dataChanged.emit()

    def getTimeRange(self):
        t = []
        if len(self._reports) > 0:
            t.append(self._reports[0]['time'])
            t.append(self._reports[-1]['time'])
        else:
            t.append(time.localtime())
            t.append(time.localtime())
        return t


def showTrace(messageBox=False, title="Crash report"):
    exc_type, exc_value, exc_tb = sys.exc_info()[:3]
    sys.stderr.write(str(exc_type) + '\n')
    sys.stderr.write(str(exc_value) + '\n')
    err = ''
    for s in traceback.extract_stack()[:-1]:
        err = err + '  ' + str(s[0]) + ',  line ' + str(s[1]) + '\n'
        err = err + '    ' + str(s[3]) + '\n'
    print(err)
    if messageBox and GRAPHICAL():
        from PySide6 import QtGui, QtWidgets
        QtWidgets.QMessageBox.critical(self, title, title + ":\n"+str(exc_type) + '\n' + str(exc_value) + '\n' + err)

