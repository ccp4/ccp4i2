from __future__ import print_function

"""
     CCP4PluginThread.py: CCP4 GUI Project
     Copyright (C) 2013 STFC

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

'''
   Liz Potterton -Jan 2013 - A thread to run plugins
'''
    
from core import CCP4PluginScript
from PySide6 import QtCore

class CRunPluginThread(QtCore.QThread):

    finished = QtCore.Signal()

    def __init__(self,parent=None,ccp4i2Path=None,comFilePath=None):
        QtCore.QThread.__init__(self,parent)    
        self.runPlugin = CCP4PluginScript.CRunPlugin(parent=self,ccp4i2Path=ccp4i2Path,comFilePath=comFilePath)
        self.runPlugin.finished.connect(self.handleFinishSignal)

    def run(self):   
        self.runPlugin.run()
        self.exec_()

    @QtCore.Slot()
    def handleFinishSignal(self):
        self.finished.emit()

    def checkIfFinished(self):
        from core import CCP4Modules
        if len(CCP4Modules.PROCESSMANAGER().activeProcesses())>0:
            print('checkIfFinished',len(CCP4Modules.PROCESSMANAGER().activeProcesses()))
        else:
            self.finished.emit()
