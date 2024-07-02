"""
     core/CCP4Peferences.py: CCP4 Gui Project
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

'''
    Liz Potterton - create class to maintain GUIPreferences - Sept 2011
'''
import os
import shutil
from PySide6 import QtCore
from core import CCP4Container

class CPreferences(CCP4Container.CContainer):

    preferencesSaved = QtCore.Signal()

    insts = None

    def __init__(self):
        from core import CCP4Modules
        CCP4Container.CContainer.__init__(self)
        CPreferences.insts = self
        defFile = CCP4Modules.TASKMANAGER().searchDefFile('guipreferences')
        self.loadContentsFromXml(defFile)
        self.load()

    def load(self):
        prefFile = self.preferencesFile()
        if os.path.exists(prefFile):
            self.loadDataFromXml(prefFile)

    def preferencesFile(self):
        from core import CCP4Utils
        return  os.path.join(CCP4Utils.getDotDirectory(), 'configs', 'guipreferences.params.xml')

    def save(self):
        prefFile = self.preferencesFile()
        if os.path.exists(prefFile):
            shutil.copyfile(prefFile, prefFile + '.bak')
        self.saveDataToXml(prefFile)
        self.preferencesSaved.emit()

    def guiUpdateEnabled(self):
        return bool(self.get('BZR_DOWNLOAD'))
