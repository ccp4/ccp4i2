"""
Copyright (C) 2011 University of York
Liz Potterton - create class to maintain GUIPreferences - Sept 2011
"""

import os
import shutil

from PySide2 import QtCore

from . import CCP4Container
from . import CCP4Utils


def PREFERENCES():
    if CPreferences.insts is None:
        CPreferences()
    return CPreferences.insts


class CPreferences(CCP4Container.CContainer):

    preferencesSaved = QtCore.Signal()

    insts = None

    def __init__(self):
        from . import CCP4TaskManager
        super().__init__()
        CPreferences.insts = self
        defFile = CCP4TaskManager.TASKMANAGER().searchDefFile('guipreferences')
        self.loadContentsFromXml(defFile)
        self.load()

    def load(self):
        prefFile = self._preferencesFile()
        if os.path.exists(prefFile):
            self.loadDataFromXml(prefFile)

    def save(self):
        prefFile = self._preferencesFile()
        if os.path.exists(prefFile):
            shutil.copyfile(prefFile, prefFile + '.bak')
        self.saveDataToXml(prefFile)
        self.preferencesSaved.emit()

    def _preferencesFile(self):
        return  os.path.join(CCP4Utils.getDotDirectory(), 'configs', 'guipreferences.params.xml')
