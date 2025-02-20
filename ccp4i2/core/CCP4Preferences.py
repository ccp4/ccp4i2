"Class to maintain GUIPreferences"

import os
import shutil

from PySide2 import QtCore

from . import CCP4Container
from . import CCP4TaskManager
from . import CCP4Utils


def PREFERENCES():
    if CPreferences.insts is None:
        CPreferences()
    return CPreferences.insts


class CPreferences(CCP4Container.CContainer):

    preferencesSaved = QtCore.Signal()

    insts = None

    def __init__(self):
        CCP4Container.CContainer.__init__(self)
        CPreferences.insts = self
        defFile = CCP4TaskManager.TASKMANAGER().searchDefFile('guipreferences')
        self.loadContentsFromXml(defFile)
        self.load()

    def load(self):
        prefFile = self.preferencesFile()
        if os.path.exists(prefFile):
            self.loadDataFromXml(prefFile)

    def preferencesFile(self):
        return  os.path.join(CCP4Utils.getDotDirectory(), 'configs', 'guipreferences.params.xml')

    def save(self):
        prefFile = self.preferencesFile()
        if os.path.exists(prefFile):
            shutil.copyfile(prefFile, prefFile + '.bak')
        self.saveDataToXml(prefFile)
        self.preferencesSaved.emit()
