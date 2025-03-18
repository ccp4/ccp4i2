from pathlib import Path
import os
import shutil

from .. import I2_TOP
from ..core import CCP4Container
from ..core import CCP4Utils


def SERVERSETUP():
    if CServerSetup.insts is None:
        CServerSetup()
    return CServerSetup.insts


class CServerSetup(CCP4Container.CContainer):
    insts = None

    def __init__(self):
        super().__init__(name="SERVER_SETUP")
        CServerSetup.insts = self
        self.__dict__["source"] = None
        from ..core.CCP4TaskManager import TASKMANAGER
        defFile = TASKMANAGER().searchDefFile("serverSetup")
        self.loadContentsFromXml(defFile)

    def load(self, source=None):
        prefFile, source = self.preferencesFile(source)
        if prefFile.exists():
            # Beware if params file has >1 serverGroup we need to add the 2+ groups to the
            # the container as they are not in the def file
            from . import CCP4File, CCP4Annotation
            fObj = CCP4File.CI2XmlDataFile(prefFile)
            for sGEle in fObj.getBodyEtree():
                if self.get(str(sGEle.tag)) is None:
                    self.setContents(
                        {str(sGEle.tag): {"class": CCP4Annotation.CServerGroup}}
                    )
            self.loadDataFromXml(prefFile)
            self.__dict__["source"] = source

    def writeAccess(self, source):
        parent = self.preferencesFile(source)[0].parent
        while not parent.exists():
            parent = parent.parent
        return os.access(parent, os.W_OK | os.X_OK)

    def preferencesFile(self, source: str = None):
        if source is None or source == "user":
            dotDir = Path(CCP4Utils.getDotDirectory()).resolve()
            filename = dotDir / "configs" / "serverSetup.params.xml"
            if filename.exists() or source == "user":
                return filename, "user"
        filename = I2_TOP / "local_setup" / "serverSetup.params.xml"
        return filename, "installation"

    def save(self, source=None):
        prefFile, source = self.preferencesFile(source)
        if prefFile.exists():
            shutil.copyfile(prefFile, str(prefFile) + ".bak")
        prefFile.parent.mkdir(parents=True, exist_ok=True)
        self.saveDataToXml(fileName=str(prefFile))
        self.__dict__["source"] = source
