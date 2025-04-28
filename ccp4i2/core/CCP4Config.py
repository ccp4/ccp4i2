from pathlib import Path
import xml.etree.ElementTree as ET

import gemmi

from .. import __version__
from .CCP4Utils import getDotDirectory


gemmi.set_leak_warnings(False)


def CONFIG(filePath: str = None, developer: bool = None, graphical: bool = None):
    if CConfig.insts is None:
        CConfig.insts = CConfig(filePath, developer, graphical)
    return CConfig.insts


class CConfig:
    insts = None

    def __init__(self, filePath, developer, graphical):
        print("ccp4i2 version", __version__)
        self.dbFile = None
        self.dbUser = None
        self.developer = True
        self.graphical = False
        self.maxRunningProcesses = 4

        if filePath is None:
            dotDir = Path(getDotDirectory()).resolve()
            filePath = dotDir / "configs" / "ccp4i2_config.params.xml"
        else:
            filePath = Path(filePath).resolve()
        if filePath.exists():
            self._loadDataFromXml(filePath)
        else:
            self._saveDataToXml(filePath)
        if developer is not None:
            self.developer = developer
        if graphical is not None:
            self.graphical = graphical

    def _loadDataFromXml(self, filePath: Path):
        root = ET.parse(filePath).getroot()
        if (tag := root.find(".//graphical")):
            self.graphical = tag.text.lower() == "true"
        if (tag := root.find(".//developer")):
            self.developer = tag.text.lower() == "true"
        if (tag := root.find(".//dbFile")):
            self.dbFile = None if tag.text.lower() == "none" else tag.text
        if (tag := root.find(".//dbUser")):
            self.dbUser = None if tag.text.lower() == "none" else tag.text
        if (tag := root.find(".//maxRunningProcesses")):
            self.maxRunningProcesses = 1 if tag.text.lower() == "none" else int(tag.text)

    def _saveDataToXml(self, filePath: Path):
        configs = ET.Element("configs")
        for tag in [
            "dbFile",
            "dbUser",
            "developer",
            "graphical",
            "maxRunningProcesses",
        ]:
            config = ET.Element(tag)
            config.text = str(getattr(self, tag))
            configs.append(config)
        tree = ET.ElementTree(configs)
        tree.write(filePath)
