from os import environ
from pathlib import Path
from core.CCP4PluginScript import CPluginScript
from core.CCP4ModelData import CPdbDataFile


class coot1(CPluginScript):
    TASKNAME = "coot1"
    TASKMODULE = "model_building"
    TASKTITLE = "Coot 1"
    TASKCOMMAND = str(Path(environ["CCP4"], "coot_py3", "bin", "coot"))
    TASKVERSION = 0.1
    DESCRIPTION = "Interactive model building with Coot 1)"
    ASYNCHRONOUS = True
    TIMEOUT_PERIOD = 9999999.9
    MAINTAINER = "paul.bond@york.ac.uk"
    WHATNEXT = ["prosmart_refmac", "coot_rebuild", "modelcraft"]

    ERROR_CODES = {}

    def makeCommandAndScript(self, container=None):
        inputData = self.container.inputData

        if inputData.XYZIN_LIST.isSet():
            for path in inputData.XYZIN_LIST:
                self.appendCommandLine(["--coords", path])

        if inputData.DICT.isSet():
            self.appendCommandLine(["--dictionary", inputData.DICT])

        script = ["import coot"]
        if inputData.FPHIIN_LIST.isSet():
            for path in inputData.FPHIIN_LIST:
                script.append(f"coot.read_mtz('{path}', 'F', 'PHI', '', False, False)")
        if inputData.DELFPHIIN_LIST.isSet():
            for path in inputData.DELFPHIIN_LIST:
                script.append(f"coot.read_mtz('{path}', 'F', 'PHI', '', False, True)")
        if inputData.DELFPHIINANOM_LIST.isSet():
            for path in inputData.DELFPHIINANOM_LIST:
                script.append(
                    f"imap = coot.read_mtz('{path}', 'F', 'PHI', '', False, True)"
                )
                script.append("coot.set_map_colour(imap, 0.75, 0.9, 0.75)")
        scriptPath = Path(self.getWorkDirectory(), "coot_script.py")
        with scriptPath.open("w", encoding="utf-8") as stream:
            stream.write("\n".join(script))
        self.appendCommandLine(["--script", scriptPath])

        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        xyzout = self.container.outputData.XYZOUT
        workDir = Path(self.getWorkDirectory())
        index = 0
        for pattern, contentFlag in [
            ("*.cif", CPdbDataFile.CONTENT_FLAG_MMCIF),
            ("*.pdb", CPdbDataFile.CONTENT_FLAG_PDB),
        ]:
            for path in workDir.glob(pattern):
                xyzout[index].setFullPath(str(path))
                xyzout[index].annotation.set(f"Coot output: {path.name}")
                xyzout[index].subType.set(CPdbDataFile.SUBTYPE_MODEL)
                xyzout[index].contentFlag.set(contentFlag)
                index += 1
        self.container.outputData.XYZOUT.set(xyzout[:index])
        return CPluginScript.SUCCEEDED
