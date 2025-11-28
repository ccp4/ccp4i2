from functools import reduce
from math import gcd
from os import environ, rename
from pathlib import Path
import gemmi
from core.CCP4ErrorHandling import SEVERITY_WARNING
from core.CCP4PluginScript import CPluginScript

class phaser_tng_picard(CPluginScript):
    TASKNAME = "phaser_tng_picard"
    TASKMODULE = "molecular_replacement"
    TASKTITLE = "Molecular replacement with Phaser TNG Picard"
    TASKCOMMAND = str(Path(environ["CCP4"], "libexec", "phasertng.picard"))
    TASKVERSION = 0.1
    MAINTAINER = "paul.bond@york.ac.uk"
    WHATNEXT = ["prosmart_refmac", "coot_rebuild", "modelcraft"]

    ERROR_CODES = {}

    def processInputFiles(self):
        ifdir = Path(self.getWorkDirectory(), "input-files")
        ifdir.mkdir(exist_ok=True)
        miniMtzs = ["OBSIN"]
        if self.container.inputData.FREERFLAG.isSet():
            miniMtzs.append("FREERFLAG")
        hklin, columns, error = self.makeHklin0(miniMtzs)
        if error.maxSeverity() > SEVERITY_WARNING:
            return CPluginScript.FAILED
        rename(hklin, ifdir / "hklin.mtz")
        asu = self.container.inputData.ASUIN
        copies = [int(s.nCopies) for s in asu.fileContent.seqList]
        divisor = reduce(gcd, copies)
        with (ifdir / "seqin.fasta").open("w", encoding="utf-8") as f:
            for seq in asu.fileContent.seqList:
                for i in range(int(seq.nCopies) // divisor):
                    f.write(f">*{seq.polymerType}* {seq.name} - Copy {i + 1}\n")
                    f.write(f"{str(seq.sequence)}\n")
        for i, path in enumerate(self.container.inputData.XYZIN_LIST):
            st = gemmi.read_structure(str(path), format=gemmi.CoorFormat.Detect)
            st.write_pdb(str(ifdir / f"model_{i}.pdb"))
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, **kw):
        self.appendCommandLine("directory=input-files")
        self.appendCommandLine("software=refmac")
        return CPluginScript.SUCCEEDED
