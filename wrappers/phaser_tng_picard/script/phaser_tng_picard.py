from functools import reduce
from math import gcd
from os import environ
from pathlib import Path
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
        self.hklin, self.columns, error = self.makeHklin0(["OBSIN", "FREERFLAG"])
        if error.maxSeverity() > SEVERITY_WARNING:
            return CPluginScript.FAILED
        self.seqin = str(Path(self.getWorkDirectory(), "seqin.fasta"))
        asu = self.container.inputData.ASUIN
        copies = [int(s.nCopies) for s in asu.fileContent.seqList]
        divisor = reduce(gcd, copies)
        with open(self.seqin, "w", encoding="utf-8") as f:
            for seq in asu.fileContent.seqList:
                for i in range(int(seq.nCopies) // divisor):
                    f.write(f">{seq.name} - Copy {i + 1} - {seq.polymerType}\n")
                    f.write(f"{str(seq.sequence)}\n")
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, **kw):
        self.appendCommandLine([self.hklin, self.seqin])
        for path in self.container.inputData.XYZIN_LIST:
            self.appendCommandLine(path)
        self.appendCommandLine("software=refmac")
        return CPluginScript.SUCCEEDED
