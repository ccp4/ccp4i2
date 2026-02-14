import os
import re

from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CObsDataFile


class zanuda(CPluginScript):
    TASKTITLE = "Zanuda"
    TASKNAME = "zanuda"
    TASKCOMMAND = "zanuda"
    PERFORMANCECLASS = "CRefinementPerformance"

    def processInputFiles(self):
        self.xyzin = os.path.join(self.getWorkDirectory(), "model.xyz")
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.xyzin)
        miniMtzs = [["F_SIGF", CObsDataFile.CONTENT_FLAG_FMEAN], "FREERFLAG"]
        self.hklin, _, errorReport = self.makeHklin0(miniMtzs)
        return errorReport

    def makeCommandAndScript(self):
        self.appendCommandLine(["xyzin", self.xyzin])
        self.appendCommandLine(["hklin", self.hklin])
        self.appendCommandLine(["xyzout", "zanuda.pdb"])
        self.appendCommandLine(["hklout", "zanuda.mtz"])
        if self.container.controlParameters.AVERAGE:
            self.appendCommandLine(["aver"])
        self.appendCommandLine(["tmpdir", self.workDirectory / "tmpdir"])

    def processOutputFiles(self):
        out = self.container.outputData
        out.XYZOUT.fullPath = self.workDirectory / "zanuda.pdb"
        files = ["FPHIOUT", "DIFFPHIOUT"]
        columns = ["FWT,PHWT", "DELFWT,PHDELWT"]
        error = self.splitHklout(files, columns, self.workDirectory / "zanuda.mtz")
        if error.maxSeverity() > SEVERITY_WARNING:
            return error
        log = self.workDirectory / "log.txt"
        if log.is_file():
            p2 = "^ *[|] *([<> ]{2}) *([0-9]+) *[|]"
            p2 += " *([A-Z][0-9 ]*[0-9]) *[|]"
            p2 += 4 * " *(--|[0-9.]+|error) *[|]" + " *$"
            with log.open(encoding="utf-8") as istream:
                rr = re.findall(p2, istream.read(), re.M)
            if rr and rr[-1][0] == "<<":
                out.PERFORMANCE.RFactor = rr[-1][-2]
                out.PERFORMANCE.RFree = rr[-1][-1]
