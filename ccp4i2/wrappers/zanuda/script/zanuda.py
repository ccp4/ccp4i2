import os
import re

from ....core.CCP4ErrorHandling import Severity
from ....core.CCP4ModelData import CPdbDataFile
from ....core.CCP4PluginScript import CPluginScript
from ....core.CCP4XtalData import CMapCoeffsDataFile, CObsDataFile


class zanuda(CPluginScript):
    TASKMODULE = "refinement"
    TASKTITLE = "Zanuda"
    TASKNAME = "zanuda"
    TASKVERSION = 0.1
    TASKCOMMAND = "zanuda"
    MAINTAINER = "andrey.lebedev@stfc.ac.uk"
    ERROR_CODES = {
        201: {"description": "Failed to analyse output files"},
        202: {"description": "Failed applying selection to PDB file"},
    }
    PURGESEARCHLIST = [["hklin.mtz", 0], ["log_mtzjoin.txt", 0]]
    PERFORMANCECLASS = "CRefinementPerformance"

    def __init__(self, *args, **kws):
        super(zanuda, self).__init__(*args, **kws)

    def processInputFiles(self):
        miniMtzs = [
            ["F_SIGF", CObsDataFile.CONTENT_FLAG_FMEAN],
            ["FREERFLAG", None],
        ]
        self.hklin, self.columns, error = self.makeHklin0(miniMtzs)
        if error.maxSeverity() > Severity.WARNING:
            return CPluginScript.FAILED
        self.model = os.path.join(self.getWorkDirectory(), "model.xyz")
        self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.model)
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self, **kw):
        self.appendCommandLine(["xyzin", self.model])
        self.appendCommandLine(["hklin", self.hklin])
        self.appendCommandLine(["xyzout", "zanuda.pdb"])
        self.appendCommandLine(["hklout", "zanuda.mtz"])
        params = self.container.controlParameters
        if params.AVERAGE:
          self.appendCommandLine(["aver"])
        tmpdir = os.path.join(self.getWorkDirectory(), "tmpdir")
        self.appendCommandLine(["tmpdir", tmpdir])
        # for tests (faster):
        # self.appendCommandLine(["notwin"])

    def processOutputFiles(self):
        directory = self.getWorkDirectory()
        zanuda_pdb = os.path.join(directory, "zanuda.pdb")
        zanuda_mtz = os.path.join(directory, "zanuda.mtz")
        outputData = self.container.outputData
        outputData.XYZOUT.setFullPath(zanuda_pdb)
        outputData.XYZOUT.annotation.set("Zanuda model")
        outputData.XYZOUT.subType.set(CPdbDataFile.SUBTYPE_MODEL)
        outputData.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_PDB)
        outputData.FPHIOUT.annotation.set("Zanuda best map")
        outputData.FPHIOUT.subType.set(CMapCoeffsDataFile.SUBTYPE_NORMAL)
        outputData.FPHIOUT.contentFlag.set(CMapCoeffsDataFile.CONTENT_FLAG_FPHI)
        outputData.DIFFPHIOUT.annotation.set("Zanuda difference map")
        outputData.DIFFPHIOUT.subType.set(CMapCoeffsDataFile.SUBTYPE_DIFFERENCE)
        outputData.DIFFPHIOUT.contentFlag.set(CMapCoeffsDataFile.CONTENT_FLAG_FPHI)
        files = ["FPHIOUT", "DIFFPHIOUT"]
        columns = ["FWT,PHWT", "DELFWT,PHDELWT"]
        error = self.splitHklout(files, columns, zanuda_mtz)
        if error.maxSeverity() > Severity.WARNING:
            return CPluginScript.FAILED
        log = os.path.join(self.getWorkDirectory(), 'log.txt')
        if os.path.isfile(log):
            p2 = '^ *[|] *([<> ]{2}) *([0-9]+) *[|]'
            p2 += ' *([A-Z][0-9 ]*[0-9]) *[|]'
            p2 += 4* ' *(--|[0-9.]+|error) *[|]' + ' *$'
            with open(log) as istream:
                rr = re.findall(p2, istream.read(), re.M)
            if rr and rr[-1][0] == '<<':
                outputData.PERFORMANCE.RFactor.set(rr[-1][-2])
                outputData.PERFORMANCE.RFree.set(rr[-1][-1])
        return CPluginScript.SUCCEEDED
