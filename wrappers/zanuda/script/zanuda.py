"""
    zanuda.py: CCP4 GUI Project

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

import os
import re
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core.CCP4ModelData import CPdbDataFile
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CMapCoeffsDataFile, CObsDataFile, CPhsDataFile


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
        if error.maxSeverity() > SEVERITY_WARNING:
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
        if error.maxSeverity() > SEVERITY_WARNING:
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

