"""
    modelcraft.py: CCP4 GUI Project

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

import json
import os
import shutil
from core.CCP4ErrorHandling import SEVERITY_WARNING
from core.CCP4ModelData import CPdbDataFile
from core.CCP4PluginScript import CPluginScript
from core.CCP4XtalData import CMapCoeffsDataFile, CObsDataFile, CPhsDataFile


class modelcraft(CPluginScript):
    TASKMODULE = "model_building"
    TASKTITLE = "ModelCraft"
    TASKNAME = "modelcraft"
    TASKVERSION = 0.1
    TASKCOMMAND = "modelcraft"
    MAINTAINER = "paul.bond@york.ac.uk"
    ERROR_CODES = {
        201: {"description": "Failed to analyse output files"},
        202: {"description": "Failed applying selection to PDB file"},
    }
    PURGESEARCHLIST = [["hklin.mtz", 0], ["log_mtzjoin.txt", 0]]
    PERFORMANCECLASS = "CRefinementPerformance"
    WHATNEXT = ["coot_rebuild"]

    def __init__(self, *args, **kws):
        super(modelcraft, self).__init__(*args, **kws)

    def processInputFiles(self):
        params = self.container.controlParameters
        miniMtzs = [
            ["F_SIGF", CObsDataFile.CONTENT_FLAG_FMEAN],
            ["FREERFLAG", None],
        ]
        if not params.USE_MODEL_PHASES:
            miniMtzs.append(["PHASES", CPhsDataFile.CONTENT_FLAG_HL])
        self.hklin, self.columns, error = self.makeHklin0(miniMtzs)
        if error.maxSeverity() > SEVERITY_WARNING:
            return CPluginScript.FAILED
        self.seqin = os.path.join(self.getWorkDirectory(), "contents.json")
        self.writeContentsJson()
        if self.container.inputData.XYZIN.isSet():
            self.model = os.path.join(self.getWorkDirectory(), "model.xyz")
            self.container.inputData.XYZIN.getSelectedAtomsPdbFile(self.model)
        return CPluginScript.SUCCEEDED

    def writeContentsJson(self):
        params = self.container.controlParameters
        contents = {"copies": 1}
        asu = self.container.inputData.ASUIN
        for seqObj in asu.fileContent.seqList:
            polymer = {
                "sequence": str(seqObj.sequence),
                "stoichiometry": int(seqObj.nCopies),
            }
            if seqObj.polymerType == "PROTEIN" and params.SELENOMET:
                polymer["modifications"] = ["M->MSE"]
            key = {"PROTEIN": "proteins", "RNA": "rnas", "DNA": "dnas"}[
                seqObj.polymerType
            ]
            contents.setdefault(key, []).append(polymer)
        with open(self.seqin, "w") as stream:
            json.dump(contents, stream, indent=4)

    def makeCommandAndScript(self, **kw):
        params = self.container.controlParameters
        self.appendCommandLine(["xray"])
        self.appendCommandLine(["--contents", self.seqin])
        self.appendCommandLine(["--data", self.hklin])
        split_columns = self.columns.split(",")
        fsigf_columns = ",".join(split_columns[:2])
        freer_column = split_columns[2]
        self.appendCommandLine(["--observations", fsigf_columns])
        self.appendCommandLine(["--freerflag", freer_column])
        if not params.USE_MODEL_PHASES:
            abcd_columns = ",".join(split_columns[3:])
            self.appendCommandLine(["--phases", abcd_columns])
            if params.UNBIASED:
                self.appendCommandLine(["--unbiased"])
        if self.container.inputData.XYZIN.isSet():
            self.appendCommandLine(["--model", self.model])
        self.appendCommandLine(["--cycles", params.CYCLES])
        if params.AUTO_STOP:
            self.appendCommandLine(["--auto-stop-cycles", params.STOP_CYCLES])
        else:
            self.appendCommandLine(["--auto-stop-cycles", 0])
        if params.BASIC:
            self.appendCommandLine(["--basic"])
        if params.TWINNED:
            self.appendCommandLine(["--twinned"])
        if not params.SHEETBEND:
            self.appendCommandLine(["--disable-sheetbend"])
        self.appendCommandLine(["--directory", "modelcraft"])
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        directory = os.path.join(self.getWorkDirectory(), "modelcraft")
        modelcraft_cif = os.path.join(directory, "modelcraft.cif")
        modelcraft_mtz = os.path.join(directory, "modelcraft.mtz")
        modelcraft_json = os.path.join(directory, "modelcraft.json")
        XYZOUT_path = os.path.join(self.getWorkDirectory(), "XYZOUT.cif")
        shutil.copy(modelcraft_cif, XYZOUT_path)
        outputData = self.container.outputData
        outputData.XYZOUT.setFullPath(XYZOUT_path)
        outputData.XYZOUT.annotation.set("ModelCraft model")
        outputData.XYZOUT.subType.set(CPdbDataFile.SUBTYPE_MODEL)
        outputData.XYZOUT.contentFlag.set(CPdbDataFile.CONTENT_FLAG_MMCIF)
        outputData.FPHIOUT.annotation.set("ModelCraft best map")
        outputData.FPHIOUT.subType.set(CMapCoeffsDataFile.SUBTYPE_NORMAL)
        outputData.FPHIOUT.contentFlag.set(CMapCoeffsDataFile.CONTENT_FLAG_FPHI)
        outputData.DIFFPHIOUT.annotation.set("Modelcraft difference map")
        outputData.DIFFPHIOUT.subType.set(CMapCoeffsDataFile.SUBTYPE_DIFFERENCE)
        outputData.DIFFPHIOUT.contentFlag.set(CMapCoeffsDataFile.CONTENT_FLAG_FPHI)
        outputData.ABCDOUT.annotation.set("ModelCraft phases")
        outputData.ABCDOUT.subType.set(CPhsDataFile.SUBTYPE_BIASED)
        outputData.ABCDOUT.contentFlag.set(CPhsDataFile.CONTENT_FLAG_HL)
        files = ["FPHIOUT", "DIFFPHIOUT", "ABCDOUT"]
        columns = ["FWT,PHWT", "DELFWT,PHDELWT", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB"]
        error = self.splitHklout(files, columns, modelcraft_mtz)
        if error.maxSeverity() > SEVERITY_WARNING:
            return CPluginScript.FAILED
        with open(modelcraft_json) as stream:
            result = json.load(stream)
        outputData.PERFORMANCE.RFactor.set(result["final"]["r_work"])
        outputData.PERFORMANCE.RFree.set(result["final"]["r_free"])
        return CPluginScript.SUCCEEDED
