import json
import os
import shutil
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CObsDataFile, CPhsDataFile


class modelcraft(CPluginScript):
    TASKMODULE = "model_building"
    TASKTITLE = "AutoBuild with ModelCraft, Buccaneer and Nautilus"
    SHORTTASKTITLE = "ModelCraft"
    TASKNAME = "modelcraft"
    TASKCOMMAND = "modelcraft"
    MAINTAINER = "paul.bond@york.ac.uk"
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
            # Convert CString to str for dictionary key lookup
            key = {"PROTEIN": "proteins", "RNA": "rnas", "DNA": "dnas"}[
                str(seqObj.polymerType)
            ]
            contents.setdefault(key, []).append(polymer)
        with open(self.seqin, "w") as stream:
            json.dump(contents, stream, indent=4)

    def makeCommandAndScript(self):
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
        if not params.BASIC and not params.PRUNING:
            self.appendCommandLine(["--disable-pruning"])
        if not params.PARROT:
            self.appendCommandLine(["--disable-parrot"])
        if not params.BASIC and not params.DUMMY_ATOMS:
            self.appendCommandLine(["--disable-dummy-atoms"])
        if not params.BASIC and not params.WATERS:
            self.appendCommandLine(["--disable-waters"])
        if not params.BASIC and not params.SIDE_CHAIN_FIXING:
            self.appendCommandLine(["--disable-side-chain-fixing"])
        self.appendCommandLine(["--directory", "modelcraft"])
        return CPluginScript.SUCCEEDED

    def processOutputFiles(self):
        directory = os.path.join(self.getWorkDirectory(), "modelcraft")
        modelcraft_cif = os.path.join(directory, "modelcraft.cif")
        modelcraft_mtz = os.path.join(directory, "modelcraft.mtz")
        modelcraft_json = os.path.join(directory, "modelcraft.json")
        outputData = self.container.outputData
        shutil.copy(modelcraft_cif, str(outputData.XYZOUT))
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
