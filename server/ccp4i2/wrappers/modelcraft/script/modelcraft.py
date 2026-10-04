import json
import os
import shutil
from ccp4i2.core.CCP4ErrorHandling import SEVERITY_WARNING
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CObsDataFile, CPhsDataFile


def program_xml(result):
    """ModelCraft's report (its modelcraft.json) as program.xml.

    The final model's cycle, residues, waters and R factors, every cycle,
    the resolution and data completeness refinement saw, and why it
    stopped: what says whether the build worked. A value ModelCraft did not
    report is left out.
    """
    from lxml import etree

    root = etree.Element("ModelCraft")

    def add(parent, tag, value):
        if value is not None:
            etree.SubElement(parent, tag).text = str(value)

    add(root, "Version", result.get("version"))
    add(root, "TerminationReason", result.get("termination_reason"))
    fields = ("cycle", "residues", "protein", "nucleic", "waters", "r_work", "r_free")
    final = result.get("final")
    if final:
        node = etree.SubElement(root, "Final")
        for field in fields:
            add(node, field, final.get(field))
    cycles = etree.SubElement(root, "Cycles")
    for cycle in result.get("cycles") or []:
        node = etree.SubElement(cycles, "Cycle")
        for field in fields:
            add(node, field, cycle.get(field))
    refinements = [job for job in result.get("jobs") or [] if "rfree" in job]
    if refinements:
        first = refinements[0]
        add(root, "ResolutionHigh", first.get("resolution_high"))
        add(root, "DataCompleteness", first.get("data_completeness"))
        if "--model" in (result.get("args") or []):
            # Given a model, ModelCraft refines it first: that refinement's
            # starting R is the model's. From phases alone, the first
            # refinement is of what was built, and says nothing of an input.
            node = etree.SubElement(root, "InputModel")
            add(node, "r_work", first.get("initial_rwork"))
            add(node, "r_free", first.get("initial_rfree"))
    return root


class modelcraft(CPluginScript):
    TASKNAME = "modelcraft"
    TASKCOMMAND = "modelcraft"
    ERROR_CODES = {
        201: {"description": "ModelCraft stopped before it produced a model"},
        210: {"description": "Phases are given but ignored while USE_MODEL_PHASES is set"},
    }
    PERFORMANCECLASS = "CRefinementPerformance"
    WHATNEXT = ["coot_rebuild", "coot1"]

    def __init__(self, *args, **kws):
        super(modelcraft, self).__init__(*args, **kws)

    def validity(self):
        error = super(modelcraft, self).validity()
        # PHASES are passed only with USE_MODEL_PHASES off (processInputFiles);
        # a job filled in from an experimental-phasing job has both set, and
        # looks as if it builds on those phases when it does not.
        if (self.container.inputData.PHASES.isSet()
                and bool(self.container.controlParameters.USE_MODEL_PHASES)):
            error.append(
                klass=self.TASKNAME, code=210,
                details="These phases are ignored: ModelCraft takes its phases from the "
                        "model while 'use model phases' is set. Turn it off to build on "
                        "the phases given.",
                name=f"{self.TASKNAME}.container.inputData.PHASES",
                severity=SEVERITY_WARNING)
        return error

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
            self.model = self.container.inputData.XYZIN.getSelectedAtomsFile(
                "model", self.getWorkDirectory())
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
        result = {}
        if os.path.exists(modelcraft_json):
            with open(modelcraft_json) as stream:
                result = json.load(stream)
            from lxml import etree
            with open(self.makeFileName("PROGRAMXML"), "wb") as stream:
                stream.write(etree.tostring(program_xml(result), pretty_print=True))
        if not os.path.exists(modelcraft_cif) or "final" not in result:
            # ModelCraft ends early (no residues built, an incompatible cell)
            # with a reason in its report and exit status 0.
            reason = result.get("termination_reason") or "no report written"
            self.appendErrorReport(201, reason)
            return CPluginScript.FAILED
        shutil.copy(modelcraft_cif, str(outputData.XYZOUT))
        files = ["FPHIOUT", "DIFFPHIOUT", "ABCDOUT"]
        columns = ["FWT,PHWT", "DELFWT,PHDELWT", "HLACOMB,HLBCOMB,HLCCOMB,HLDCOMB"]
        error = self.splitHklout(files, columns, modelcraft_mtz)
        if error.maxSeverity() > SEVERITY_WARNING:
            return CPluginScript.FAILED
        outputData.PERFORMANCE.RFactor.set(result["final"]["r_work"])
        outputData.PERFORMANCE.RFree.set(result["final"]["r_free"])
        return CPluginScript.SUCCEEDED
