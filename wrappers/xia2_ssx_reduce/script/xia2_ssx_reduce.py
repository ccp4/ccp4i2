#
#  Copyright (C) 2024 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: Martin Maly, David Waterman
#
# TO DO
# =====
# errors when files missing
# Stuart reference path
# Stuart verdict Rmeas
# clean comments
#  ***************  new performance classes also need the keytype to be registered with the database **************************
#   See the definition of KEYTYPELIST in dbapi/CCP4DbApi.py
from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import *
import os, glob, shutil

# from core import CCP4Utils
from lxml import etree
from core import CCP4Container
from core import CCP4XtalData
import platform
import json
from math import sqrt
from dxtbx.model.experiment_list import ExperimentList
import gemmi
import glob
from pathlib import Path
import shutil


class Cxia2_ssx_reduce(CPluginScript):

    TASKTITLE = "Reduction of serial datasets using xia2.ssx_reduce"
    TASKNAME = "xia2_ssx_reduce"
    TASKCOMMAND = "xia2.ssx_reduce"
    if platform.system() == "Windows":
        TASKCOMMAND = "xia2.ssx_reduce.bat"
    TASKMODULE = "data_reduction"
    TASKVERSION = 0.0
    ERROR_CODES = {
        200: {"description": "Failed harvesting integrated data"},
    }
    PERFORMANCECLASS = "CDataReductionPerformance"
    ASYNCHRONOUS = True
    WHATNEXT = [
        "phaser_pipeline",
        "molrep_pipe",
        "prosmart_refmac"
    ]
    MAINTAINER = "martin.maly@soton.ac.uk"


    # def __init__(self, *args, **kwargs):
    #     self.reference = None
    #     CPluginScript.__init__(self, *args, **kwargs)

    # def processInputFiles(self):
    #     # if self.container.inputData.reference.isSet():
    #     #     self.reference = self.container.inputData.reference.fullPath.__str__()
    #     if self.container.controlParameters.reference.isSet():
    #         self.reference = self.container.controlParameters.reference.fullPath.__str__()

    def extract_parameters(self, container):
        """Walk through a container locating parameters that have been set
        and return a list of name, value pairs"""

        result = []
        dataOrder = container.dataOrder()
        contents = [getattr(container, name) for name in dataOrder]
        for model in contents:
            if isinstance(model, CCP4Container.CContainer):
                result.extend(self.extract_parameters(model))
            elif model.isSet():
                name = model.objectName().replace("__", ".")
                if name == "reference":                   # MM
                    # val = self.reference
                    val = self.container.controlParameters.reference.fullPath.__str__()
                    result.append((name, val))
                    continue
                if name == "dials_cosym_phil_d_min":      # MM
                    phil_file_cosym = os.path.normpath(
                        os.path.join(self.getWorkDirectory(), "cosym_i2.phil")
                    )
                    with open(phil_file_cosym, "w") as f:
                        f.write("d_min={0}\n".format(val))
                    name_xia2 = "symmetry.phil"
                    val_xia2 = "cosym_i2.phil"
                    result.append((name_xia2, val_xia2))
                    continue
                # ensure commas are converted to whitespace-separated lists.
                # Only whitespace appears to work correctly with PHIL multiple
                # choice definitions.
                val = str(model.get()).split()
                val = " ".join([v[:-1] if v.endswith(",") else v for v in val])
                result.append((name, val))
        return result

    def _setCommandLineCore(self, phil_filename):
        par = self.container.controlParameters
        inp = self.container.inputData

        # PHIL parameters set by the gui
        phil_file = os.path.normpath(
            os.path.join(self.getWorkDirectory(), phil_filename)
        )
        with open(phil_file, "w") as f:
            for (name, val) in self.extract_parameters(par):
                f.write(name + "={0}\n".format(val))
        self.appendCommandLine(phil_file)

        # Extract integrated DIALS files
        for refl in inp.DIALS_INTEGRATED:
            refl = str(refl)
            expt = refl.rsplit(".refl", 1)[0] + ".expt"

            self.appendCommandLine([f"experiments={expt}", f"reflections={refl}"])

        # Disable unit cell refinement by two_theta_refine, as the geometry
        # derived from MTZ files is suspect.
        # if self._disable_geometry_refinement:
        #    self.appendCommandLine(["unit_cell.refine=None",])

        return

    def makeCommandAndScript(self):

        # Create PHIL file and command line
        self._setCommandLineCore(phil_filename="xia2_ssx_reduce.phil")

        self.xmlroot = etree.Element("Xia2SsxReduce")

        self.watchFile(
            os.path.normpath(
                os.path.join(self.getWorkDirectory(), "xia2.ssx_reduce.log")
            ),
            self.handleSsxReduceLogChanged,
        )

        return CPluginScript.SUCCEEDED

    # def extract_integrated_dials_files(self, datafiles_dir):

    #     results = []

    #     # Start by looking for DIALS files from a xia2/dials run
    #     experiments = glob.glob(os.path.join(datafiles_dir, "*.expt"))
    #     reflections = glob.glob(os.path.join(datafiles_dir, "*.refl"))

    #     # Exclude scaled files and sort
    #     experiments = sorted(e for e in experiments if not e.endswith("_scaled.expt"))
    #     reflections = sorted(e for e in reflections if not e.endswith("_scaled.refl"))

    #     for expt, refl in zip(experiments, reflections):
    #         results.append((expt, refl))

    #     if results:
    #          return results

    #     # If there are no DIALS files, then import integrated MTZ to DIALS files
    #     self._disable_geometry_refinement = False
    #     # XXX TODO

    #     return results

    @staticmethod
    def _extract_data_from_json(json_txt):
        """Get basic data from a LogFiles/dials.merge.json text"""
        dataset = json.loads(json_txt)
        json_root = list(dataset)[0] # usually wavelength
        stat = dataset[json_root]["merging_stats"]["overall"] # MM
        results = {}
        results["Rmeas overall"] = stat["r_meas"]
        results["CC1/2 overall"] = stat["cc_one_half"]
        results["High resolution limit"] = 1 / sqrt(stat["d_star_sq_min"])
        # results["space group"] is get from merged.mtz
        # results["cell"]        is get from merged.mtz
        return results

    @staticmethod
    def _extract_col_name_type(mtz_filename):
        """Get the list of column name and type pairs from an MTZ file"""
        mtz = CCP4XtalData.CMtzDataFile(mtz_filename)
        return [
            (str(e.columnLabel), str(e.columnType))
            for e in mtz.fileContent.getListOfColumns()
        ]

    def processOutputFiles(self):

        # Check for exit status of the program
        from core.CCP4Modules import PROCESSMANAGER

        exitStatus = PROCESSMANAGER().getJobData(
            pid=self.getProcessId(), attribute="exitStatus"
        )
        if exitStatus != CPluginScript.SUCCEEDED:
            element = etree.SubElement(self.xmlroot, "Xia2SsxReduceError")
            element.text = "Unknown xia2.ssx_reduce error"
            return exitStatus

        # Remove data_reduction directory to save space
        dataReductionPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "data_reduction")
        )
        if os.path.isdir(dataReductionPath):
            shutil.rmtree(dataReductionPath)

        # Read xia2.ssx_reduce.log
        xia2SsxReduceLogPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "xia2.ssx_reduce.log")
        )
        if os.path.isfile(xia2SsxReduceLogPath):
            with open(xia2SsxReduceLogPath, "r") as xia2SsxReduceLogFile:
                element = etree.SubElement(self.xmlroot, "Xia2SsxReduceLog")
                element.text = etree.CDATA(xia2SsxReduceLogFile.read())
                print(element.text)

        # Locate LogFiles/dials.merge*.{json,log} files
        DialsMergeJsonFilesPath = []
        DialsMergeLogFilesPath = []
        LogFilesPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "LogFiles")
        )
        for JsonFilePath in glob.iglob(str(Path(LogFilesPath) / "*.json")):
            DialsMergePrefix = os.path.splitext(os.path.basename(JsonFilePath))[0]
            if "dials.merge" in DialsMergePrefix:
                DialsMergeJsonFilePath = os.path.normpath(
                    os.path.join(self.getWorkDirectory(), "LogFiles", DialsMergePrefix + ".json")
                )
                if os.path.isfile(DialsMergeJsonFilePath):
                    DialsMergeJsonFilesPath.append(DialsMergeJsonFilePath)
                DialsMergeLogFilePath = os.path.normpath(
                    os.path.join(self.getWorkDirectory(), "LogFiles", DialsMergePrefix + ".log")
                )
                if os.path.isfile(DialsMergeLogFilePath):
                    DialsMergeLogFilesPath.append(DialsMergeLogFilePath)

        # Read LogFiles/dials.merge*.log files
        element_master = etree.SubElement(self.xmlroot, "DialsMergeLogMaster")
        for DialsMergeLogPath in DialsMergeLogFilesPath:
            if os.path.isfile(DialsMergeLogPath):
                with open(DialsMergeLogPath, "r") as DialsMergeLogFile:
                    element = etree.SubElement(element_master, "DialsMergeLog")
                    element.text = etree.CDATA(DialsMergeLogFile.read())
                    print(element.text)

        # Read the 1st LogFiles/dials.merge*.json file to read performance
        xia2SsxReduceJsonPath = DialsMergeJsonFilesPath[0]
        if os.path.isfile(xia2SsxReduceJsonPath):
            with open(xia2SsxReduceJsonPath, "r") as xia2SsxReduceJsonFile:
                json_txt = xia2SsxReduceJsonFile.read()
                run_data = self._extract_data_from_json(json_txt)
        else:
            run_data = {}
            run_data["space group"] = ""
            run_data["unit cell"] = ""
            # run_data["Rmeas overall"] = 0
            run_data["CC1/2 overall"] = 0
            run_data["High resolution limit"] = 0

        # Read LogFiles/dials.cosym_reindex.log - if exists - it appears when indexing ambiguity
        DialsCosymLogPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "LogFiles", "dials.cosym_reindex.log")
        )
        if os.path.isfile(DialsCosymLogPath):
            with open(DialsCosymLogPath, "r") as DialsCosymLogFile:
                element = etree.SubElement(self.xmlroot, "DialsCosymLog")
                element.text = etree.CDATA(DialsCosymLogFile.read())
                print(element.text)

        obsOut = self.container.outputData.HKLOUT

        DataFilesPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "DataFiles")
        )
        merged = []
        merged.extend(glob.glob(os.path.normpath(os.path.join(DataFilesPath, "*.mtz"))))

        # Find MTZs
        #search = [
        #    dirpath
        #    for dirpath, _, _ in os.walk(
        #        self.getWorkDirectory(),
        #    )
        #]
        #candidates = []
        #for pth in search:
        #    candidates.extend(glob.glob(os.path.normpath(os.path.join(pth, "*.mtz"))))
        #merged = []
        #for mtz in candidates:
        #    _, srcFilename = os.path.split(mtz)
        #    if srcFilename[0].isdigit():
        #        # exclude files like 3_scaled_unmerged.mtz, which are considered intermediate output
        #        continue
        #    # if mtz.endswith("unmerged.mtz"): # MM
        #    #     unmerged.append(mtz)  # MM
        #    else:
        #        merged.append(mtz)

        #if not unmerged:
        #    element = etree.SubElement(self.xmlroot, "Xia2SsxReduceError")
        #    element.text = "Unable to find unmerged MTZs"
        #    return CPluginScript.FAILED
        #for candidate in unmerged:

        #    srcDirectory, srcFilename = os.path.split(candidate)
        #    unmergedOut.append(unmergedOut.makeItem())
        #    if os.path.normpath(srcDirectory) != os.path.normpath(
        #        self.getWorkDirectory()
        #    ):
        #        dest = os.path.join(
        #            self.getWorkDirectory(),
        #            os.path.basename(srcDirectory) + "_" + srcFilename,
        #        )
        #        shutil.copy(candidate, dest)
        #        unmergedOut[-1].setFullPath(dest)
        #        unmergedOut[-1].annotation.set(
        #            "Unmerged reflections: "
        #            + os.path.basename(srcDirectory)
        #            + "_"
        #            + srcFilename
        #        )
        #    else:
        #        unmergedOut[-1].setFullPath(candidate)
        #        unmergedOut[-1].annotation.set("Unmerged reflections: " + srcFilename)

        if not merged:
            element = etree.SubElement(self.xmlroot, "Xia2SsxReduceError")
            element.text = "Unable to find merged MTZs"
            return CPluginScript.FAILED

        # Read space group and cell from the 1st merged MTZ to complete performance info
        MergedMtzPath = merged[0]
        if os.path.isfile(MergedMtzPath):
            mtz = gemmi.read_mtz_file(MergedMtzPath)
            run_data["space group"] = mtz.spacegroup.hm
            cell = [mtz.cell.a, mtz.cell.b, mtz.cell.c, mtz.cell.alpha, mtz.cell.beta, mtz.cell.gamma]
            run_data[
                "unit cell"
            ] = "{:.2f}, {:.2f}, {:.2f}<br/>{:.2f}, {:.2f}, {:.2f}".format(*cell)

        # Also store these in the XML for the report
        element = etree.SubElement(self.xmlroot, "Xia2SsxReduceSG")
        element.text = etree.CDATA(run_data["space group"])
        element = etree.SubElement(self.xmlroot, "Xia2SsxReduceCell")
        element.text = etree.CDATA(run_data["unit cell"])

        for srcPath in merged:
            srcDirectory, srcFilename = os.path.split(srcPath)

            srcFilename = os.path.split(srcPath)[1]

            # Prepare for splitMtz
            col_names, col_types = zip(*self._extract_col_name_type(srcPath))
            anomalous = "K" in col_types
            if anomalous:
                colin = "I(+),SIGI(+),I(-),SIGI(-)"
                colout = "Iplus,SIGIplus,Iminus,SIGIminus"
            else:
                colin = "IMEAN,SIGIMEAN"
                colout = "I,SIGI"

            logFile = os.path.join(self.getWorkDirectory(), "cmtzsplit.log")
            out_triplet = [srcPath, colin, colout]
            status = self.splitMtz(srcPath, [out_triplet], logFile)
            if status == CPluginScript.SUCCEEDED:
                obsOut.append(obsOut.makeItem())

                if os.path.normpath(srcDirectory) != os.path.normpath(
                    self.getWorkDirectory()
                ):
                    dest = os.path.join(
                        self.getWorkDirectory(),
                        os.path.basename(srcDirectory) + "_" + srcFilename,
                    )
                    shutil.copy(srcPath, dest)
                    obsOut[-1].setFullPath(dest)
                    obsOut[-1].annotation.set(
                        "Reflections: "
                        + os.path.basename(srcDirectory)
                        + "_"
                        + srcFilename
                    )
                else:
                    obsOut[-1].setFullPath(srcPath)
                    obsOut[-1].annotation.set("Reflections: " + srcFilename)

                if anomalous:
                    flag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR
                    obsOut[-1].contentFlag = flag
                else:
                    flag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN
                    obsOut[-1].contentFlag = flag
            else:
                self.appendErrorReport(200)
                return CPluginScript.FAILED

            self.flushXML()

        # Populate the performance indicator
        spGp = self.container.outputData.PERFORMANCE.spaceGroup
        spGp.set(spGp.fix(str(run_data["space group"]).strip()))

        self.container.outputData.PERFORMANCE.rMeas.set(
            float(str(run_data["Rmeas overall"]).strip())
        )
        #self.container.outputData.PERFORMANCE.ccHalf.set(
        #    float(str(run_data["CC1/2 overall"]).strip())
        #)

        self.container.outputData.PERFORMANCE.highResLimit.set(
            float(str(run_data["High resolution limit"]).strip())
        )

        return CPluginScript.SUCCEEDED

    def handleSsxReduceLogChanged(self, filename):
        # remove xia2.ssx_reduce.log nodes
        for Xia2SsxReduceLogNode in self.xmlroot.xpath("Xia2SsxReduceLog"):
            self.xmlroot.remove(Xia2SsxReduceLogNode)
        xia2SsxReduceLogNode = etree.SubElement(self.xmlroot, "Xia2SsxReduceLog")
        with open(filename, "r") as xia2SsxReduceLogFile:
            xia2SsxReduceLogNode.text = etree.CDATA(xia2SsxReduceLogFile.read())
        self.flushXML()

    def flushXML(self):
        tmpFilename = self.makeFileName("PROGRAMXML") + "_tmp"
        with open(tmpFilename, "wb") as xmlFile:
            xmlFile.write(etree.tostring(self.xmlroot, pretty_print=True))
        if os.path.exists(self.makeFileName("PROGRAMXML")):
            os.remove(self.makeFileName("PROGRAMXML"))
        os.rename(tmpFilename, self.makeFileName("PROGRAMXML"))
