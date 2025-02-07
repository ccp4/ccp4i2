#
#  Copyright (C) 2024 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: Martin Maly, David Waterman
#
from core.CCP4PluginScript import CPluginScript
from core.CCP4ErrorHandling import *
import os, glob, shutil

from lxml import etree
from core import CCP4Container
from core import CCP4XtalData
import platform
import json
from math import sqrt
from dxtbx.model.experiment_list import ExperimentList
import gemmi
gemmi.set_leak_warnings(False)
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
        225: {"description": "Unable to find merged MTZ files"},
        226: {"description": "A process in the process pool was terminated abruptly while the future was running or pending.\nOverload of RAM memory occured likely.\nRe-run the job with lower batch size and number of processes."},
    }
    PERFORMANCECLASS = "CDataReductionCCPerformance"
    ASYNCHRONOUS = True
    WHATNEXT = [
        "xia2_ssx_reduce",
        "phaser_pipeline",
        "molrep_pipe",
        "prosmart_refmac"
    ]
    MAINTAINER = "martin.maly@soton.ac.uk"


    # def __init__(self, *args, **kwargs):
    #     self.reference = None
    #     CPluginScript.__init__(self, *args, **kwargs)

    # def processInputFiles(self):
    #     if self.container.inputData.reference.isSet():
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
                if name == "dials_cosym_phil_d_min":
                    val = self.container.controlParameters.dials_cosym_phil_d_min.__str__()
                    phil_file_cosym = os.path.normpath(
                        os.path.join(self.getWorkDirectory(), "cosym_i2.phil")
                    )
                    with open(phil_file_cosym, "w") as f:
                        f.write("d_min={0}\n".format(val))
                    name_xia2 = "symmetry.phil"
                    val_xia2 = "cosym_i2.phil"
                    result.append((name_xia2, val_xia2))
                    continue
                elif name == "MEDIAN_CELL":
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
            # reference and grouping are in self.container.inputData
            if self.container.inputData.reference.isSet():
                name = "reference"
                val = self.container.inputData.reference.fullPath.__str__()
                f.write(name + "={0}\n".format(val))
            if self.container.inputData.grouping.isSet():
                name = "grouping"
                val = self.container.inputData.grouping.fullPath.__str__()
                f.write(name + "={0}\n".format(val))
        self.appendCommandLine(phil_file)

        # Extract integrated DIALS files
        for refl in inp.DIALS_INTEGRATED:
            refl = str(refl)
            expt = refl.rsplit(".refl", 1)[0] + ".expt"

            self.appendCommandLine([f"experiments={expt}", f"reflections={refl}"])
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

        # Locate DataFiles/*.mtz files
        xia2SsxReduceLogPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "xia2.ssx_reduce.log")
        )
        DataFilesPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "DataFiles")
        )
        merged = []
        merged.extend(glob.glob(os.path.normpath(os.path.join(DataFilesPath, "*.mtz"))))
        if not merged:
            element = etree.SubElement(self.xmlroot, "Xia2SsxReduceError")
            with open(xia2SsxReduceLogPath, "r") as xia2SsxReduceLogFile:
                xia2SsxReduceLog = xia2SsxReduceLogFile.read()
                text_terminated_abruptly = \
                    "A process in the process pool was terminated abruptly while the future was running or pending."
            if text_terminated_abruptly in xia2SsxReduceLog:
                self.appendErrorReport(226)
                element.text = text_terminated_abruptly
                element.text += " Overload of RAM memory occured likely. Re-run the job with lower batch size and number of processors."
            else:
                self.appendErrorReport(225)
                element.text = "Unable to find merged MTZ files"
            return CPluginScript.FAILED

        # Remove data_reduction directory to save space
        dataReductionPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "data_reduction")
        )
        if os.path.isdir(dataReductionPath):
            shutil.rmtree(dataReductionPath)

        # Read xia2.ssx_reduce.log
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
        if DialsMergeJsonFilesPath:
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

        # Read LogFiles/dials.cosym_reindex.log or dials.scale.scaled_batch1.log - if exists - it appears when indexing ambiguity
        DialsCosymLogPath1 = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "LogFiles", "dials.cosym_reindex.log")
        )
        DialsCosymLogPath2 = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "LogFiles", "dials.scale.scaled_batch1.log")
        )
        DialsCosymLogPath = None
        if os.path.isfile(DialsCosymLogPath1):
            DialsCosymLogPath = DialsCosymLogPath1
        elif os.path.isfile(DialsCosymLogPath2):
            DialsCosymLogPath = DialsCosymLogPath2
        if DialsCosymLogPath:
            with open(DialsCosymLogPath, "r") as DialsCosymLogFile:
                element = etree.SubElement(self.xmlroot, "DialsCosymLog")
                element.text = etree.CDATA(DialsCosymLogFile.read())
                print(element.text)

        obsOut = self.container.outputData.HKLOUT
        scaledOut = self.container.outputData.DIALS_INTEGRATED

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

        # Locate scaled .refl and .expt files
        scaledPaths = []
        print("\n\n\nglob.iglob(str(Path(DataFilesPath)")
        print(glob.iglob(str(Path(DataFilesPath) / "*.expt")))
        for expt in glob.iglob(str(Path(DataFilesPath) / "*.expt")):
            if not "scaled" in os.path.splitext(os.path.basename(expt))[0]:
                continue
            # Only keep files for which a prefix.{expt,refl} pair exists
            prefix = os.path.splitext(expt)[0]
            refl = prefix + ".refl"
            if not os.path.exists(refl):
                continue
            # print("--- It has .refl " + str(expt))
            # Load the experiments
            el = ExperimentList.from_file(expt, check_format=False)
            if any((not e.scaling_model for e in el)):
                continue
            # Otherwise, looks good. Keep it
            scaledPaths.append(prefix)
        print(str(scaledPaths))
        for scaledPath in scaledPaths:
            scaledDirectory, scaledFilename = os.path.split(scaledPath)
            scaledOut.append(scaledOut.makeItem())
            dest = os.path.join(
                self.getWorkDirectory(),
                os.path.basename(scaledDirectory) + "_" + scaledFilename,
            )
            shutil.copy(scaledPath + ".refl", dest + ".refl")
            shutil.copy(scaledPath + ".expt", dest + ".expt")
            os.remove(scaledPath + ".refl")
            scaledOut[-1].setFullPath(dest + ".refl")
            scaledOut[-1].annotation.set(
                "Reflections: "
                + os.path.basename(scaledDirectory) + "_" + scaledFilename + ".refl"
            )
            self.flushXML()

        # Populate the performance indicator
        spGp = self.container.outputData.PERFORMANCE.spaceGroup
        spGp.set(spGp.fix(str(run_data["space group"]).strip()))

        # self.container.outputData.PERFORMANCE.rMeas.set(
        #     float(str(run_data["Rmeas overall"]).strip())
        # )
        self.container.outputData.PERFORMANCE.ccHalf.set(
            float(str(run_data["CC1/2 overall"]).strip())
        )

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
