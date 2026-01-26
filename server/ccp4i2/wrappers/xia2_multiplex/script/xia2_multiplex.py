import glob
import json
import os
import shutil
from math import sqrt

from dxtbx.model.experiment_list import ExperimentList

from lxml import etree

from ccp4i2.core import CCP4Container, CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript


class Cxia2_multiplex(CPluginScript):

    TASKTITLE = "Data set combination with xia2.multiplex"
    TASKNAME = "xia2_multiplex"
    TASKCOMMAND = "xia2.multiplex"
    TASKMODULE = "data_reduction"
    TASKVERSION = 0.0
    ERROR_CODES = {
        200: {"description": "Failed harvesting integrated data"},
        205: {"description": "Failed parsing xia2.json"},
    }
    PERFORMANCECLASS = "CDataReductionPerformance"
    ASYNCHRONOUS = True
    WHATNEXT = [
        "phaser_pipeline",
        "molrep_pipe",
        "crank2",
        "ShelxCD",
        "ShelxCDE",
    ]
    MAINTAINER = "ccp4@stfc.ac.uk"

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
        self._setCommandLineCore(phil_filename="xia2_multiplex.phil")

        self.xmlroot = etree.Element("Xia2Multiplex")

        self.watchFile(
            os.path.normpath(
                os.path.join(self.getWorkDirectory(), "xia2.multiplex.log")
            ),
            self.handleMultiplexLogChanged,
        )

        return CPluginScript.SUCCEEDED

    def extract_integrated_dials_files(self, datafiles_dir):

        results = []

        # Start by looking for DIALS files from a xia2/dials run
        experiments = glob.glob(os.path.join(datafiles_dir, "*.expt"))
        reflections = glob.glob(os.path.join(datafiles_dir, "*.refl"))

        # Exclude scaled files and sort
        experiments = sorted(e for e in experiments if not e.endswith("_scaled.expt"))
        reflections = sorted(e for e in reflections if not e.endswith("_scaled.refl"))

        for expt, refl in zip(experiments, reflections):
            results.append((expt, refl))

        if results:
            return results

        # If there are no DIALS files, then import integrated MTZ to DIALS files
        self._disable_geometry_refinement = False
        # XXX TODO

        return results

    @staticmethod
    def _extract_data_from_json(json_txt):
        """Get basic data from a xia2.multiplex.json text"""
        d = json.loads(json_txt)

        datasets = d["datasets"]
        stat = datasets["All data"]["merging_stats"]["overall"]

        results = {}
        results["space group"] = "TODO-get from scaled.expt?"
        results["Rmeas overall"] = stat["r_meas"]
        results["High resolution limit"] = 1 / sqrt(stat["d_star_sq_min"])

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
        from ccp4i2.core.CCP4Modules import PROCESSMANAGER

        exitStatus = PROCESSMANAGER().getJobData(
            pid=self.getProcessId(), attribute="exitStatus"
        )
        if exitStatus != CPluginScript.SUCCEEDED:
            element = etree.SubElement(self.xmlroot, "Xia2MultiplexError")
            element.text = "Unknown xia2.multiplex error"
            return exitStatus

        # Read xia2.multiplex.log
        xia2MultiplexLogPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "xia2.multiplex.log")
        )
        if os.path.isfile(xia2MultiplexLogPath):
            with open(xia2MultiplexLogPath, "r") as xia2MultiplexLogFile:
                element = etree.SubElement(self.xmlroot, "Xia2MultiplexLog")
                element.text = etree.CDATA(xia2MultiplexLogFile.read())

        # Read xia2.multiplex.json to read performance
        xia2MultiplexJsonPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "xia2.multiplex.json")
        )
        if os.path.isfile(xia2MultiplexJsonPath):
            with open(xia2MultiplexJsonPath, "r") as xia2MultiplexJsonFile:
                json_txt = xia2MultiplexJsonFile.read()
                run_data = self._extract_data_from_json(json_txt)

        # Read space group and cell from scaled experiments
        scaled = ExperimentList.from_file(
            os.path.normpath(os.path.join(self.getWorkDirectory(), "scaled.expt")),
            check_format=False,
        )
        run_data["space group"] = str(scaled.crystals()[0].get_space_group().info())
        cell = scaled.crystals()[0].get_unit_cell().parameters()
        run_data[
            "unit cell"
        ] = "{:.2f}, {:.2f}, {:.2f}<br/>{:.2f}, {:.2f}, {:.2f}".format(*cell)

        # Also store these in the XML for the report
        element = etree.SubElement(self.xmlroot, "Xia2MultiplexSG")
        element.text = etree.CDATA(run_data["space group"])
        element = etree.SubElement(self.xmlroot, "Xia2MultiplexCell")
        element.text = etree.CDATA(run_data["unit cell"])

        unmergedOut = self.container.outputData.UNMERGEDOUT
        obsOut = self.container.outputData.HKLOUT

        # Find MTZs
        search = [
            dirpath
            for dirpath, _, _ in os.walk(
                self.getWorkDirectory(),
            )
        ]
        candidates = []
        for pth in search:
            candidates.extend(glob.glob(os.path.normpath(os.path.join(pth, "*.mtz"))))
        merged, unmerged = [], []
        for mtz in candidates:
            _, srcFilename = os.path.split(mtz)
            if srcFilename[0].isdigit():
                # exclude files like 3_scaled_unmerged.mtz, which are considered intermediate output
                continue
            if mtz.endswith("unmerged.mtz"):
                unmerged.append(mtz)
            else:
                merged.append(mtz)

        if not unmerged:
            element = etree.SubElement(self.xmlroot, "Xia2MultiplexError")
            element.text = "Unable to find unmerged MTZs"
            return CPluginScript.FAILED
        for candidate in unmerged:

            srcDirectory, srcFilename = os.path.split(candidate)
            unmergedOut.append(unmergedOut.makeItem())
            if os.path.normpath(srcDirectory) != os.path.normpath(
                self.getWorkDirectory()
            ):
                dest = os.path.join(
                    self.getWorkDirectory(),
                    os.path.basename(srcDirectory) + "_" + srcFilename,
                )
                shutil.copy(candidate, dest)
                unmergedOut[-1].setFullPath(dest)
                unmergedOut[-1].annotation.set(
                    "Unmerged reflections: "
                    + os.path.basename(srcDirectory)
                    + "_"
                    + srcFilename
                )
            else:
                unmergedOut[-1].setFullPath(candidate)
                unmergedOut[-1].annotation.set("Unmerged reflections: " + srcFilename)

        if not merged:
            element = etree.SubElement(self.xmlroot, "Xia2MultiplexError")
            element.text = "Unable to find merged MTZs"
            return CPluginScript.FAILED
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

        self.container.outputData.PERFORMANCE.highResLimit.set(
            float(str(run_data["High resolution limit"]).strip())
        )

        return CPluginScript.SUCCEEDED

    def handleMultiplexLogChanged(self, filename):
        # remove xia2.multiplex.log nodes
        for Xia2MultiplexLogNode in self.xmlroot.xpath("Xia2MultiplexLog"):
            self.xmlroot.remove(Xia2MultiplexLogNode)
        xia2MultiplexLogNode = etree.SubElement(self.xmlroot, "Xia2MultiplexLog")
        with open(filename, "r") as xia2MultiplexLogFile:
            xia2MultiplexLogNode.text = etree.CDATA(xia2MultiplexLogFile.read())
        self.flushXML()

    def flushXML(self):
        tmpFilename = self.makeFileName("PROGRAMXML") + "_tmp"
        with open(tmpFilename, "wb") as xmlFile:
            xmlFile.write(etree.tostring(self.xmlroot, pretty_print=True))
        if os.path.exists(self.makeFileName("PROGRAMXML")):
            os.remove(self.makeFileName("PROGRAMXML"))
        os.rename(tmpFilename, self.makeFileName("PROGRAMXML"))
