#
#  Copyright (C) 2016 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#  Acknowledgements: based on code by Graeme Winter and Martin Noble.
#

import glob
import json
import os
import platform
import re
import shutil

from lxml import etree

from ....core import CCP4Container
from ....core import CCP4XtalData
from ....core.CCP4Modules import PROCESSMANAGER
from ....core.CCP4PluginScript import CPluginScript


class Cxia2_dials(CPluginScript):

    TASKTITLE = "Data processing with xia2/dials"
    TASKNAME = "xia2_dials"
    TASKCOMMAND = "xia2"
    if platform.system() == "Windows":
        TASKCOMMAND = "xia2.exe"
    TASKMODULE = "data_processing"
    TASKVERSION = 0.0
    ERROR_CODES = {
        200: {"description": "Failed harvesting integrated data"},
        201: {"description": "Failed scaled data"},
        202: {"description": "Failed harvesting pointless XML"},
        203: {"description": "Failed harvesting aimless xml"},
        204: {"description": "Failed harvesting truncate xml"},
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
    ]  # , 'dials_image', 'dials_rlattice']
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
        """Parts of makeCommandAndScript shared by xia2/dials and xia2/xds"""
        par = self.container.controlParameters
        inp = self.container.inputData

        # PHIL parameters set by the gui
        phil_file = os.path.normpath(
            os.path.join(self.getWorkDirectory(), phil_filename)
        )
        with open(phil_file, "w") as f:
            for (name, val) in self.extract_parameters(par):
                f.write(name + "={0}\n".format(val))
        self.appendCommandLine([phil_file])

        # Data location as image files
        for e in inp.IMAGE_FILE:
            im_file = str(e.imageFile).strip()
            start = e.imageStart
            end = e.imageEnd
            postfix = ""
            if start and end:
                postfix = ":{0}:{1}".format(start, end)
            if im_file:
                self.appendCommandLine('image="{0}{1}"'.format(im_file, postfix))

        # Data location as a directory
        if inp.IMAGE_DIRECTORY:
            # FIXME: quoting the image directory in order to deal with
            # whitespace in paths causes xia2 to fail! It is not obvious how to
            # work around this, so for now remove this change:

            # self.appendCommandLine(['"%s"' % str(inp.IMAGE_DIRECTORY)])

            # and replace with the unquoted version, which at least only fails
            # when the path does contain whitespace
            self.appendCommandLine(["%s" % str(inp.IMAGE_DIRECTORY)])

        return

    def makeCommandAndScript(self):
        par = self.container.controlParameters
        inp = self.container.inputData

        # Set xia2 switches
        self.appendCommandLine(
            [
                "pipeline=dials",
            ]
        )

        # Create PHIL file and command line
        self._setCommandLineCore(phil_filename="xia2_dials.phil")

        self.xmlroot = etree.Element("Xia2Dials")

        self.watchFile(
            os.path.normpath(os.path.join(self.getWorkDirectory(), "xia2.txt")),
            self.handleXia2DotTxtChanged,
        )

        return CPluginScript.SUCCEEDED

    @staticmethod
    def _get_annotation(prefix, suffix):
        """Form suitable annotation strings"""
        return prefix + " from DIALS integration of " + suffix

    @staticmethod
    def _extract_data_from_json(json_txt):
        """Get basic data from a xia2.json text"""
        d = json.loads(json_txt)
        xls = d["_crystals"]
        crystal_names = [str(e) for e in xls.keys()]
        results = []
        for name in crystal_names:
            dic = {"crystal_name": name}
            xl = xls[name]
            wl = xl["_wavelengths"]
            dic["wavelengths"] = [str(e) for e in wl.keys()]
            dic["space group"] = xl["_scaler"]["_scalr_likely_spacegroups"][0]
            stats = list(xl["_scaler"]["_scalr_statistics"].values())[0]
            dic["Rmeas overall"] = stats["Rmeas(I)"][0]
            dic["High resolution limit"] = stats["High resolution limit"][0]
            results.append(dic)
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
        exitStatus = PROCESSMANAGER().getJobData(
            pid=self.getProcessId(), attribute="exitStatus"
        )
        if exitStatus != CPluginScript.SUCCEEDED:
            element = etree.SubElement(self.xmlroot, "Xia2Error")
            element.text = "Failed to locate XIA2"
            return exitStatus

        # Read xia2.txt
        xia2TxtPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "xia2.txt")
        )
        if os.path.isfile(xia2TxtPath):
            with open(xia2TxtPath, "r") as xia2TxtFile:
                element = etree.SubElement(self.xmlroot, "Xia2Txt")
                element.text = etree.CDATA(xia2TxtFile.read())

        # Infer if xia2 gave an error by virtue of xia2.error existing
        xia2ErrorPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "xia2.error")
        )
        if os.path.isfile(xia2ErrorPath):
            with open(xia2ErrorPath, "r") as xia2ErrorFile:
                element = etree.SubElement(self.xmlroot, "Xia2Error")
                element.text = etree.CDATA(xia2ErrorFile.read())
                self.flushXML()
            return CPluginScript.SUCCEEDED

        # Read xia2.json to extract crystal and wavelength names. The
        # wavelength names will end up appended to column names if >1
        # wavelength. Will need for splitMtz. Also read performance
        # statistics from this file.
        wavelength_names = []
        crystal_name = "DEFAULT"
        performance_stats = None
        xia2JsonPath = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "xia2.json")
        )
        if os.path.isfile(xia2JsonPath):
            with open(xia2JsonPath, "r") as xia2JsonFile:
                json_txt = xia2JsonFile.read()
                json_data = self._extract_data_from_json(json_txt)
                if len(json_data) != 1:
                    self.appendErrorReport(200)
                    return CPluginScript.FAILED
                crystal_name = json_data[0]["crystal_name"]
                wavelength_names = json_data[0]["wavelengths"]
        element = etree.SubElement(self.xmlroot, "Xia2CrystalName")
        element.text = etree.CDATA(crystal_name)

        par = self.container.controlParameters
        tmp = par.xia2.xia2__settings.xia2__settings__input
        anomalous = (
            tmp.xia2__settings__input__atom.isSet()
            or str(tmp.xia2__settings__input__anomalous) == "True"
        )

        unmergedOut = self.container.outputData.UNMERGEDOUT
        obsOut = self.container.outputData.HKLOUT
        freerOut = self.container.outputData.FREEROUT

        # Grab integrated data
        candidates = glob.glob(
            os.path.normpath(
                os.path.join(self.getWorkDirectory(), "DataFiles", "*INTEGRATE.mtz")
            )
        )
        for candidateIntegratedFile in candidates:
            srdDirectory, srcFilename = os.path.split(candidateIntegratedFile)
            destPath = os.path.normpath(
                os.path.join(self.getWorkDirectory(), srcFilename)
            )
            shutil.copyfile(candidateIntegratedFile, destPath)
            unmergedOut.append(unmergedOut.makeItem())
            unmergedOut[-1].fullPath = destPath
            anno = self._get_annotation("Unmerged reflections", srcFilename[:-13])
            unmergedOut[-1].annotation = anno

        # Grab merged files
        pattern = os.path.normpath(
            os.path.join(self.getWorkDirectory(), "DataFiles", "*free.mtz")
        )
        possibleFilesToCopy = glob.glob(pattern)
        for srcPath in possibleFilesToCopy:
            srcDirectory, srcFilename = os.path.split(srcPath)

            if len(wavelength_names) > 1:
                obsPath_list = []
                for w in wavelength_names:
                    obsPath = os.path.join(
                        self.getWorkDirectory(),
                        srcFilename[:-9] + "_obs_{0}.mtz".format(w),
                    )
                    obsPath_list.append(obsPath)
            else:
                obsPath_list = [
                    os.path.join(self.getWorkDirectory(), srcFilename[:-9] + "_obs.mtz")
                ]
            freerPath = os.path.join(
                self.getWorkDirectory(), srcFilename[:-9] + "_freer.mtz"
            )

            # Prepare for splitMtz
            col_names, col_types = zip(*self._extract_col_name_type(srcPath))
            anomalous = "K" in col_types
            if anomalous:
                colin_base = "I(+){0},SIGI(+){0},I(-){0},SIGI(-){0}"
                if len(wavelength_names) > 1:
                    colin_list = [colin_base.format("_" + w) for w in wavelength_names]
                else:
                    colin_list = [colin_base.format("")]
                colout = "Iplus,SIGIplus,Iminus,SIGIminus"
            else:
                colin_list = ["IMEAN,SIGIMEAN"]
                colout = "I,SIGI"
            colfree = "FreeR_flag"
            colfreeout = "FREER"
            logFile = os.path.join(self.getWorkDirectory(), "cmtzsplit.log")
            out_triplets = [
                [obspth, colin, colout]
                for obspth, colin in zip(obsPath_list, colin_list)
            ]
            out_triplets.append([freerPath, colfree, colfreeout])
            status = self.splitMtz(srcPath, out_triplets, logFile)
            if status == CPluginScript.SUCCEEDED:
                for w, obsPath in zip(wavelength_names, obsPath_list):
                    obsOut.append(obsOut.makeItem())
                    obsOut[-1].fullPath = obsPath
                    anno = self._get_annotation(
                        "Reflections: {0}".format(w), srcFilename[:-8]
                    )
                    obsOut[-1].annotation = anno
                    if anomalous:
                        flag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_IPAIR
                        obsOut[-1].contentFlag = flag
                    else:
                        flag = CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN
                        obsOut[-1].contentFlag = flag
                freerOut.append(freerOut.makeItem())
                freerOut[-1].fullPath = freerPath
                anno = self._get_annotation("FreeR", srcFilename[:-8])
                freerOut[-1].annotation = anno
            else:
                self.appendErrorReport(200)
                return CPluginScript.FAILED

            self.flushXML()

            # Grab all log files
            allLogs = glob.glob(
                os.path.normpath(
                    os.path.join(self.getWorkDirectory(), "LogFiles", "*.log")
                )
            )
            for logFilePath in allLogs:
                destLogPath = os.path.normpath(
                    os.path.join(self.getWorkDirectory(), os.path.split(logFilePath)[1])
                )
                shutil.copyfile(logFile, destLogPath)

        self._collect_pickles_and_jsons()

        # Populate the performance indicator
        spGp = self.container.outputData.PERFORMANCE.spaceGroup
        spGp.set(spGp.fix(str(json_data[0]["space group"]).strip()))

        self.container.outputData.PERFORMANCE.rMeas.set(
            float(str(json_data[0]["Rmeas overall"]).strip())
        )

        self.container.outputData.PERFORMANCE.highResLimit.set(
            float(str(json_data[0]["High resolution limit"]).strip())
        )

        return CPluginScript.SUCCEEDED

    def _collect_pickles_and_jsons(self):
        """Copy DIALS pickle and json files for the viewers"""
        # KJS : Get the name of the right directories (extracted from the names of the logfiles)
        outContJ = self.container.outputData.DIALSJOUT
        outContP = self.container.outputData.DIALSPOUT
        idDirs = glob.glob(
            os.path.normpath(os.path.join(self.getWorkDirectory(), "*INDEX.log"))
        )
        if len(idDirs) == 0:
            print(
                "i2 Xia2/Dials : Failed to find Xia2/Dials logfiles (using pattern *INDEX.log)"
            )
            return
        searchOutputDirs = []
        for idDir in idDirs:
            _, fltmp = os.path.split(idDir)
            splitD = fltmp.split("_")
            if len(splitD) > 4:
                searchOutputDirs.append(
                    os.path.normpath(
                        os.path.join(
                            self.getWorkDirectory(), splitD[1], splitD[2], splitD[3]
                        )
                    )
                )
        SearchPoss = [
            [
                "index",
                "*_" + splitD[3] + "*strong.expt",
                "*strong.refl",
                "(spot finding)",
                "spotfinder",
            ],
            ["index", "*_indexed.expt", "*_indexed.refl", "(indexing)", "indexed"],
            ["refine", "*_refined.expt", "*_refined.refl", "(refinement)", "refined"],
            [
                "integrate",
                "*_integrated.expt",
                "*_integrated.refl",
                "(integration)",
                "integrated",
            ],
        ]
        # Loop over the spot-finding, indexing, refinement & integration folders looking for the most relevant files.
        for OutputDir in searchOutputDirs:
            for SetType in SearchPoss:
                JsonCands = glob.glob(
                    os.path.normpath(os.path.join(OutputDir, SetType[0], SetType[1]))
                )
                PickCands = glob.glob(
                    os.path.normpath(os.path.join(OutputDir, SetType[0], SetType[2]))
                )
                tmp1 = []
                tmp2 = []
                jnam = None
                pnam = None
                for JsonCan in JsonCands:
                    _, fltmp = os.path.split(JsonCan)
                    ctmp1 = re.match("[0-9]*", fltmp).group(0)
                    if ctmp1:
                        tmp1.append(int(ctmp1))
                    if not jnam:
                        jnam = re.sub("[0-9]*", "", fltmp, 1)
                for PickCan in PickCands:
                    _, fltmp = os.path.split(PickCan)
                    ctmp2 = re.match("[0-9]*", fltmp).group(0)
                    if ctmp2:
                        tmp2.append(int(ctmp2))
                    if not pnam:
                        pnam = re.sub("[0-9]*", "", fltmp, 1)
                matchlist = list(set(tmp1).intersection(tmp2))
                if matchlist:
                    selNumber = max(matchlist)
                    jfile_act = os.path.normpath(
                        os.path.join(OutputDir, SetType[0], str(selNumber) + jnam)
                    )
                    pfile_act = os.path.normpath(
                        os.path.join(OutputDir, SetType[0], str(selNumber) + pnam)
                    )
                    use_jnam = os.path.normpath(
                        os.path.join(self.getWorkDirectory(), SetType[4] + ".expt")
                    )
                    use_pnam = os.path.normpath(
                        os.path.join(self.getWorkDirectory(), SetType[4] + ".refl")
                    )
                    # ------- Copy the files over & log them.
                    shutil.copyfile(jfile_act, use_jnam)
                    shutil.copyfile(pfile_act, use_pnam)
                    outContJ.append(outContJ.makeItem())
                    outContJ[-1].fullPath = use_jnam
                    outContJ[-1].annotation = "Dials json file " + SetType[3]
                    outContP.append(outContP.makeItem())
                    outContP[-1].fullPath = use_pnam
                    outContP[-1].annotation = "Dials pick file " + SetType[3]

    def handleXia2DotTxtChanged(self, filename):
        # remove xia2.txt nodes
        for xia2TxtNode in self.xmlroot.xpath("Xia2Txt"):
            self.xmlroot.remove(xia2TxtNode)
        xia2TxtNode = etree.SubElement(self.xmlroot, "Xia2Txt")
        with open(filename, "r") as xia2DotTxtFile:
            xia2TxtNode.text = etree.CDATA(xia2DotTxtFile.read())
        self.flushXML()

    def flushXML(self):
        tmpFilename = self.makeFileName("PROGRAMXML") + "_tmp"
        with open(tmpFilename, "wb") as xmlFile:
            xmlFile.write(etree.tostring(self.xmlroot, pretty_print=True))
        if os.path.exists(self.makeFileName("PROGRAMXML")):
            os.remove(self.makeFileName("PROGRAMXML"))
        os.rename(tmpFilename, self.makeFileName("PROGRAMXML"))
