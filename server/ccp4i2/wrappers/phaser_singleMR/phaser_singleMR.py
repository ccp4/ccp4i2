import os
import xml.etree.ElementTree as ET

from ccp4i2.core import CCP4XtalData
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.smartie import smartie


class phaser_singleMR(CPluginScript):
    TASKNAME = 'phaser_singleMR'
    TASKCOMMAND = 'phaser'
    PERFORMANCECLASS = 'CExpPhasPerformance'
    ASYNCHRONOUS = True

    def __init__(self, *args, **kwargs):
        self.hklin1 = None
        self.bFData = None
        self.bIData = None
        self.xmlout = None
        CPluginScript.__init__(self, *args, **kwargs)

    def processInputFiles(self):
        cols1 = []
        self.container.inputData.F_SIGF.loadFile()
        self.container.inputData.F_SIGF.setContentFlag()
        self.bFData = self.container.inputData.F_SIGF.contentFlag == 2 or self.container.inputData.F_SIGF.contentFlag == 4
        self.bIData = self.container.inputData.F_SIGF.contentFlag == 1 or self.container.inputData.F_SIGF.contentFlag == 3
        if self.bIData:
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN])
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_IMEAN])
        if self.bFData:
            cols1.append(['F_SIGF', CCP4XtalData.CObsDataFile.CONTENT_FLAG_FMEAN])
        self.hklin1, _, error1 = self.makeHklin0(cols1)
        # Need to join up the previous i2-only mini-mtz's to get back the prior mtz file ...
        self.seqin = os.path.join(self.workDirectory, 'input_seq.fasta')
        self.container.inputData.ASUFILE.writeFasta(self.seqin)
        return error1

    def processOutputFiles(self):
        num_sol = 1
        for i in range(1, num_sol + 1):
            xyzout = os.path.join(self.getWorkDirectory(), "SingleMR." + str(i) + ".pdb")
            if os.path.exists(xyzout):
                self.container.outputData.XYZOUT.append(self.container.outputData.XYZOUT.makeItem())
                self.container.outputData.XYZOUT[-1].setFullPath(xyzout)
                self.container.outputData.XYZOUT[-1].annotation = 'Positioned coordinates for solution ' + str(i)
            else:
                self.appendErrorReport(201, xyzout)
                return CPluginScript.FAILED
            hklout = os.path.join(self.getWorkDirectory(), "SingleMR." + str(i) + ".mtz")
            if os.path.exists(hklout):
                self.container.outputData.HKLOUT.append(self.container.outputData.HKLOUT.makeItem())
                self.container.outputData.HKLOUT[-1].setFullPath(hklout)
            else:
                self.appendErrorReport(201, hklout)
                return CPluginScript.FAILED
        self.splitHkloutList(miniMtzsOut=['MAPOUT', 'ABCDOUT'],
                             programColumnNames=['FWT,PHWT', 'HLA,HLB,HLC,HLD'],
                             outputBaseName=['MAPOUT', 'ABCDOUT'],
                             outputContentFlags=[1, CCP4XtalData.CPhsDataFile.CONTENT_FLAG_HL],
                             infileList=self.container.outputData.HKLOUT)
        for indx in range(len(self.container.outputData.MAPOUT)):
            self.container.outputData.MAPOUT[indx].annotation = 'Map for solution ' + str(indx + 1)
            self.container.outputData.MAPOUT[indx].contentFlag = 1
            self.container.outputData.MAPOUT[indx].subType = 1
            self.container.outputData.ABCDOUT[indx].annotation = 'H-L Co-efficients' + str(indx + 1)
        self.parseLogfile()
        return CPluginScript.SUCCEEDED

    def makeCommandAndScript(self):
        # nb. the -xml flag doesn't work in phaser & the html flag does not work in smartie.
        self.appendCommandScript("TITLE SingleMR")
        self.appendCommandScript("MODE MR_ATOM")
        self.appendCommandScript("ROOT SingleMR")
        self.appendCommandScript("HKLIN \"%s\" "%str(self.hklin1))
        if self.bIData:
            self.appendCommandScript("LABIN I=F_SIGF_I SIGI=F_SIGF_SIGI")
        if self.bFData:
            self.appendCommandScript("LABIN F=F_SIGF_F SIGF=F_SIGF_SIGF")  # Should be I or F..
        if not self.bIData and not self.bFData:
            print("phaser_singleMR : Something has gone horribly wrong with the input") # Should really trigger an exception here...
        if self.container.inputData.RESOL_ON:
            self.appendCommandScript("RESOLUTION HIGH %s LOW %s"%(str(self.container.inputData.RESOL_HI),
                                                                  str(self.container.inputData.RESOL_LO)))
        # self.appendCommandScript("RESOLUTION")  # The resolution.
        self.appendCommandScript("SGALTERNATIVE SELECT NONE")  # Need to fix the space group so user can run it manually.
        # Define the search
        self.appendCommandScript("ENSEMBLE SINGLEATOM ATOM %s"%str(self.container.inputData.SINGLE_ATOM_TYPE)) # Needs the atom & number
        self.appendCommandScript("SEARCH ENSEMBLE SINGLEATOM NUMBER %s"%str(self.container.inputData.SINGLE_ATOM_NUM))
        self.appendCommandScript("LLGCOMPLETE SCATTERER %s"%str(self.container.inputData.LLG_COMPL_ATOMTYP)) # Needs the LLG atom
        if self.container.inputData.LLG_COMPL_ATMSEP_ON:
            self.appendCommandScript("LLGCOMPLETE CLASH %s"%str(self.container.inputData.LLG_COMPL_ATMSEP))
        if self.container.inputData.LLG_COMPL_SIGCO_ON:
            self.appendCommandScript("LLGCOMPLETE SIGMA %s"%str(self.container.inputData.LLG_COMPL_SIGCO))
        if self.container.inputData.LLG_COMPL_MAXCYC_ON:
            self.appendCommandScript("LLGCOMPLETE NCYC %s"%str(self.container.inputData.LLG_COMPL_MAXCYC))
        # Need the ASU Unit to be defined
        self.appendCommandScript("COMPOSITION BY ASU")   # Obviously can be something else.
        self.appendCommandScript("COMPOSITION PROTEIN SEQ \"%s\" NUMBER 1"%str(self.seqin))
        # Expert switches
        if self.container.inputDataExp.EXP_PACKCT_ON:
            if str(self.container.inputDataExp.EXP_PACKCT_TYPE) != "ALL":
                self.appendCommandScript("PACK SELECT %s CUTOFF %s"%(str(self.container.inputDataExp.EXP_PACKCT_TYPE),
                                                                     str(self.container.inputDataExp.EXP_PACKCT_AMT)))
                if self.container.inputDataExp.EXP_QKPACK_ON:
                    self.appendCommandScript("PACK QUICK ON")
                else:
                    self.appendCommandScript("PACK QUICK OFF")
            else:
                self.appendCommandScript("PACK SELECT ALL")
        if self.container.inputDataExp.EXP_TRAN_SRCHPK_ON:
            if str(self.container.inputDataExp.EXP_TRAN_SRCHPK_TYPE) != "ALL":
                self.appendCommandScript("PEAKS TRA SELECT %s CUTOFF %s"%(str(self.container.inputDataExp.EXP_TRAN_SRCHPK_TYPE),
                                                                          str(self.container.inputDataExp.EXP_TRAN_SRCHPK_AMT)))
            else:
                self.appendCommandScript("PEAKS TRA SELECT ALL")
        # PURGE TRA ENABLE ON  & then PURGE TRA PERCENT 75  & also   PURGE TRA NUMBER 52
        if self.container.inputDataExp.EXP_PURGE_TRANPK_ON:
            self.appendCommandScript("PURGE TRA ENABLE ON")
            self.appendCommandScript("PURGE TRA PERCENT %s"%str(self.container.inputDataExp.EXP_PURGE_TRANPK_PER))
            self.appendCommandScript("PURGE TRA NUMBER %s"%str(self.container.inputDataExp.EXP_PURGE_TRANPK_NUM))
        # RESOLUTION AUTO HIGH 0.95 LOW 11.0
        if self.container.inputDataExp.EXP_RESRAN_HRREFINE_ON:
            self.appendCommandScript("RESOLUTION AUTO HIGH %s LOW %s"%(str(self.container.inputDataExp.EXP_RESRAN_HRREFINE_HI),
                                                                       str(self.container.inputDataExp.EXP_RESRAN_HRREFINE_LO)))
        self.xmlout = self.makeFileName('PROGRAMXML')
        return CPluginScript.SUCCEEDED

    def parseLogfile(self):
        logfile = self.makeFileName("LOG")
        smin = smartie.parselog(logfile)
        xmltree, xrttree = _convert_log(smin)
        xmltree.write(str(self.workDirectory / "program.xml"))
        with (self.workDirectory / "program.xrt").open("w") as xrto:
            xrto.write((ET.tostring(xrttree.getroot()).decode()).replace("ns0", "xrt"))


def _convert_log(logfile):
    job_tag = "Job"
    job_path = "/Job"
    xrtns = "{http://www.ccp4.ac.uk/xrt}%s"

    xrttree = ET.ElementTree(ET.Element("report"))
    e0 = xrttree.getroot()
    e0 = ET.SubElement(e0, xrtns % "results")

    xmltree = ET.ElementTree(ET.Element(job_tag))
    f0 = xmltree.getroot()

    e1 = ET.SubElement(e0, xrtns % "title", select=job_path + "/Title")

    f1 = ET.SubElement(f0, "Title")
    f1.text = None

    for ip in range(logfile.nprograms()):
        program = logfile.program(ip)
        program_tag = f"SubJob_{ip:03d}"
        program_path = job_path + "/" + program_tag

        all_attributes = program.attributes()

        progname = ""
        if "name" in all_attributes:
            progname = program.get_attribute("name")

        elif "termination_name" in all_attributes:
            progname = program.get_attribute("termination_name")

        progtitle = "Run"
        if progname:
            progtitle += f" of {progname}"

        rundate = ""
        if "rundate" in all_attributes:
            rundate = program.get_attribute("rundate")
            progtitle += f" on {rundate}"

        runtime = ""
        if "runtime" in all_attributes:
            runtime = program.get_attribute("runtime")
            progtitle += f" at {runtime}"

        e1 = ET.SubElement(e0, xrtns % "section", title=progtitle)

        f1 = ET.SubElement(f0, program_tag)

        attributes = program.attributes()
        if "name" in attributes:
            e1.attrib["progname"] = program.get_attribute("name")

        if "rundate" in attributes:
            info_tag = "RunDate"
            e1.attrib["rundate"] = program_path + "/" + info_tag
            f2 = ET.SubElement(f1, info_tag)
            f2.text = program.get_attribute("rundate")

        if "runtime" in attributes:
            info_tag = "RunTime"
            e1.attrib["runtime"] = program_path + "/" + info_tag
            f2 = ET.SubElement(f1, info_tag)
            f2.text = program.get_attribute("runtime")

        e2 = None
        for ik in range(program.nkeytexts()):
            keytext = program.keytext(ik)
            keytext_tag = f"KeyText_{ik + 1:03d}"
            keytext_path = program_path + "/" + keytext_tag

            e2 = ET.SubElement(e1, xrtns % "keytext")
            e2.attrib["folded"] = "false"
            e2.attrib["name"] = keytext.name()
            e2.attrib["select"] = keytext_path

            f2 = ET.SubElement(f1, keytext_tag)
            f2.text = "\n" + keytext.message().strip() + "\n"

        for xrt_table_tag in ("graph", "table"):
            jt = 0
            e1a = e1

            if xrt_table_tag == "graph" and program.ntables() > 0:
                e1a = ET.SubElement(e1, xrtns % "graph")

            for it in range(program.ntables()):
                jt += 1
                table = program.table(it)
                table_tag = f"Table_{jt:03d}"
                table_path = program_path + "/" + table_tag

                e2 = ET.SubElement(e1a, xrtns % "table", select=table_path)
                e2.attrib["type"] = "plain"
                e2.attrib["title"] = table.title().strip()
                if xrt_table_tag == "table":
                    e2.attrib["folded"] = "true"

                for ic in range(table.ncolumns()):
                    column = table.table_column(ic)

                    e3 = ET.SubElement(e2, xrtns % "data", title=column.title())

                if xrt_table_tag == "graph":
                    f2 = ET.SubElement(f1, table_tag)
                    f2.text = table.data().rstrip() + "\n"

                    for ig in range(table.ngraphs()):
                        graph = table.table_graph(ig)
                        columns = graph.columns()

                        e3 = ET.SubElement(e2, xrtns % "plot")
                        e3.attrib["title"] = graph.title().strip()

                        for column in columns[1:]:
                            ET.SubElement(
                                e3,
                                "plotline",
                                xcol=str(columns[0]),
                                ycol=str(column),
                                colour="auto",
                            )

                        scaling = graph.scaling()
                        if scaling == "N":
                            e3.attrib["ymin"] = "0"

                        elif scaling and scaling != "A":
                            xrange, sep, yrange = scaling.partition("x")
                            if sep == "x":
                                xmin, sep, xmax = xrange.partition("|")
                                if sep == "|":
                                    e3.attrib["xmin"] = xmin.strip()
                                    e3.attrib["xmax"] = xmax.strip()

                                ymin, sep, ymax = yrange.partition("|")
                                if sep == "|":
                                    e3.attrib["ymin"] = ymin.strip()
                                    e3.attrib["ymax"] = ymax.strip()

    ET.SubElement(f0, "Files")

    return xmltree, xrttree
