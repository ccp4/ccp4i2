import json
import shutil
from pathlib import Path

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CObsDataFile


class slicendice(CPluginScript):
    TASKTITLE = "SliceNDice"
    TASKNAME = "slicendice"
    TASKCOMMAND = "slicendice"
    MAINTAINER = "ronan.keegan@stfc.ac.uk"
    PERFORMANCECLASS = "CRefinementPerformance"
    ERROR_CODES = {
        19121: {
            "description": "SliceNDice, Json Data file not found. "
            "Please check the SliceNDice log file for details."
        },
        19122: {
            "description": "SliceNDice, No solution found in json file. "
            "Please check the SliceNDice log file for details."
        },
    }

    def processInputFiles(self):
        dataObjects = [["F_SIGF", CObsDataFile.CONTENT_FLAG_FMEAN], "FREERFLAG"]
        self.hklin, errorReport = self.makeHklin(dataObjects)
        return errorReport

    def makeCommandAndScript(self):
        inp = self.container.inputData
        par = self.container.controlParameters
        mod = self.container.modelParameters
        self.appendCommandLine(["--xyzin", inp.XYZIN])
        self.appendCommandLine(["--hklin", self.hklin])
        seqFile = self.workDirectory / "SEQIN.fasta"
        inp.ASUIN.writeFasta(fileName=str(seqFile))
        self.appendCommandLine(["--seqin", seqFile])
        self.appendCommandLine(["--bfactor_column", mod.BFACTOR_TREATMENT])
        if mod.BFACTOR_TREATMENT == "plddt":
            self.appendCommandLine(["--plddt_threshold", mod.PLDDT_THRESHOLD])
        elif mod.BFACTOR_TREATMENT == "rms":
            self.appendCommandLine(["--rms_threshold", mod.RMS_THRESHOLD])
        self.appendCommandLine(["--min_splits", mod.MIN_SPLITS])
        self.appendCommandLine(["--max_splits", mod.MAX_SPLITS])
        self.appendCommandLine(["--nproc", par.NPROC])
        self.appendCommandLine(["--ncyc_refmac", par.NCYC])
        self.appendCommandLine(["--no_mols", par.NO_MOLS])

    def processOutputFiles(self):
        out = self.container.outputData
        # Load Json
        try:
            jsfloc = self.workDirectory / "slicendice_0" / "slicendice_results.json"
            with jsfloc.open() as jfi:
                jdd = json.load(jfi)
        except:
            # Failed to find a solution in the json file.
            self.appendErrorReport(19121)
            print("SlicenDice: NO json output found.")
            return CPluginScript.FAILED
        # Get the best soln from the json file
        try:
            lowest_rfree = 1.0
            best_split = None
            for split in jdd["dice"].keys():
                if jdd["dice"][split]["final_r_free"] <= lowest_rfree:
                    lowest_rfree = jdd["dice"][split]["final_r_free"]
                    best_split = split
        except:
            # Failed to find a solution in the json file.
            self.appendErrorReport(19122)
            print("SlicenDice: NO solution found in the json outfile.")
            return CPluginScript.FAILED
        xyz = Path(jdd["dice"][best_split]["xyzout"]).resolve()
        hkl = Path(jdd["dice"][best_split]["hklout"]).resolve()
        xyzout = self.workDirectory / xyz.name
        hklout = self.workDirectory / hkl.name

        if xyz.is_file():
            shutil.copy2(xyz, xyzout)
            out.XYZOUT = xyzout
        if hkl.is_file():
            shutil.copy2(hkl, hklout)
            out.HKLOUT = hklout

        # Split out data objects that have been generated. Do this after applying the annotation, and flagging
        # above, since splitHklout needs to know the ABCDOUT contentFlag
        outputFiles = ["FPHIOUT", "DIFFPHIOUT"]
        outputColumns = ["FWT,PHWT", "DELFWT,PHDELWT"]
        errorReport = self.splitHklout(outputFiles, outputColumns, infile=hklout)
        if errorReport.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return errorReport

        # Set performance indicators
        bid = str(best_split.split("_")[-1])
        rwork = str(jdd["dice"][best_split]["final_r_fact"])
        rfree = str(jdd["dice"][best_split]["final_r_free"])
        out.PERFORMANCEINDICATOR.RFactor = rwork
        out.PERFORMANCEINDICATOR.RFree = rfree

        # xml info
        rootNode = etree.Element("SliceNDice")
        xmlRI = etree.SubElement(rootNode, "RunInfo")
        xmlbcyc = etree.SubElement(xmlRI, "Best")
        etree.SubElement(xmlbcyc, "bid").text = bid
        etree.SubElement(xmlbcyc, "R").text = rwork
        etree.SubElement(xmlbcyc, "RFree").text = rfree
        # Get solns & save
        for key in jdd["dice"].keys():
            xmlcyc = etree.SubElement(xmlRI, "Sol")
            etree.SubElement(xmlcyc, "SolID").text = str(key.split("_")[-1])
            etree.SubElement(xmlcyc, "llg").text = str(jdd["dice"][key]["phaser_llg"])
            etree.SubElement(xmlcyc, "tfz").text = str(jdd["dice"][key]["phaser_tfz"])
            etree.SubElement(xmlcyc, "srf").text = str(jdd["dice"][key]["final_r_fact"])
            etree.SubElement(xmlcyc, "sre").text = str(jdd["dice"][key]["final_r_free"])
        # Save xml
        xmlfile = open(self.makeFileName("PROGRAMXML"), "wb")
        xmlString = etree.tostring(rootNode, pretty_print=True)
        xmlfile.write(xmlString)
        xmlfile.close()
