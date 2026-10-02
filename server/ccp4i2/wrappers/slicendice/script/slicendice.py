import json
import re
import shutil
from pathlib import Path

from lxml import etree

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4XtalData import CObsDataFile


# SliceNDice's own test of a placement (slicendice.dice.submit_jobs_MX):
# after refinement, R and R-free both below 0.45.
SOLVED_BELOW = 0.45


def split_ranges(results):
    """{split: ["26-111", "247-275, 360-490"]}: each split's search models, as
    residue ranges, from slicendice_results.json's "slice" section."""
    out = {}
    for split, info in results.get("slice", {}).items():
        models = []
        for _, segments in sorted(info.get("residues_ranges", {}).items()):
            models.append(", ".join(
                "%s-%s" % (a.split(":")[-1], b.split(":")[-1]) for a, b in segments))
        out[split] = models
    return out


def phaser_components(log_text):
    """Each search model's placement in Phaser's best solution, in the order
    they were placed: [(cluster, TFZ, LLG after placing it, clashes)].

    SliceNDice keeps one TFZ and LLG per split, those of the last model
    placed; a split that placed one lobe clearly (TFZ 26) and the other not
    at all (TFZ 4.5, LLG +10) showed only the 4.5. Phaser's annotation
    history has them all: "SOLU SET RFZ=.. TFZ=.. PAK=.. LLG=.. RFZ=.. ..",
    one RFZ group per model, followed by a SOLU 6DIM line naming each."""
    blocks = re.findall(r"SOLU SET\s+(.*)\n((?:.*\n){0,40}?)\s*(?:Solution #|$)", log_text)
    if not blocks:
        return []
    def final_llg(history):
        values = re.findall(r"LLG=(-?[\d.]+)", history)
        return float(values[-1]) if values else float("-inf")
    history, rest = max(blocks, key=lambda b: final_llg(b[0]))
    clusters = re.findall(r"SOLU 6DIM ENSE \S*?cluster_(\d+)", rest)
    out = []
    for n, group in enumerate(g for g in re.split(r"(?=RFZ=)", history) if g.startswith("RFZ=")):
        tfz = re.search(r"(?<!=)TFZ=(-?[\d.]+)", group)
        llg = re.search(r"LLG=(-?[\d.]+)", group)
        pak = re.search(r"PAK=(\d+)", group)
        out.append((clusters[n] if n < len(clusters) else str(n),
                    tfz.group(1) if tfz else "", llg.group(1) if llg else "",
                    pak.group(1) if pak else ""))
    return out


def best_solution(results):
    """(split, solved): the split with the lowest R-free after refinement, and
    whether it meets SliceNDice's own test of a solution."""
    dice = results.get("dice", {})
    if not dice:
        return None, False
    split = min(dice, key=lambda k: float(dice[k]["final_r_free"]))
    best = dice[split]
    solved = (float(best["final_r_fact"]) < SOLVED_BELOW
              and float(best["final_r_free"]) < SOLVED_BELOW)
    return split, solved


class slicendice(CPluginScript):
    TASKNAME = "slicendice"
    TASKCOMMAND = "slicendice"
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
        # The best placement, and whether it is a solution at all. (The
        # lowest R-free was reported as "the best MR solution" whatever it
        # was: R-free 0.555, no solution, read as a result.)
        best_split, solved = best_solution(jdd)
        if best_split is None:
            self.appendErrorReport(19122)
            print("SlicenDice: NO solution found in the json outfile.")
            return CPluginScript.FAILED
        xyz = Path(jdd["dice"][best_split]["xyzout"]).resolve()
        hkl = Path(jdd["dice"][best_split]["hklout"]).resolve()
        xyzout = self.workDirectory / xyz.name
        hklout = self.workDirectory / hkl.name

        # setFullPath, not assignment: `out.XYZOUT = path` replaced the file
        # object with a bare Path, so nothing could annotate it and the
        # gleaner saw no output file.
        if xyz.is_file():
            shutil.copy2(xyz, xyzout)
            out.XYZOUT.setFullPath(str(xyzout))
        if hkl.is_file():
            shutil.copy2(hkl, hklout)
            out.HKLOUT.setFullPath(str(hklout))

        # Split out data objects that have been generated. Do this after applying the annotation, and flagging
        # above, since splitHklout needs to know the ABCDOUT contentFlag
        outputFiles = ["FPHIOUT", "DIFFPHIOUT"]
        outputColumns = ["FWT,PHWT", "DELFWT,PHDELWT"]
        errorReport = self.splitHklout(outputFiles, outputColumns, infile=hklout)
        if errorReport.maxSeverity() > CCP4ErrorHandling.SEVERITY_WARNING:
            return errorReport

        ranges = split_ranges(jdd)
        n_splits = best_split.split("_")[-1]
        rfree_text = "%.3f" % float(jdd["dice"][best_split]["final_r_free"])
        # A split can pass on R-free with one of its pieces never placed:
        # say so when any piece's TFZ is below 8 (Phaser's clear placement).
        best_log = jdd["dice"][best_split].get("phaser_logfile")
        best_parts = (phaser_components(Path(best_log).read_text(encoding="utf-8", errors="replace"))
                      if best_log and Path(best_log).is_file() else [])
        partial = solved and any(t and float(t) < 8 for _, t, _, _ in best_parts)
        what = ("SliceNDice partial solution" if partial else
                "SliceNDice solution" if solved else "SliceNDice, no solution")
        out.XYZOUT.annotation = "%s: %s split%s, R-free %s" % (
            what, n_splits, "" if n_splits == "1" else "s", rfree_text)
        out.HKLOUT.annotation = "%s: refined data and map coefficients" % what
        out.FPHIOUT.annotation = "%s: 2Fo-Fc map coefficients" % what
        out.DIFFPHIOUT.annotation = "%s: Fo-Fc map coefficients" % what

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
        etree.SubElement(xmlbcyc, "Solved").text = str(solved)
        etree.SubElement(xmlbcyc, "Partial").text = str(partial)
        for split, models in sorted(ranges.items()):
            xmlsplit = etree.SubElement(xmlRI, "Split", id=split.split("_")[-1])
            for model in models:
                etree.SubElement(xmlsplit, "Model").text = model
        # Get solns & save
        for key in jdd["dice"].keys():
            xmlcyc = etree.SubElement(xmlRI, "Sol")
            etree.SubElement(xmlcyc, "SolID").text = str(key.split("_")[-1])
            etree.SubElement(xmlcyc, "llg").text = str(jdd["dice"][key]["phaser_llg"])
            etree.SubElement(xmlcyc, "tfz").text = str(jdd["dice"][key]["phaser_tfz"])
            etree.SubElement(xmlcyc, "srf").text = str(jdd["dice"][key]["final_r_fact"])
            etree.SubElement(xmlcyc, "sre").text = str(jdd["dice"][key]["final_r_free"])
            log = jdd["dice"][key].get("phaser_logfile")
            if log and Path(log).is_file():
                for cluster, tfz, llg, pak in phaser_components(
                        Path(log).read_text(encoding="utf-8", errors="replace")):
                    etree.SubElement(xmlcyc, "Component", cluster=cluster, tfz=tfz,
                                     llg=llg, clashes=pak)
        # Save xml
        xmlfile = open(self.makeFileName("PROGRAMXML"), "wb")
        xmlString = etree.tostring(rootNode, pretty_print=True)
        xmlfile.write(xmlString)
        xmlfile.close()
        # A placement that is not a solution keeps its files (to look at) but
        # does not finish as a success.
        return CPluginScript.SUCCEEDED if solved else CPluginScript.UNSATISFACTORY
