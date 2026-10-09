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


def space_group(path):
    """The space group (gemmi) of an MTZ or coordinate file, or None."""
    try:
        import gemmi
        path = str(path)
        if path.lower().endswith(".mtz"):
            return gemmi.read_mtz_file(path).spacegroup
        return gemmi.find_spacegroup_by_name(gemmi.read_structure(path).spacegroup_hm)
    except Exception:
        return None


def relabel_mtz(source, destination, group):
    """Copy an MTZ with its space group set to ``group`` (a gemmi
    SpaceGroup), indices untouched: what Phaser's choice of another group
    of the point group means for the data, since it never reindexes."""
    import gemmi
    mtz = gemmi.read_mtz_file(str(source))
    mtz.spacegroup = group
    mtz.write_to_file(str(destination))


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
        # Phaser tests the point group's space groups unless told not to:
        # pinned to the data's group (SliceNDice's default, none), the
        # BAD1330 lobes reached TFZ 8.5 and 7.5 with 28 clashes in P 21 2 21,
        # POINTLESS's choice at confidence 0.23, and 15.8 and 17.6 in P 21 21 21.
        self.appendCommandLine(["--sgalternative", mod.SGALTERNATIVE if mod.SGALTERNATIVE.isSet() else "all"])
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
        # Phaser may solve in another space group of the point group
        # (SGALTERNATIVE all). SliceNDice 0.1.3 then reindexes the data for
        # REFMAC with an operator that permutes the axes while Phaser's
        # model stays in the data's setting (BAD1330 merged as P 21 2 21,
        # solved in P 21 21 21: refinement of mismatched frames, R-free
        # 0.562, and an XYZOUT whose CRYST1 no longer matches its
        # coordinates). Phaser itself never reindexes: its model and maps
        # are in the data's setting with the new label, so those are the
        # outputs, with the project's reflections relabelled to match.
        log = jdd["dice"][best_split].get("phaser_logfile")
        phaser_xyz = Path(log).with_suffix(".pdb") if log else None
        phaser_mtz = Path(log).parent / "phaser_mr_output.1.mtz" if log else None
        sg_in = space_group(self.hklin) if getattr(self, "hklin", None) else None
        sg_out = space_group(phaser_xyz) if phaser_xyz and phaser_xyz.is_file() else None
        changed = bool(sg_in and sg_out and sg_in.number != sg_out.number)
        if changed and phaser_xyz.is_file() and phaser_mtz.is_file():
            xyz, hkl = phaser_xyz, phaser_mtz
            inp = self.container.inputData
            for item, output in ((inp.F_SIGF, out.F_SIGF_OUT), (inp.FREERFLAG, out.FREERFLAG_OUT)):
                if item.isSet() and Path(str(item.fullPath)).is_file():
                    target = self.workDirectory / ("%s_%s.mtz" % (output.objectName(), sg_out.short_name()))
                    relabel_mtz(str(item.fullPath), target, sg_out)
                    output.setFullPath(str(target))
                    output.annotation = "%s relabelled %s, Phaser's choice (indices unchanged)" % (
                        item.objectName(), sg_out.hm)
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
        outputColumns = ["FWT,PHWT", "DELFWT,PHDELWT"]  # REFMAC's, or Phaser's after a group change
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
        # Each piece judged by its own search TFZ (Phaser's clear placement,
        # 8), whatever the refinement said: BAD1330's two lobes placed at
        # 15.8 and 17.6 and were the structure (built to R-free 0.253), yet
        # R-free after ten cycles was 0.509 and the 0.45 rule called it no
        # solution; Lck passed the rule with its N-lobe at 6.0, 19 clashes.
        tfzs = [float(t) for _, t, _, _ in best_parts if t]
        placed = bool(tfzs) and all(t >= 8 for t in tfzs)
        partial = bool(tfzs) and any(t < 8 for t in tfzs) and any(t >= 8 for t in tfzs)
        what = ("SliceNDice partial solution" if solved and partial else
                "SliceNDice solution" if solved else
                "SliceNDice placement, not yet a solution" if placed else
                "SliceNDice partial placement" if partial else
                "SliceNDice, no solution")
        pieces = ", ".join("%.1f" % t for t in tfzs)
        if changed:
            # The refinement was of mismatched frames: its R-free says nothing
            solved = False
            what = ("SliceNDice placement in %s, not the data's group" % sg_out.hm if placed else
                    "SliceNDice partial placement in %s" % sg_out.hm if partial else
                    "SliceNDice, no solution (searched %s)" % sg_out.hm)
            out.XYZOUT.annotation = "%s: %s split%s, %s; Phaser's model, refine against F_SIGF_OUT" % (
                what, n_splits, "" if n_splits == "1" else "s",
                "piece TFZ %s" % pieces if pieces else "no piece placed")
        else:
            out.XYZOUT.annotation = "%s: %s split%s, %s, R-free %s" % (
                what, n_splits, "" if n_splits == "1" else "s",
                "piece TFZ %s" % pieces if pieces else "no piece placed", rfree_text)
        out.HKLOUT.annotation = "%s: %s" % (
            what, "Phaser's data and map coefficients" if changed else "refined data and map coefficients")
        out.FPHIOUT.annotation = "%s: 2Fo-Fc map coefficients" % what
        out.DIFFPHIOUT.annotation = "%s: Fo-Fc map coefficients" % what

        # Set performance indicators
        bid = str(best_split.split("_")[-1])
        rwork = str(jdd["dice"][best_split]["final_r_fact"])
        rfree = str(jdd["dice"][best_split]["final_r_free"])
        if not changed:  # after a group change the refinement numbers mean nothing
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
        etree.SubElement(xmlbcyc, "Placed").text = str(placed)
        etree.SubElement(xmlbcyc, "Partial").text = str(partial)
        if sg_out is not None:
            etree.SubElement(xmlbcyc, "SpaceGroup").text = sg_out.hm
        if sg_in is not None:
            etree.SubElement(xmlbcyc, "SpaceGroupInput").text = sg_in.hm
        if sg_in is not None and sg_out is not None:
            etree.SubElement(xmlbcyc, "SpaceGroupChanged").text = str(changed)
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
        # Every piece placed by Phaser's cutoff is a success to build from,
        # whether or not ten cycles took R-free below 0.45; a placement with
        # a piece unplaced, or none, keeps its files (to look at) but does
        # not finish as a success.
        return CPluginScript.SUCCEEDED if (solved or placed) else CPluginScript.UNSATISFACTORY
