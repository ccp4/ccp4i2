"""SliceNDice end to end: the AlphaFold model of the Src-family kinase Lck,
trimmed to its kinase domain, against the Lck kinase domain crystal (4c3f,
1.72 A), split into its two lobes.

There was no test, and the task could not run at all: its def.xml required
an ENSEMBLES list nothing filled, so every job failed validity. Then its
outputs were assigned as bare paths (``out.XYZOUT = path``), which left
nothing to annotate or glean. About ten minutes.

It is a partial solution, and the test says so: the C-lobe places (TFZ 26.4,
0.46 A from 4c3f after csymmatch), the N-lobe does not (TFZ 6.0, 18 A off),
yet the refined R-free (0.427) passes SliceNDice's own test.
"""
import json
import xml.etree.ElementTree as ET

import gemmi

from .urls import pdbe_fasta, redo_mtz
from .utils import download, i2run


def test_lck_kinase_two_lobes(tmp_path):
    with download("https://alphafold.ebi.ac.uk/api/prediction/P06239") as info:
        with open(info, encoding="utf-8") as f:
            entry = json.load(f)[0]
    with download(entry["pdbUrl"]) as full, download(redo_mtz("4c3f")) as mtz, \
            download(pdbe_fasta("4c3f")) as fasta:
        structure = gemmi.read_structure(str(full))
        for chain in structure[0]:
            for i in reversed(range(len(chain))):
                if not 225 <= chain[i].seqid.num <= 509:
                    del chain[i]
        model = tmp_path / "lck_kinase.pdb"
        structure.write_pdb(str(model))

        args = ["slicendice"]
        args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
        args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FREE]"]
        args += ["--ASUIN", f"seqFile={fasta}"]
        args += ["--XYZIN", str(model)]
        args += ["--BFACTOR_TREATMENT", "plddt", "--SEARCH_PDB", "False", "--SEARCH_AFDB", "False"]
        # One number of splits: SliceNDice 0.1.3 tries only one of a range.
        args += ["--NO_MOLS", "1", "--MIN_SPLITS", "2", "--MAX_SPLITS", "2"]
        with i2run(args) as job:
            best = ET.parse(job / "program.xml").find(".//RunInfo/Best")
            assert best.findtext("Solved") == "True", ET.tostring(best)
            assert float(best.findtext("RFree")) < 0.45
            from ccp4i2.db import models
            record = models.Job.objects.filter(number=job.name.replace("job_", "")).first()
            xyz = models.File.objects.filter(job=record, job_param_name="XYZOUT").first()
            # The C-lobe places (TFZ 26), the N-lobe does not (TFZ about 6):
            # R-free passes SliceNDice's test, but it is a partial solution.
            assert xyz is not None and xyz.annotation.startswith("SliceNDice partial solution: 2 splits"), \
                xyz and xyz.annotation
            parts = ET.parse(job / "program.xml").findall(".//Sol/Component")
            assert float(parts[0].get("tfz")) > 8 > float(parts[1].get("tfz")), \
                [p.attrib for p in parts]
