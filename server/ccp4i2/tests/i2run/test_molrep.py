import re
import xml.etree.ElementTree as ET
import gemmi
from .utils import demoData, i2run


# TODO: Test long ligand names (e.g. 8xfm)


def test_molrep():
    args = ["molrep_pipe"]
    args += ["--inputData.F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--inputData.FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        for name in ["XYZOUT_MOLREP", "XYZOUT_SHEETBEND", "XYZOUT"]:
            gemmi.read_pdb(str(job / f"{name}.pdb"))
        xml = ET.parse(job / "program.xml")
        rworks = [float(e.text) for e in xml.findall(".//Cycle/r_factor")]
        rfrees = [float(e.text) for e in xml.findall(".//Cycle/r_free")]
        assert rworks[-1] < rworks[0]
        assert rfrees[-1] < rfrees[0]
        assert rworks[-1] < 0.26
        assert rfrees[-1] < 0.28
        # MOLREP's own result, not the placeholders once written in its place
        # (Gamma gives score 0.79, z-score 18.9).
        assert float(xml.findtext(".//MR_TF/mr_score")) > 0.5
        assert float(xml.findtext(".//MR_TF/mr_zscore")) > 10


def test_molrep_uses_the_selected_sequence(tmp_path):
    """One sequence selected from a two-sequence AU file reaches MOLREP.

    MolRep takes exactly one sequence (selectionMode 1). The selection was
    lost between the job's params and the MOLREP sub-job, so MOLREP got both
    and failed with "more than one sequence" (issue #669).
    """
    gamma = open(demoData("gamma", "gamma.asu.xml")).read()
    seq = re.search(r"\s*<CAsuContentSeq>.*?</CAsuContentSeq>", gamma, re.S).group(0)
    decoy = seq.replace("<name>Gamma</name>", "<name>Decoy</name>")
    asu = tmp_path / "two.asu.xml"
    asu.write_text(gamma.replace(seq, decoy + seq))

    args = ["molrep_pipe"]
    args += ["--inputData.F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--inputData.FREERFLAG", demoData("gamma", "freeR.mtz")]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--ASUIN", f"fullPath={asu}", "selection/Gamma=True", "selection/Decoy=False"]
    with i2run(args) as job:
        fastas = list(job.glob("job_*/SEQIN.fasta"))
        assert fastas, "MOLREP was given no sequence"
        for fasta in fastas:
            names = [line[1:].strip() for line in fasta.read_text().splitlines()
                     if line.startswith(">")]
            assert names == ["Gamma"]
        gemmi.read_pdb(str(job / "XYZOUT.pdb"))


def test_molrep_takes_an_mmcif_model(tmp_path):
    """An mmCIF search model reaches MOLREP as PDB.

    The wrapper converted it to PDB for the run, then wrote the original to
    molrep_input.pdb unconverted; MOLREP refused mmCIF under a .pdb name
    ("CIF_file must has extension .crd or .cif"), so every mmCIF model failed.
    """
    model = tmp_path / "gamma_model.cif"
    gemmi.read_structure(demoData("gamma", "gamma_model.pdb")).make_mmcif_document().write_file(str(model))
    args = ["molrep_mr"]
    args += ["--F_SIGF", demoData("gamma", "merged_intensities_Xe.mtz")]
    args += ["--XYZIN", str(model)]
    with i2run(args) as job:
        assert (job / "molrep_input.pdb").read_text().startswith(("HEADER", "CRYST1", "ATOM", "REMARK"))
        gemmi.read_pdb(str(job / "XYZOUT.pdb"))


def test_molrep_finds_the_space_group_and_uses_it(tmp_path):
    """Data labelled with the wrong space group of the right Laue group: the
    first MOLREP tries them all, the pipeline reindexes into the best, and the
    second MOLREP and the refinement must use the reindexed data (#785).
    Copying the inputs and then setting the reindexed file's path kept the
    original file's database id, which wins over the path, so the second MOLREP
    read the unreindexed data in the wrong group."""
    relabelled = {}
    for name, source in (("F_SIGF", "merged_intensities_Xe.mtz"), ("FREERFLAG", "freeR.mtz")):
        mtz = gemmi.read_mtz_file(demoData("gamma", source))
        mtz.spacegroup = gemmi.SpaceGroup("P 2 2 2")  # gamma is P 21 21 21
        path = tmp_path / f"p222_{source}"
        mtz.write_to_file(str(path))
        relabelled[name] = str(path)
    args = ["molrep_pipe"]
    args += ["--inputData.F_SIGF", relabelled["F_SIGF"]]
    args += ["--inputData.FREERFLAG", relabelled["FREERFLAG"]]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--SG_OPTIONS", "laue"]
    with i2run(args) as job:
        xml = ET.parse(job / "program.xml")
        chosen = xml.findtext(".//MolrepSpaceGroup/laue_group_alternatives/selected") or ""
        assert gemmi.SpaceGroup(chosen).number == 19, f"chose {chosen!r}"
        scores = [float(e.text) for e in xml.findall(".//MR_TF/mr_score")]
        assert scores and scores[-1] > 0.5, scores
        rfrees = [float(e.text) for e in xml.findall(".//Cycle/r_free")]
        assert rfrees and rfrees[-1] < 0.3, rfrees
