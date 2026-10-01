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
