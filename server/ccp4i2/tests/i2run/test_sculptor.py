import os
import subprocess
from shutil import which
from tempfile import TemporaryDirectory

import gemmi
from pytest import fixture, mark

from .utils import demoData, i2run

CCP4 = os.environ.get("CCP4", "")
CLUSTALW2 = os.path.join(CCP4, "libexec", "clustalw2") if CCP4 else None


def _has_clustalw2():
    return CLUSTALW2 and os.path.isfile(CLUSTALW2)


def _extract_sequence(pdb_path):
    """Extract amino-acid sequence from the first chain of a PDB file."""
    st = gemmi.read_structure(pdb_path)
    chain = st[0][0]
    return "".join(
        gemmi.find_tabulated_residue(res.name).one_letter_code
        for res in chain
        if gemmi.find_tabulated_residue(res.name).is_amino_acid()
    )


@fixture(scope="module")
def alignment_file():
    """Create a ClustalW alignment from gamma_model.pdb and gamma.pir."""
    pdb_path = demoData("gamma", "gamma_model.pdb")
    pir_path = demoData("gamma", "gamma.pir")

    model_seq = _extract_sequence(pdb_path)

    # Read target sequence from PIR file (skip header, strip trailing *)
    with open(pir_path) as f:
        lines = f.readlines()
    target_seq = "".join(l.strip() for l in lines if not l.startswith(">") and l.strip())
    target_seq = target_seq.rstrip("*")

    with TemporaryDirectory() as tmpdir:
        # Write FASTA input for clustalw2
        fasta = os.path.join(tmpdir, "input.fasta")
        with open(fasta, "w") as f:
            f.write(">target\n%s\n>model\n%s\n" % (target_seq, model_seq))

        aln_path = os.path.join(tmpdir, "input.aln")
        subprocess.run(
            [CLUSTALW2, "-INFILE=" + fasta, "-OUTFILE=" + aln_path, "-OUTPUT=CLUSTAL"],
            check=True, capture_output=True,
        )
        yield aln_path


@mark.skipif(not _has_clustalw2(), reason="clustalw2 not available")
def test_gamma_sculptor(alignment_file):
    """Test sculptor model trimming with gamma alignment."""
    args = ["sculptor"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--ALIGNMENTORSEQUENCEIN", "ALIGNMENT"]
    args += ["--ALIGNIN", alignment_file]
    with i2run(args) as job:
        # sculptor outputs a list of PDB files; check at least one exists
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"


@fixture(scope="module")
def two_copy_model(tmp_path_factory):
    """gamma_model.pdb with a second copy, chain B, moved 40 A along x."""
    st = gemmi.read_structure(demoData("gamma", "gamma_model.pdb"))
    copy = st[0][0].clone()
    copy.name = "B"
    for residue in copy:
        for atom in residue:
            atom.pos = gemmi.Position(atom.pos.x + 40.0, atom.pos.y, atom.pos.z)
    st[0].add_chain(copy)
    st.setup_entities()
    path = tmp_path_factory.mktemp("models") / "two_copies.pdb"
    st.write_pdb(str(path))
    return str(path)


@mark.skipif(not _has_clustalw2(), reason="clustalw2 not available")
def test_sculptor_trims_only_the_selected_chain(alignment_file, two_copy_model):
    """The atom selection is applied: Sculptor was given the whole file
    whatever it said (MDM2 job 23 asked for chain A of four copies and got
    all four). program.xml says how many chains came out, and the identity."""
    import xml.etree.ElementTree as ET
    args = ["sculptor"]
    args += ["--XYZIN", f"fullPath={two_copy_model}", "selection/text=A/"]
    args += ["--ALIGNMENTORSEQUENCEIN", "ALIGNMENT"]
    args += ["--ALIGNIN", alignment_file]
    with i2run(args) as job:
        root = ET.parse(job / "program.xml").getroot()
        assert root.findtext("selection_applied") == "True"
        assert [o.findtext("chains") for o in root.findall("output")] == ["1"]
        identities = root.findall("identity")
        assert len(identities) == 1 and float(identities[0].text) > 0
