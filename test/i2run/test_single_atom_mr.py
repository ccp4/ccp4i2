from gemmi import read_mtz_file, read_pdb
from pytest import fixture
from .urls import pdbe_fasta, redo_mtz
from .utils import download, i2run


@fixture(scope="module", name="fasta")
def fasta_fixture():
    with download(pdbe_fasta("3njw")) as path:
        yield path


@fixture(scope="module", name="mtz")
def mtz_fixture():
    with download(redo_mtz("3njw")) as path:
        yield path


def test_single_atom_mr(fasta, mtz):
    args = ["phaser_singleMR"]
    args += ["--F_SIGF", f"fullPath={mtz}", "columnLabels=/*/*/[FP,SIGFP]"]
    args += ["--FREERFLAG", f"fullPath={mtz}", "columnLabels=/*/*/[FREE]"]
    args += ["--ASUFILE", f"seqFile={fasta}"]
    args += ["--SINGLE_ATOM_NUM", "2"]
    with i2run(args) as job:
        structure = read_pdb(str(job / "SingleMR.1.pdb"))
        atoms = list(structure[0].all())
        assert len(atoms) > 0, "No atoms in the output structure"
        for name in ["ABCDOUT_1", "MAPOUT_1", "SingleMR.1"]:
            read_mtz_file(str(job / f"{name}.mtz"))
