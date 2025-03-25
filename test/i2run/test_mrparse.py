import gemmi
from .utils import demoData, i2run


def test_mrparse():
    args = ["mrparse"]
    args += ["--SEQIN", demoData("gamma", "gamma.pir")]
    args += ["--DATABASE", "PDB"]
    args += ["--USEAPI", "False"]
    with i2run(args) as directory:
        outputs = list(directory.glob("*.pdb"))
        assert len(outputs) > 0
        for output in outputs:
            gemmi.read_structure(str(output))
