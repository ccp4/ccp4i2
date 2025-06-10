import json
from pytest import approx
import gemmi
from .utils import demoData, i2run


def test_mrparse():
    args = ["mrparse"]
    args += ["--SEQIN", demoData("gamma", "gamma.pir")]
    args += ["--DATABASE", "PDB"]
    args += ["--USEAPI", "False"]
    with i2run(args) as job:
        outputs = list(job.glob("*.pdb"))
        assert len(outputs) > 0
        for output in outputs:
            gemmi.read_structure(str(output))
        with (job / "mrparse_0" / "homologs.json").open() as f:
            homologues = json.load(f)
            assert max(h["seq_ident"] for h in homologues) == approx(1.0)
