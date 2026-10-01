import re

import pytest
from .utils import demoData, i2run


@pytest.mark.slow
def test_gamma_lattice():
    """Test SIMBAD lattice search with gamma native data.

    Marked slow because SIMBAD lattice search queries a database.
    """
    args = ["SIMBAD"]
    args += ["--F_SIGF", demoData("gamma", "gamma_native.mtz")]
    args += ["--SIMBAD_SEARCH_LEVEL", "Lattice"]
    args += ["--SIMBAD_NPROC", "1"]
    with i2run(args) as job:
        assert (job / "diagnostic.xml").exists()
        # The gamma cell matches 1gyu best, and its placement refines.
        assert (job / "1gyu_refinement_output.pdb").exists()
        # Each output names the structure it came from.
        params = (job / "params.xml").read_text()
        assert "SIMBAD hit 1: 1gyu, placed and refined" in params
        assert not re.search(r"PDB #\d from REFMAC", params)
