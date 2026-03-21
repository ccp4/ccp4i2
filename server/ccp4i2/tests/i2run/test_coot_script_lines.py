import pytest
from .utils import demoData, i2run


@pytest.mark.skip(reason="Requires coot --no-graphics binary; may not be available in CI")
def test_rename_chain():
    """Test scripted model building with coot."""
    args = ["coot_script_lines"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--SCRIPT", "rename_chain(0, 'A', 'B')"]
    with i2run(args) as job:
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"
