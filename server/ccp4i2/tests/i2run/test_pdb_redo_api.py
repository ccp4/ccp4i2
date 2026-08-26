import pytest
from .utils import i2run


@pytest.mark.skip(reason="NOT IMPLEMENTED: no test body. Would also need the live PDB-REDO service, so it could not run offline; see the docstring.")
def test_pdb_redo():
    """Test PDB-REDO web services task.

    Submits a structure to the PDB-REDO server for re-refinement.
    Requires network access and may take several minutes.
    """
    pass
