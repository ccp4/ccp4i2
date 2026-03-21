import pytest
from .utils import i2run


@pytest.mark.skip(reason="Requires external PDB-REDO web service; not suitable for offline CI")
def test_pdb_redo():
    """Test PDB-REDO web services task.

    Submits a structure to the PDB-REDO server for re-refinement.
    Requires network access and may take several minutes.
    """
    pass
