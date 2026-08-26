import pytest
from .utils import demoData, i2run


@pytest.mark.skip(reason="NOT IMPLEMENTED: no test body. Would need a prior refinement job in the same project; see the docstring.")
def test_adding_stats():
    """Test adding_stats_to_mmcif_i2 deposition preparation.

    This task requires XYZIN from a refinement job, ASUCONTENT,
    REFMACINPUTPARAMSXML, and reflection data.  It traces lineage
    back through the refinement → scaling job chain, making it hard
    to test in isolation without running a full pipeline first.
    """
    pass
