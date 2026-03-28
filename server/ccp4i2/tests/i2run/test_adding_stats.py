# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
import pytest
from .utils import demoData, i2run


@pytest.mark.skip(reason="Requires a prior refinement job in the same project (complex setup)")
def test_adding_stats():
    """Test adding_stats_to_mmcif_i2 deposition preparation.

    This task requires XYZIN from a refinement job, ASUCONTENT,
    REFMACINPUTPARAMSXML, and reflection data.  It traces lineage
    back through the refinement → scaling job chain, making it hard
    to test in isolation without running a full pipeline first.
    """
    pass
