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
import xml.etree.ElementTree as ET
from .utils import i2run


def test_1h1s_selfrot(mtz1h1s):
    """Test molrep self-rotation function with 1h1s data."""
    args = ["molrep_selfrot"]
    args += ["--F_SIGF", f"fullPath={mtz1h1s}", "columnLabels=/*/*/[FP,SIGFP]"]

    with i2run(args) as job:
        # Check PostScript output exists
        ps_files = list(job.glob("*.ps"))
        assert len(ps_files) >= 1, f"No PS output: {list(job.iterdir())}"

        # Check program.xml has Patterson peaks
        tree = ET.parse(job / "program.xml")
        root = tree.getroot()
        peaks = root.findall(".//Patterson/Peak")
        assert len(peaks) >= 10, f"Expected >= 10 Patterson peaks, got {len(peaks)}"

        # Verify peak structure has expected elements
        first_peak = peaks[0]
        for field in ("No", "Xfrac", "Yfrac", "Zfrac", "Dens", "Dens_sigma"):
            assert first_peak.find(field) is not None, f"Missing '{field}' in Patterson peak"
