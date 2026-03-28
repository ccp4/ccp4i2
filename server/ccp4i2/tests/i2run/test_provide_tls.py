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
from .utils import demoData, i2run


def test_auto_from_coords():
    """Test ProvideTLS generating TLS groups from coordinates."""
    args = ["ProvideTLS"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    with i2run(args) as job:
        tlsout = job / "TLSFILE.tls"
        assert tlsout.exists(), f"No TLS file: {list(job.iterdir())}"
        content = tlsout.read_text()
        assert "TLS" in content, "TLS file has no TLS keyword"
