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


_GAMMA_FASTA = (
    ">gamma_chymotrypsinogen\n"
    "MKYLLPTAAAGLLLLAAQPAMAMTQKIIMKDGLPSDSKLIQFEWNNPDQ"
    "FQKDRISCGQVSIPNNGDCIYGDIAFKAGAWEALFNAAGATFESAERRQ"
    "ADKEGCYREATTCLPLFTIFCNAFMPQIHGAVEGYHWDKYGGMCDDRYC"
    "GTIMGGMYAGGCQAVKSAAADATFNIKNIFESRKIDGLMSQMMCGAAYPF"
    "YVKDGENRCCQAANFYKSGAVILTHPGSEPNLTYWKF\n"
)


def test_freetext_sequence():
    """Test ProvideSequence with pasted sequence text."""
    args = ["ProvideSequence"]
    args += ["--SEQUENCETEXT", _GAMMA_FASTA]
    with i2run(args) as job:
        asu = job / "CASUCONTENTOUT.asu.xml"
        assert asu.exists(), f"No ASU content file: {list(job.iterdir())}"
        tree = ET.parse(asu)
        seqs = tree.findall(".//sequence")
        assert len(seqs) >= 1, "No sequences found in ASU output"
