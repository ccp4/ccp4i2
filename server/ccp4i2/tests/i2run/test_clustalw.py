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


def test_cdk2_cdk6_alignment(fasta_cdk2, fasta_cdk6):
    """Test clustalw pairwise alignment of CDK2 and CDK6."""
    args = ["clustalw"]
    args += ["--SEQUENCELISTORALIGNMENT", "SEQUENCELIST"]
    args += ["--SEQIN", fasta_cdk2]
    args += ["--SEQIN", fasta_cdk6]

    with i2run(args) as job:
        # Check alignment output exists
        aln_files = list(job.glob("*.aln"))
        assert len(aln_files) >= 1, f"No alignment output: {list(job.iterdir())}"

        # Check program XML has pairwise scores
        tree = ET.parse(job / "program.xml")
        root = tree.getroot()
        scores = root.findall(".//PairwiseScore")
        assert len(scores) >= 1, "No pairwise scores in program.xml"
