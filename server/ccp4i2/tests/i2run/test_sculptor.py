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
import os
import subprocess
from shutil import which
from tempfile import TemporaryDirectory

import gemmi
from pytest import fixture, mark

from .utils import demoData, i2run

CCP4 = os.environ.get("CCP4", "")
CLUSTALW2 = os.path.join(CCP4, "libexec", "clustalw2") if CCP4 else None


def _has_clustalw2():
    return CLUSTALW2 and os.path.isfile(CLUSTALW2)


def _extract_sequence(pdb_path):
    """Extract amino-acid sequence from the first chain of a PDB file."""
    st = gemmi.read_structure(pdb_path)
    chain = st[0][0]
    return "".join(
        gemmi.find_tabulated_residue(res.name).one_letter_code
        for res in chain
        if gemmi.find_tabulated_residue(res.name).is_amino_acid()
    )


@fixture(scope="module")
def alignment_file():
    """Create a ClustalW alignment from gamma_model.pdb and gamma.pir."""
    pdb_path = demoData("gamma", "gamma_model.pdb")
    pir_path = demoData("gamma", "gamma.pir")

    model_seq = _extract_sequence(pdb_path)

    # Read target sequence from PIR file (skip header, strip trailing *)
    with open(pir_path) as f:
        lines = f.readlines()
    target_seq = "".join(l.strip() for l in lines if not l.startswith(">") and l.strip())
    target_seq = target_seq.rstrip("*")

    with TemporaryDirectory() as tmpdir:
        # Write FASTA input for clustalw2
        fasta = os.path.join(tmpdir, "input.fasta")
        with open(fasta, "w") as f:
            f.write(">target\n%s\n>model\n%s\n" % (target_seq, model_seq))

        aln_path = os.path.join(tmpdir, "input.aln")
        subprocess.run(
            [CLUSTALW2, "-INFILE=" + fasta, "-OUTFILE=" + aln_path, "-OUTPUT=CLUSTAL"],
            check=True, capture_output=True,
        )
        yield aln_path


@mark.skipif(not _has_clustalw2(), reason="clustalw2 not available")
def test_gamma_sculptor(alignment_file):
    """Test sculptor model trimming with gamma alignment."""
    args = ["sculptor"]
    args += ["--XYZIN", demoData("gamma", "gamma_model.pdb")]
    args += ["--ALIGNMENTORSEQUENCEIN", "ALIGNMENT"]
    args += ["--ALIGNIN", alignment_file]
    with i2run(args) as job:
        # sculptor outputs a list of PDB files; check at least one exists
        pdb_files = list(job.glob("*.pdb"))
        assert len(pdb_files) > 0, f"No PDB output: {list(job.iterdir())}"
