"""
Thin test of the failure mode that crashed servalcat/refmac validation:
clipper reading a gemmi-written mmCIF into MMDB.

iris_validation reads models via clipper.MMDBfile, whose MMDB cannot parse some
gemmi-written mmCIF and aborts with clipper::Message_fatal -- a C++ terminate,
not a Python exception, so it sails through Iris's own try/except and kills the
whole (in-process) pipeline. validate_protein.clipper_can_read probes that read
in a subprocess so the abort is isolated and Iris can be skipped gracefully.

This exercises exactly that: a gemmi-written mmCIF that clipper aborts on (2acy is
a deterministic crasher; mmCIF readability by this MMDB is content-dependent, so a
known-crashing model is used rather than an arbitrary one), and asserts the probe
returns False WITHOUT taking the test process down, and True for a readable PDB.

Runs only where clipper is present (CCP4); skipped on a slim interpreter.
"""
import subprocess
import sys

import pytest
import gemmi

from ccp4i2 import I2_TOP
from ccp4i2.wrappers.validate_protein.script.validate_protein import clipper_can_read

pytest.importorskip("clipper", reason="clipper not present (slim interpreter)")

PDB = I2_TOP / "demo_data" / "hypf" / "2acy.pdb"


@pytest.fixture
def gemmi_mmcif(tmp_path):
    """2acy written as gemmi-style mmCIF -- a confirmed clipper/MMDB crasher."""
    out = tmp_path / "2acy_gemmi.cif"
    gemmi.read_structure(str(PDB)).make_mmcif_document().write_file(str(out))
    return str(out)


def test_clipper_aborts_on_gemmi_mmcif(gemmi_mmcif):
    """Document the bug: clipper's MMDB read of gemmi-written mmCIF aborts."""
    probe = ("import sys, clipper\n"
             "f = clipper.MMDBfile(); f.read_file(sys.argv[1])\n"
             "m = clipper.MiniMol(); f.import_minimol(m)\n")
    rc = subprocess.run([sys.executable, "-c", probe, gemmi_mmcif],
                        capture_output=True).returncode
    assert rc != 0, "expected clipper/MMDB to abort reading gemmi-written mmCIF"


def test_probe_isolates_the_abort(gemmi_mmcif):
    """The fix's probe returns False for the crasher WITHOUT aborting this
    process, and True for a readable PDB."""
    # If isolation failed, the abort would kill this test process (no assert reached).
    assert clipper_can_read(gemmi_mmcif) is False
    assert clipper_can_read(str(PDB)) is True


def test_probe_handles_missing_file():
    assert clipper_can_read("/no/such/file.cif") is False
