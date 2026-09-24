"""Parity test: the pure-Python Scalepack reader vs the CCP4 scalepack2mtz binary.

reflection_formats.read_scalepack reproduces scalepack2mtz so a merged .sca can
be imported on a CCP4-free (slim) server, where scalepack2mtz is not available.
scalepack2mtz is the reference: for a real demo .sca, the reader must produce
the same hkl set, the same I(+)/SIGI(+)/I(-)/SIGI(-) exactly, the same IMEAN/
SIGIMEAN (to float32 precision — those are computed then rounded), and the same
missing-mate pattern. Runs only where scalepack2mtz is on PATH (under CCP4);
skipped on a slim interpreter.
"""
import re
import shutil
import subprocess
import tempfile
from pathlib import Path

import pytest
import gemmi
import numpy as np

from ccp4i2 import I2_TOP
from ccp4i2.lib.utils.files.reflection_formats import read_scalepack

SCALEPACK2MTZ = shutil.which("scalepack2mtz")
pytestmark = pytest.mark.skipif(
    SCALEPACK2MTZ is None, reason="scalepack2mtz binary not on PATH"
)

SCA = (
    I2_TOP / "demo_data" / "baz2b" / "BAZ2BA_x828.xia2" / "3daii-run"
    / "DataFiles" / "nt5073v16_xBAZ2BAx8281_scaled.sca"
)


def _header_cell_sg(sca_path: Path):
    """Cell (6 floats) + space group from the .sca 3-line header (line 3)."""
    line3 = sca_path.read_text(errors="replace").splitlines()[2]
    cell = tuple(float(x) for x in re.findall(r"[-\d.]+", line3)[:6])
    sg_name = line3[60:].strip() or line3.split(str(int(cell[5])))[-1].strip()
    # the SG symbol trails the six cell numbers; grab everything after them
    m = re.match(r"\s*(?:[-\d.]+\s+){6}(.+)", line3)
    sg_name = m.group(1).strip() if m else sg_name
    return cell, gemmi.SpaceGroup(sg_name).number


@pytest.fixture(scope="module")
def reference(tmp_path_factory):
    cell, sgnum = _header_cell_sg(SCA)
    out = tmp_path_factory.mktemp("scalepack") / "ref.mtz"
    script = (
        f"CELL {cell[0]} {cell[1]} {cell[2]} {cell[3]} {cell[4]} {cell[5]}\n"
        f"SYMMETRY {sgnum}\n"
        "NAME PROJECT p CRYSTAL c DATASET d\n"
        "ANOMALOUS YES\n"
        "END\n"
    )
    subprocess.run(
        [SCALEPACK2MTZ, "HKLIN", str(SCA), "HKLOUT", str(out)],
        input=script, text=True, capture_output=True, check=True,
    )
    return gemmi.read_mtz_file(str(out)), cell, sgnum


def _by_hkl(mtz):
    arr = np.array(mtz, copy=False)
    labels = [c.label for c in mtz.columns]
    return {
        (int(r[0]), int(r[1]), int(r[2])): dict(zip(labels, r)) for r in arr
    }


def test_scalepack_reader_matches_binary(reference):
    ref_mtz, cell, sgnum = reference
    mine = read_scalepack(SCA, cell, sgnum, anomalous=True)

    ref = _by_hkl(ref_mtz)
    got = _by_hkl(mine)

    # Same set of reflections.
    assert set(got) == set(ref), (
        f"hkl set differs: only-mine={len(set(got)-set(ref))} "
        f"only-ref={len(set(ref)-set(got))}"
    )

    # Raw anomalous columns must be exact; the derived mean to float32 ULP.
    exact = ["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"]
    derived = ["IMEAN", "SIGIMEAN"]
    for hkl, rrow in ref.items():
        grow = got[hkl]
        for lab in exact + derived:
            x, y = grow[lab], rrow[lab]
            xn, yn = x != x, y != y  # NaN checks
            assert xn == yn, f"{hkl} {lab}: missing-pattern mismatch ({x} vs {y})"
            if xn:
                continue
            if lab in exact:
                assert x == y, f"{hkl} {lab}: {x} != {y} (must be exact)"
            else:
                # relative tolerance at float32 precision
                assert abs(x - y) <= 1e-4 + 1e-5 * abs(y), (
                    f"{hkl} {lab}: {x} vs {y}"
                )
