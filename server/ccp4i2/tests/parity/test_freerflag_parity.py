"""
Parity test: gemmi/numpy-native free-R assignment vs the CCP4 freerflag binary.

freerflag's exact byte stream is *not* reproducible (it bottoms out in the
compiler's RNG, which is gfortran-version-dependent), so this checks *statistical*
parity rather than flag-for-flag: native and binary must land in the same regime
(valid range, same free fraction, ~same number of segments). Skipped when the
freerflag binary is absent (i.e. on a slim CCP4-free interpreter).
"""
import os
import shutil
import subprocess
import tempfile

import numpy as np
import gemmi
import pytest

from ccp4i2 import I2_TOP
from ccp4i2.wrappers.freerflag.script.freerflag import assign_class_flags

FREERFLAG = shutil.which("freerflag")
pytestmark = pytest.mark.skipif(FREERFLAG is None, reason="freerflag binary not on PATH")

MTZ = I2_TOP / "demo_data" / "gamma" / "HKLOUT_unmerged.mtz"
IRFRAC = 20


def _binary_flags(mtz_path, frac=0.05):
    out = tempfile.mktemp(suffix=".mtz")
    subprocess.run([FREERFLAG, "HKLIN", str(mtz_path), "HKLOUT", out],
                   input=f"FREERFRAC {frac}\nEND\n", text=True,
                   capture_output=True, timeout=60)
    m = gemmi.read_mtz_file(out)
    os.unlink(out)
    col = next(c for c in m.columns
               if c.type == 'I' and 'free' in c.label.lower())
    return np.array(col, dtype=int)


def test_native_and_binary_same_regime():
    m = gemmi.read_mtz_file(str(MTZ))
    hkl = np.array(m, copy=False)[:, :3].astype(int)
    native, _ = assign_class_flags(hkl, m.spacegroup, IRFRAC)
    binary = _binary_flags(MTZ)

    # Both must produce valid flags in [0, IRFRAC-1]
    assert 0 <= native.min() and native.max() <= IRFRAC - 1
    assert 0 <= binary.min() and binary.max() <= IRFRAC - 1

    # Both must hold out ~1/IRFRAC as the free set, and agree with each other
    fn = (native == 0).mean()
    fb = (binary == 0).mean()
    assert 0.03 <= fn <= 0.07, f"native free fraction {fn}"
    assert 0.03 <= fb <= 0.07, f"binary free fraction {fb}"
    assert abs(fn - fb) < 0.02, f"native {fn} vs binary {fb} free fraction differ too much"

    # Both should use the full set of segments (0..IRFRAC-1)
    assert len(np.unique(binary)) == IRFRAC
    assert len(np.unique(native)) == IRFRAC
