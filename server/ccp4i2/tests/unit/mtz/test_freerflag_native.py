"""
Unit tests for the gemmi/numpy-native free-R assignment (freerflag.assign_class_flags),
which replaced the CCP4 freerflag binary. These run CCP4-free and assert the
*semantics* we reproduce from freerflag (binning, free set = 0, one draw per
symmetry-equivalence class so equivalents share a flag, reproducibility).
Statistical parity with the binary is covered by tests/parity/test_freerflag_parity.py.
"""
from collections import defaultdict

import numpy as np
import gemmi
import pytest

from ccp4i2 import I2_TOP
from ccp4i2.wrappers.freerflag.script.freerflag import assign_class_flags

# unmerged demo data: has symmetry-equivalent reflections, so the sharing
# invariant is non-trivial here.
MTZ = I2_TOP / "demo_data" / "gamma" / "HKLOUT_unmerged.mtz"
IRFRAC = 20  # 1/0.05


@pytest.fixture
def mtz():
    return gemmi.read_mtz_file(str(MTZ))


def _hkl(mtz):
    return np.array(mtz, copy=False)[:, :3].astype(int)


def test_flags_in_valid_range(mtz):
    flags, _ = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC)
    assert flags.min() >= 0
    assert flags.max() <= IRFRAC - 1


def test_free_fraction_near_target(mtz):
    flags, _ = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC)
    frac = (flags == 0).mean()
    assert 0.03 <= frac <= 0.07, f"free fraction {frac} not near 1/{IRFRAC}"


def test_symmetry_equivalents_share_a_flag(mtz):
    flags, keys = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC)
    by_class = defaultdict(set)
    for f, k in zip(flags, keys):
        by_class[k].add(int(f))
    # every equivalence class has exactly one flag
    assert all(len(s) == 1 for s in by_class.values())
    # and there really are multi-reflection classes here (else the test is vacuous)
    sizes = defaultdict(int)
    for k in keys:
        sizes[k] += 1
    assert max(sizes.values()) > 1, "expected symmetry-equivalent reflections in unmerged data"


def test_reproducible_for_fixed_seed(mtz):
    a, _ = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC)
    b, _ = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC)
    assert np.array_equal(a, b)


def test_seed_changes_assignment(mtz):
    a, _ = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC, seed=1)
    b, _ = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC, seed=2)
    assert not np.array_equal(a, b)
