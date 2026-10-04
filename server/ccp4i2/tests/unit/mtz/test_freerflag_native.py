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


# --- twin laws ---------------------------------------------------------------
#
# freerflag gives twin-related reflections the same flag by default
# (merohedral and pseudo-merohedral laws, obliquity 5 degrees), so that if the
# crystal is twinned the test set holds no reflection whose twin mate was
# refined against. The native port grouped by the space group alone.

TWINNABLE = I2_TOP / "demo_data" / "beta_blip" / "beta_blip_P3221.mtz"  # law -x,-y,z


def twin_mates_share(flags, hkl, mtz):
    """(pairs checked, pairs with different flags) over each reflection and its twin mates."""
    gops = mtz.spacegroup.operations()
    asu = gemmi.ReciprocalAsu(mtz.spacegroup)
    flag_of = {tuple(asu.to_asu(tuple(int(i) for i in h), gops)[0]): int(f)
               for h, f in zip(hkl, flags)}
    checked = differ = 0
    for op in gemmi.find_twin_laws(mtz.cell, mtz.spacegroup, 5.0, False):
        for h, f in zip(hkl, flags):
            mate = tuple(asu.to_asu(op.apply_to_hkl([int(i) for i in h]), gops)[0])
            if mate in flag_of and mate != tuple(asu.to_asu(tuple(int(i) for i in h), gops)[0]):
                checked += 1
                differ += flag_of[mate] != int(f)
    return checked, differ


def test_twin_mates_share_a_flag():
    twinnable = gemmi.read_mtz_file(str(TWINNABLE))
    hkl = _hkl(twinnable)
    flags, _ = assign_class_flags(hkl, twinnable.spacegroup, IRFRAC, cell=twinnable.cell)
    checked, differ = twin_mates_share(flags, hkl, twinnable)
    assert checked > 1000 and differ == 0
    # Grouped by the space group alone, as before, they mostly do not
    old, _ = assign_class_flags(hkl, twinnable.spacegroup, IRFRAC)
    assert twin_mates_share(old, hkl, twinnable)[1] > checked // 2
    assert 0.03 <= (flags == 0).mean() <= 0.07


def test_no_twin_law_no_change(mtz):
    # Gamma (P 21 21 21, distinct axes) has no twin law: the cell changes nothing
    assert not gemmi.find_twin_laws(mtz.cell, mtz.spacegroup, 5.0, False)
    a, _ = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC)
    b, _ = assign_class_flags(_hkl(mtz), mtz.spacegroup, IRFRAC, cell=mtz.cell)
    assert (a == b).all()
