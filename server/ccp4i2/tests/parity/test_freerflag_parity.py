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
from ccp4i2.wrappers.freerflag.script.freerflag import (
    assign_class_flags,
    complete_class_flags,
)

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


# --- completing an existing free set ----------------------------------------
#
# The test above covers generating a free set from nothing. Completing one --
# what import_merged does whenever the input already carries flags -- has its
# own contract, and it is the one that matters scientifically: a reflection
# that was used in refinement must never end up in the test set.
#
# Inferring the segment count from max(existing) instead made the result depend
# on how the input was labelled: an mmCIF status column gives a binary free
# set, so the output was binary where CCP4 produces twenty segments.

COMPLETE_FRACTION = 0.05


def _mtz_with_partial_free_set(tmp_path, missing=0.02):
    """A *merged* MTZ carrying a binary free set with some reflections unflagged.

    Binary because that is what an mmCIF `status` column becomes, and it is the
    case that exposed the defect. Merged because freerflag uniquifies an
    unmerged list, and then the rows it returns are not the rows it was given.
    """
    m = gemmi.read_mtz_file(str(MTZ))
    data = np.array(m, copy=False)

    # one row per unique Miller index, in ASU
    seen = {}
    for row in data:
        key = (int(row[0]), int(row[1]), int(row[2]))
        seen.setdefault(key, row)
    hkl = np.array(sorted(seen), dtype=int)
    intensity = np.array([seen[tuple(h)][3] for h in hkl], dtype=float)
    sigma = np.array([seen[tuple(h)][4] for h in hkl], dtype=float)

    flags, _ = assign_class_flags(hkl, m.spacegroup, IRFRAC)
    binary = (flags != 0).astype(float)          # 0 = free, 1 = working
    rng = np.random.default_rng(11)
    binary[rng.random(len(binary)) < missing] = np.nan

    out = gemmi.Mtz(with_base=True)
    out.spacegroup = m.spacegroup
    out.set_cell_for_all(m.cell)
    out.add_dataset("free")
    for label, ctype in (("I", "J"), ("SIGI", "Q"), ("FREER", "I")):
        out.add_column(label, ctype)
    out.set_data(np.column_stack([hkl.astype(float), intensity, sigma, binary]))
    path = tmp_path / "partial_free.mtz"
    out.write_to_file(str(path))
    return str(path), {tuple(h): f for h, f in zip(hkl, binary)}


def _flags_by_index(mtz_path):
    m = gemmi.read_mtz_file(str(mtz_path))
    data = np.array(m, copy=False)
    idx = [c.label for c in m.columns].index("FREER")
    return {(int(r[0]), int(r[1]), int(r[2])): r[idx] for r in data}


def _binary_completed(mtz_path):
    """Complete *mtz_path*'s free set with the CCP4 binary; flags by index."""
    out = tempfile.mktemp(suffix=".mtz")
    subprocess.run([FREERFLAG, "HKLIN", str(mtz_path), "HKLOUT", out],
                   input="COMPLETE FREE=FREER\nEND\n", text=True,
                   capture_output=True, timeout=60)
    flags = _flags_by_index(out)
    os.unlink(out)
    return flags


def test_completion_keeps_every_free_reflection_free(tmp_path):
    """The contract: the test set survives completion intact."""
    path, before = _mtz_with_partial_free_set(tmp_path)
    after = _binary_completed(path)

    was_free = [h for h, f in before.items() if f == 0]
    assert len(was_free) > 100, "need a free set worth checking"
    kept = [h for h in was_free if after.get(h) == 0]
    assert len(kept) == len(was_free), (
        f"{len(was_free) - len(kept)} of {len(was_free)} free reflections "
        "stopped being free"
    )


def test_completion_never_moves_a_working_reflection_into_the_test_set(tmp_path):
    """The one that would quietly invalidate an R-free."""
    path, before = _mtz_with_partial_free_set(tmp_path)
    after = _binary_completed(path)

    moved = [h for h, f in before.items() if f == 1 and after.get(h) == 0]
    assert not moved, (
        f"{len(moved)} reflections already used in refinement became free"
    )


def test_completion_segments_the_working_set(tmp_path):
    """CCP4 re-segments rather than preserving a binary labelling: the output
    spans many segments even though the input held only 0 and 1."""
    path, _ = _mtz_with_partial_free_set(tmp_path)
    after = _binary_completed(path)

    values = np.array([v for v in after.values() if not np.isnan(v)])
    distinct = len(np.unique(values))
    assert distinct > 2, (
        f"binary in, {distinct} distinct values out -- the working set is "
        "re-segmented, so a later re-partition has segments to work with"
    )
    free_fraction = float(np.mean(values == 0))
    assert 0.03 <= free_fraction <= 0.08, f"free fraction {free_fraction}"


def _native_completed(mtz_path):
    """Complete the same file with the native implementation; flags by index."""
    m = gemmi.read_mtz_file(str(mtz_path))
    data = np.array(m, copy=False)
    hkl = data[:, :3].astype(int)
    idx = [c.label for c in m.columns].index("FREER")
    flags = complete_class_flags(hkl, m.spacegroup, data[:, idx].astype(float), IRFRAC)
    return {tuple(h): float(f) for h, f in zip(hkl, flags)}


@pytest.mark.parametrize("complete", [_binary_completed, _native_completed],
                         ids=["binary", "native"])
def test_completion_contract_holds_for_both(tmp_path, complete):
    """The same three questions, asked of the CCP4 binary and of ours.

    This is the parity that matters. The generation test above compares a fresh
    assignment; completion has its own contract, and it is the one with an
    R-free at stake.
    """
    path, before = _mtz_with_partial_free_set(tmp_path)
    after = complete(path)

    was_free = [h for h, f in before.items() if f == 0]
    was_working = [h for h, f in before.items() if f == 1]
    assert len(was_free) > 100 and len(was_working) > 100

    lost = [h for h in was_free if after.get(h) != 0]
    assert not lost, f"{len(lost)} free reflections stopped being free"

    moved = [h for h in was_working if after.get(h) == 0]
    assert not moved, f"{len(moved)} working reflections became free"

    values = np.array([v for v in after.values() if not np.isnan(v)])
    assert len(np.unique(values)) > 2, "the working set must be re-segmented"
    assert 0.03 <= float(np.mean(values == 0)) <= 0.08


def test_native_and_binary_complete_to_the_same_free_fraction(tmp_path):
    path, _ = _mtz_with_partial_free_set(tmp_path)
    binary = np.array([v for v in _binary_completed(path).values() if not np.isnan(v)])
    native = np.array([v for v in _native_completed(path).values() if not np.isnan(v)])

    fb = float(np.mean(binary == 0))
    fn = float(np.mean(native == 0))
    assert abs(fn - fb) < 0.02, f"native {fn:.4f} vs binary {fb:.4f}"
