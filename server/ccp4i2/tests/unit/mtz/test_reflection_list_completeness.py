"""No reflection may be lost when a reflection list is rebuilt.

Every routine that writes a fresh MTZ and copies columns into it first
builds the output reflection list, and until now each one built it from
``gemmi.make_miller_array`` between the input's resolution limits. Those
limits are exclusive at floating-point precision, so the reflection that
defines the input's ``resolution_low()`` was left out of the list and
``copy_column`` silently dropped its observation. Every import, split and
merge in a pipeline shed one more reflection that way, until a refinement
was handed an observed reflection with no free flag.

The contract stated here: after a split, a merge, an import or a uniqueify,
the output holds every reflection of every input, and the resolution limits
are unchanged.
"""
import gemmi
import numpy as np
import pytest

from ccp4i2.core.CCP4Utils import (
    complete_reflection_list,
    merge_mtz_files,
    split_mtz_file,
)

CELLS = [
    (58.0, 64.4, 186.2, 90.0, 90.7, 90.0),
    (57.6, 64.4, 184.7, 90.0, 91.7, 90.0),
    (30.0, 40.0, 50.0, 90.0, 90.0, 90.0),
    (91.3, 91.3, 91.3, 90.0, 90.0, 90.0),
]


def _write(path, cell, columns, sg="P 1 21 1", hmax=3, kmax=3, lmax=4):
    mtz = gemmi.Mtz(with_base=True)
    mtz.spacegroup = gemmi.find_spacegroup_by_name(sg)
    mtz.set_cell_for_all(gemmi.UnitCell(*cell))
    mtz.add_dataset("data")
    for label, ctype, _ in columns:
        mtz.add_column(label, ctype)
    hkl = [(h, k, l) for h in range(0, hmax) for k in range(0, kmax) for l in range(1, lmax)]
    rows = [[h, k, l] + [vals[i % len(vals)] for _, _, vals in columns]
            for i, (h, k, l) in enumerate(hkl)]
    mtz.set_data(np.array(rows, dtype=np.float32))
    mtz.write_to_file(str(path))
    return path


def _hkls(path):
    m = gemmi.read_mtz_file(str(path))
    m.ensure_asu()
    return set(map(tuple, m.array[:, :3].astype(int).tolist()))


def _limits(path):
    m = gemmi.read_mtz_file(str(path))
    return m.resolution_low(), m.resolution_high()


def _assert_nothing_lost(inputs, output):
    """Every input's reflections are in the output, so the output's resolution
    range covers every input's range (it may be wider, never narrower)."""
    out_lo, out_hi = _limits(output)
    for path in inputs:
        lo, hi = _limits(path)
        assert out_lo >= lo * (1 - 1e-6), f"low-resolution limit of {path.name} lost"
        assert out_hi <= hi * (1 + 1e-6), f"high-resolution limit of {path.name} lost"
        assert _hkls(path) <= _hkls(output), f"reflections of {path.name} lost"


@pytest.mark.parametrize("cell", CELLS)
def test_the_helper_keeps_every_observed_reflection(tmp_path, cell):
    data = _write(tmp_path / "d.mtz", cell, [("F", "F", [1.0])])
    m = gemmi.read_mtz_file(str(data))
    m.ensure_asu()
    full = complete_reflection_list(m.cell, m.spacegroup, m.resolution_high(),
                                    m.resolution_low(), m.array[:, :3])
    assert _hkls(data) <= set(map(tuple, full.tolist()))
    assert full.dtype.kind == "i"
    assert (np.lexsort((full[:, 2], full[:, 1], full[:, 0])) == np.arange(len(full))).all()


@pytest.mark.parametrize("cell", CELLS)
def test_split_mtz_file_loses_nothing(tmp_path, cell):
    data = _write(tmp_path / "d.mtz", cell, [("FMEAN", "F", [1.0, 2.0]), ("SIGFMEAN", "Q", [0.1])])
    out = split_mtz_file(data, tmp_path / "s.mtz", {"FMEAN": "F", "SIGFMEAN": "SIGF"})
    _assert_nothing_lost([data], out)
    merged = gemmi.read_mtz_file(str(out))
    col = merged.column_with_label("F")
    observed = _hkls(data)
    for row in merged.array:
        if tuple(row[:3].astype(int).tolist()) in observed:
            assert not np.isnan(row[col.idx])


@pytest.mark.parametrize("cell", CELLS)
def test_merge_mtz_files_loses_nothing_from_any_input(tmp_path, cell):
    data = _write(tmp_path / "d.mtz", cell, [("F", "F", [1.0, 2.0]), ("SIGF", "Q", [0.1])])
    # A FreeR set that reaches beyond the data in every direction.
    freer = _write(tmp_path / "f.mtz", cell, [("FreeR_flag", "I", [0, 1, 2, 3, 4])],
                   hmax=5, kmax=5, lmax=8)
    out = merge_mtz_files(
        [{"path": data, "column_mapping": {"F": "F", "SIGF": "SIGF"}},
         {"path": freer, "column_mapping": {"FreeR_flag": "FreeR_flag"}}],
        tmp_path / "m.mtz")
    _assert_nothing_lost([data, freer], out)
    merged = gemmi.read_mtz_file(str(out))
    col = merged.column_with_label("FreeR_flag")
    flagged = _hkls(freer)
    for row in merged.array:
        if tuple(row[:3].astype(int).tolist()) in flagged:
            assert not np.isnan(row[col.idx])


def test_a_chain_of_merges_keeps_every_observation(tmp_path):
    data = _write(tmp_path / "d.mtz", CELLS[0], [("F", "F", [1.0, 2.0]), ("SIGF", "Q", [0.1])])
    freer = _write(tmp_path / "f.mtz", CELLS[0], [("FreeR_flag", "I", [0, 1, 2])])
    current = data
    for i in range(4):
        current = merge_mtz_files(
            [{"path": current, "column_mapping": {"F": "F", "SIGF": "SIGF"}},
             {"path": freer, "column_mapping": {"FreeR_flag": "FreeR_flag"}}],
            tmp_path / f"c{i}.mtz")
    _assert_nothing_lost([data], current)


@pytest.mark.parametrize("cell", CELLS)
def test_gemmi_split_mtz_import_loses_nothing(tmp_path, cell):
    from ccp4i2.lib.utils.formats.gemmi_split_mtz import gemmi_split_mtz
    data = _write(tmp_path / "d.mtz", cell, [("FP", "F", [1.0, 2.0]), ("SIGFP", "Q", [0.1])])
    out = gemmi_split_mtz(input_file_path=data, input_column_path="/*/*/[FP,SIGFP]",
                          preferred_dest=tmp_path / "imported.mtz")
    _assert_nothing_lost([data], out)


def test_freerflag_uniqueify_keeps_every_observation(tmp_path):
    from ccp4i2.core.tasks import get_plugin_class
    plugin_class = get_plugin_class("freerflag")
    assert plugin_class is not None
    directory = tmp_path / "freerflag"
    directory.mkdir()
    plugin = plugin_class(workDirectory=str(directory), name="freerflag")
    data = _write(tmp_path / "d.mtz", CELLS[0], [("F", "F", [1.0, 2.0]), ("SIGF", "Q", [0.1])])
    mtz = gemmi.read_mtz_file(str(data))
    mtz.ensure_asu()
    arr = np.array(mtz, copy=True)
    full, rows = plugin._complete_to_unique(mtz, arr)
    kept = set(map(tuple, np.asarray(full).astype(int).tolist()))
    assert _hkls(data) <= kept
    observed = {tuple(r[:3].astype(int).tolist()): r for r in arr}
    for row in rows:
        key = tuple(row[:3].astype(int).tolist())
        if key in observed:
            assert row[3] == pytest.approx(observed[key][3])
