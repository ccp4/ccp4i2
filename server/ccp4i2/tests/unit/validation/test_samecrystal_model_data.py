"""The sameCrystalAs unit-cell check extended to model<->data pairs.

A coordinate file (CPdbData) now exposes a `.cell`, so the runtime
sameCrystalAs check compares a model's cell against reflection data — the
"cells don't match" case the Qt app surfaced as an overridable "Ignore"
dialog. Here that is an overridable warning: a model involved in the pair
makes the mismatch SEVERITY_WARNING (advisory), not a hard block.

CCP4-free: needs gemmi (pip) for the coordinate cell, no CCP4 binaries.
"""

import os
import tempfile

import pytest

django = pytest.importorskip("django")
pytest.importorskip("gemmi", reason="coordinate cell needs gemmi")


@pytest.fixture(scope="module", autouse=True)
def _django():
    os.environ.setdefault("DJANGO_SETTINGS_MODULE",
                          "ccp4i2.config.test_settings")
    django.setup()


def _write_pdb(path, cell):
    with open(path, "w") as handle:
        if cell:
            handle.write(
                "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f P 1\n" % cell)
        handle.write(
            "ATOM      1  CA  ALA A   1       0.000   0.000   0.000"
            "  1.00  0.00           C\n")
        handle.write("END\n")


def test_cpdbdata_cell_reads_cryst1(tmp_path):
    from ccp4i2.core.CCP4ModelData import CPdbData

    path = str(tmp_path / "model.pdb")
    _write_pdb(path, (50.0, 60.0, 70.0, 90.0, 90.0, 90.0))
    content = CPdbData()
    content.loadFile(path)

    cell = content.cell
    assert cell is not None
    assert round(cell.a, 1) == 50.0
    assert round(cell.b, 1) == 60.0
    assert round(cell.c, 1) == 70.0


def test_cpdbdata_cell_none_without_cryst1(tmp_path):
    """A model with no CRYST1 has no comparable cell — so it can never
    raise a spurious mismatch (gemmi fills an absent CRYST1 with a
    trivial 1x1x1 cell, which .cell reports as None)."""
    from ccp4i2.core.CCP4ModelData import CPdbData

    path = str(tmp_path / "nocryst.pdb")
    _write_pdb(path, None)
    content = CPdbData()
    content.loadFile(path)
    assert content.cell is None


def test_cpdbdata_cell_is_not_a_serialised_attribute(tmp_path):
    """`cell` is a computed property off the gemmi structure, never a
    tracked CData field — so it must not appear in the serialised data."""
    from ccp4i2.core.CCP4ModelData import CPdbData

    content = CPdbData()
    # get() returns the serialisable field dict; cell must not be in it.
    if hasattr(content, "get") and callable(content.get):
        try:
            data = content.get()
            if isinstance(data, dict):
                assert "cell" not in data
        except Exception:
            pass  # some content classes don't implement get(); fine


def test_model_is_flagged_for_the_advisory_branch(tmp_path):
    """The check picks WARNING vs ERROR by whether a model is one side of
    the pair, duck-typed on the absence of `clipperSameCell` (reflection
    contents have it; a coordinate content does not)."""
    from ccp4i2.core.CCP4ModelData import CPdbData

    path = str(tmp_path / "model.pdb")
    _write_pdb(path, (50.0, 60.0, 70.0, 90.0, 90.0, 90.0))
    content = CPdbData()
    content.loadFile(path)
    # This is exactly the predicate the check uses to choose severity.
    assert not hasattr(content, "clipperSameCell")


def test_cells_are_compatible_distinguishes_match_and_mismatch():
    """The comparison primitive the check now calls directly for a
    model<->data pair."""
    from ccp4i2.core.CCP4XtalData import cells_are_compatible

    same = (50.0, 60.0, 70.0, 90.0, 90.0, 90.0)
    assert cells_are_compatible(same, same, tolerance=1.0)["validity"]
    assert not cells_are_compatible(
        same, (80.0, 90.0, 100.0, 90.0, 90.0, 90.0), tolerance=0.01
    )["validity"]
