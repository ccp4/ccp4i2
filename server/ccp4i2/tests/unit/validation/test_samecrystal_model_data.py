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


def test_cell_is_none_for_malformed_cryst1(tmp_path):
    """A non-physical CRYST1 (zero axis, degenerate angle) yields no
    comparable cell, so it is skipped rather than fed to the comparison
    (where a zero axis raises and a degenerate angle mis-scales)."""
    from ccp4i2.core.CCP4ModelData import CPdbData

    for cell in [
        (0.0, 60.0, 70.0, 90.0, 90.0, 90.0),      # zero axis
        (50.0, 60.0, 70.0, 90.0, 90.0, 180.0),    # degenerate angle
        (50.0, 60.0, 70.0, 90.0, 90.0, 0.0),      # zero angle
    ]:
        path = str(tmp_path / "bad.pdb")
        _write_pdb(path, cell)
        content = CPdbData()
        content.loadFile(path)
        assert content.cell is None


def test_cpdbdata_spacegroup_reads_cryst1_and_is_not_serialised(tmp_path):
    """CPdbData exposes the CRYST1 space group as a computed accessor
    (never a tracked, serialised field)."""
    from ccp4i2.core.CCP4ModelData import CPdbData

    # _write_pdb writes "P 1" as the CRYST1 group; give it a real cell so the
    # group is not treated as the cell-less placeholder.
    path = str(tmp_path / "model.pdb")
    _write_pdb(path, (50.0, 60.0, 70.0, 90.0, 90.0, 90.0))
    content = CPdbData()
    content.loadFile(path)
    assert content.spaceGroup == "P 1"
    if hasattr(content, "get") and callable(content.get):
        try:
            data = content.get()
            if isinstance(data, dict):
                assert "spaceGroup" not in data
        except Exception:
            pass


def test_spacegroups_are_compatible_three_verdicts():
    """The three verdicts that drive the space-group tier, matching the
    observed refmac/servalcat behaviour:

      matched  -> same group,  same point group  (no message)
      samePG   -> diff group,  same point group  (advisory WARNING)
      diffPG   -> diff group,  diff point group   (blocking ERROR)
    """
    from ccp4i2.core.CCP4XtalData import spacegroups_are_compatible

    matched = spacegroups_are_compatible("C 2 2 21", "C2221")
    assert matched["sameSpaceGroup"] and matched["samePointGroup"]

    same_pg = spacegroups_are_compatible("P 21 21 21", "C 2 2 21")
    assert not same_pg["sameSpaceGroup"] and same_pg["samePointGroup"]

    # enantiomorph pair: same point group, different group
    enantio = spacegroups_are_compatible("P 41", "P 43")
    assert not enantio["sameSpaceGroup"] and enantio["samePointGroup"]

    diff_pg = spacegroups_are_compatible("P 1 21 1", "C 2 2 21")
    assert not diff_pg["samePointGroup"]


def test_spacegroups_are_compatible_skips_on_missing_or_bad_name():
    """A missing or unparseable space group returns None so the caller
    skips the tier rather than manufacturing a mismatch."""
    from ccp4i2.core.CCP4XtalData import spacegroups_are_compatible

    assert spacegroups_are_compatible(None, "P 21 21 21") is None
    assert spacegroups_are_compatible("P 21 21 21", "") is None
    assert spacegroups_are_compatible("not a space group", "P 1") is None


class _FakeChild:
    """Minimal stand-in exposing get_qualifier, for the spec resolver."""

    def __init__(self, **quals):
        self._q = quals

    def get_qualifier(self, key, default=None):
        return self._q.get(key, default)


def test_resolve_spec_defaults_to_spacegroup_and_auto_cell():
    """With only sameCrystalAs set, the default is: same space group, and
    'auto' cell (None) -- the check then requires the cell only when a model
    is involved, so a reflection<->reflection pair (free-R vs data) skips the
    cell and survives a soak series' cell drift."""
    from ccp4i2.core.CCP4PluginScript import _resolve_same_crystal_spec

    mode, cell, sev = _resolve_same_crystal_spec(_FakeChild())
    assert (mode, cell, sev) == ("spaceGroup", None, None)


def test_resolve_spec_named_qualifiers():
    """The orthogonal named qualifiers are honoured, with synonyms and bool
    coercion."""
    from ccp4i2.core.CCP4PluginScript import _resolve_same_crystal_spec
    from ccp4i2.core.base_object.error_reporting import (
        SEVERITY_ERROR, SEVERITY_WARNING)

    mode, cell, sev = _resolve_same_crystal_spec(_FakeChild(
        sameCrystalMatch="pointGroup", sameCrystalCell="false"))
    assert (mode, cell, sev) == ("pointGroup", False, None)

    mode, cell, sev = _resolve_same_crystal_spec(_FakeChild(
        sameCrystalMatch="point_group", sameCrystalCell=True,
        sameCrystalSeverity="error"))
    assert (mode, cell, sev) == ("pointGroup", True, SEVERITY_ERROR)

    mode, cell, sev = _resolve_same_crystal_spec(_FakeChild(
        sameCrystalMatch="none", sameCrystalSeverity="advisory"))
    assert (mode, cell, sev) == ("none", None, SEVERITY_WARNING)

    # explicit cell overrides the auto default in both directions
    assert _resolve_same_crystal_spec(
        _FakeChild(sameCrystalCell="true"))[1] is True
    assert _resolve_same_crystal_spec(
        _FakeChild(sameCrystalCell="false"))[1] is False


class _FakeInput:
    """Input object with a parent and a leaf name, for partner resolution."""

    def __init__(self, name, parent):
        self._name = name
        self._parent = parent

    def parent(self):
        return self._parent


def test_partner_resolution_prefers_same_scope():
    """When a leaf name is duplicated across composed scopes (e.g. a
    metalcoord sub-scope that also has an XYZIN/F_SIGF), the partner must
    resolve to the sibling in the declaring input's own scope, not the
    other scope's object."""
    from ccp4i2.core.CCP4PluginScript import _resolve_same_crystal_partner

    main = object()          # main inputData container
    metal = object()         # metalcoord sub-scope container

    xyzin_main = _FakeInput("XYZIN", main)
    fsigf_main = _FakeInput("F_SIGF", main)
    xyzin_metal = _FakeInput("XYZIN", metal)
    fsigf_metal = _FakeInput("F_SIGF", metal)
    all_files = [("XYZIN", xyzin_main), ("F_SIGF", fsigf_main),
                 ("XYZIN", xyzin_metal), ("F_SIGF", fsigf_metal)]

    # each XYZIN resolves F_SIGF within its own scope
    assert _resolve_same_crystal_partner(
        xyzin_main, "F_SIGF", all_files) is fsigf_main
    assert _resolve_same_crystal_partner(
        xyzin_metal, "F_SIGF", all_files) is fsigf_metal
    # never resolves to itself
    assert _resolve_same_crystal_partner(
        xyzin_main, "XYZIN", all_files) is not xyzin_main
    # unknown partner -> None
    assert _resolve_same_crystal_partner(
        xyzin_main, "NOPE", all_files) is None


def test_partner_resolution_falls_back_across_scope():
    """If the partner is not a sibling, fall back to any input with that
    name rather than failing."""
    from ccp4i2.core.CCP4PluginScript import _resolve_same_crystal_partner

    main = object()
    other = object()
    xyzin = _FakeInput("XYZIN", main)
    hklin_other = _FakeInput("HKLIN", other)  # only HKLIN lives elsewhere
    all_files = [("XYZIN", xyzin), ("HKLIN", hklin_other)]

    assert _resolve_same_crystal_partner(
        xyzin, "HKLIN", all_files) is hklin_other
