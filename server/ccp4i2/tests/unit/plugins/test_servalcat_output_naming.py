"""servalcat harvests its outputs by resolving names, not hard-coding them.

servalcat's SPA output basenames drifted between versions: 0.4 renamed the
reflection file ``refined_diffmap.mtz`` -> ``refined_maps.mtz`` and dropped the
``_diffmap`` infix from the normalised maps. The wrapper still named the old
files, so a *successful* SPA refinement then failed to split a file that did not
exist, recorded a warning, and returned FAILED anyway (CryoMapMR job 15).
``_firstExistingInWork`` resolves against the names we know so the harvest
survives the drift; these tests pin that behaviour without running servalcat.
"""

from ccp4i2.wrappers.servalcat.script.servalcat import servalcat


def _probe(tmp_path):
    p = object.__new__(servalcat)
    p.getWorkDirectory = lambda: str(tmp_path)
    return p


def test_prefers_the_current_name(tmp_path):
    (tmp_path / "refined_maps.mtz").write_bytes(b"x")
    p = _probe(tmp_path)
    assert p._firstExistingInWork(
        ["refined_maps.mtz", "refined_diffmap.mtz"]
    ) == str(tmp_path / "refined_maps.mtz")


def test_falls_back_to_the_legacy_name(tmp_path):
    # Only the old name on disk (an older servalcat) -> still found.
    (tmp_path / "refined_diffmap.mtz").write_bytes(b"x")
    p = _probe(tmp_path)
    assert p._firstExistingInWork(
        ["refined_maps.mtz", "refined_diffmap.mtz"]
    ) == str(tmp_path / "refined_diffmap.mtz")


def test_first_listed_wins_when_both_exist(tmp_path):
    (tmp_path / "refined_maps.mtz").write_bytes(b"x")
    (tmp_path / "refined_diffmap.mtz").write_bytes(b"x")
    p = _probe(tmp_path)
    assert p._firstExistingInWork(
        ["refined_maps.mtz", "refined_diffmap.mtz"]
    ).endswith("refined_maps.mtz")


def test_none_when_nothing_matches(tmp_path):
    p = _probe(tmp_path)
    assert p._firstExistingInWork(["refined_maps.mtz", "refined_diffmap.mtz"]) is None
