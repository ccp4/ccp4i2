"""SIMBAD's lattice table carries raw floats (0.0058999999999969) that push
its last columns out of the box; the report rounds them to four decimals and
leaves everything else alone."""
from ccp4i2.wrappers.SIMBAD.script.SIMBAD_report import _tidy_number


def test_long_floats_are_rounded():
    assert _tidy_number("0.0058999999999969") == "0.0059"
    assert _tidy_number(" 1.2721000000000018 ") == "1.2721"


def test_other_cells_are_untouched():
    for text in ("1gyu", "2347.033", "90.0", "0.891", "", None):
        assert _tidy_number(text) == text
