"""BUSTER is licensed separately and not installed in CI, so what can be tested
is what the wrapper builds. The free set is optional in the interface, but the
CAD keywords always named its column, so a job without one failed in CAD."""
from ccp4i2.wrappers.buster.script.buster import cad_keywords


def test_without_free_set_names_no_free_column():
    text = cad_keywords(intensities=False, free_set=False)
    assert "FREER" not in text
    assert text.splitlines()[0] == "LABIN FILE 1 E1=F_SIGF_F E2=F_SIGF_SIGF"


def test_with_free_set_and_intensities():
    lines = cad_keywords(intensities=True, free_set=True).splitlines()
    assert lines[0].endswith("E5=FREERFLAG_FREER")
    assert lines[1].endswith("E5=FREER")
