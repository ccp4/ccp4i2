"""
crank2.compute_anomalous_scattering must compute f'/f'' from pure gemmi
(Cromer-Liberman) — no CCP4 install, no crossec binary — because it is reachable
from the API via the plugin_method endpoint and so runs on the CCP4-free server.

This is the "port a request-time computation to pip gemmi" pattern (the boundary
is wider than 'no CCP4 imports'; see tests/unit/slim/test_no_ccp4_native_imports.py
for the import half).
"""

import pytest

# crank2_script pulls in `future` and the plugin base; skip cleanly if absent
# rather than error (the computation itself is what we are testing).
pytest.importorskip("future")
crank2_mod = pytest.importorskip("ccp4i2.pipelines.crank2.script.crank2_script")


def _calc():
    # The method does not touch self.container, so a bare instance is enough.
    inst = crank2_mod.crank2.__new__(crank2_mod.crank2)
    return inst.compute_anomalous_scattering


def test_selenium_matches_reference():
    # Se @ 1.0332 A (~12 keV): gemmi Cromer-Liberman fp -2.510, fpp 0.550.
    out = _calc()("SE", 1.0332)
    assert out["fp"] == pytest.approx(-2.510, abs=0.02)
    assert out["fpp"] == pytest.approx(0.550, abs=0.01)


def test_case_insensitive_element():
    assert _calc()("se", 1.0332)["fp"] == pytest.approx(-2.510, abs=0.02)


def test_sulfur_reasonable():
    out = _calc()("S", 1.5418)  # Cu K-alpha
    assert out["fpp"] == pytest.approx(0.557, abs=0.02)


def test_unknown_element_errors():
    assert "error" in _calc()("FOO", 1.0)


def test_bad_wavelength_errors():
    assert "error" in _calc()("SE", 0)
    assert "error" in _calc()("SE", "not-a-number")
