"""
Unit tests for CMtzData.matthewsCoeff() — the gemmi/numpy-native Matthews
coefficient calculation that replaced the CCP4 matthews_coef binary.

These run on a slim, CCP4-free interpreter (only gemmi needed): they assert the
internal consistency of the result rather than parity with the binary (that is
covered by tests/parity/test_matthews_coeff_parity.py).
"""
import pytest
import gemmi

from ccp4i2 import I2_TOP
from ccp4i2.core.CCP4XtalData import CMtzDataFile
from ccp4i2.core._matthews_data import SOLVENT_K

MTZ = I2_TOP / "demo_data" / "gamma" / "freeR.mtz"


@pytest.fixture
def content():
    # Keep the CMtzDataFile alive for the test (the fileContent's cell is lost if
    # the parent file object is garbage-collected), so yield via the file object.
    f = CMtzDataFile()
    f.setFullPath(str(MTZ))
    f.loadFile()
    yield f.fileContent


def test_runs_without_binary(content):
    rv = content.matthewsCoeff(molWt=8000.0, polymerMode="P")
    assert rv["cell_volume"] > 0
    assert len(rv["results"]) >= 1
    for r in rv["results"]:
        assert set(r) == {"nmol_in_asu", "matth_coef", "percent_solvent", "prob_matth"}


def test_vm_matches_definition(content):
    """Vm must equal V / (nmol * nsym * MW) exactly."""
    MW = 8000.0
    rv = content.matthewsCoeff(molWt=MW, polymerMode="P")
    V = rv["cell_volume"]
    nsym = len(gemmi.SpaceGroup(str(content.spaceGroup)).operations())
    for r in rv["results"]:
        expected = V / (r["nmol_in_asu"] * nsym * MW)
        assert r["matth_coef"] == pytest.approx(expected, rel=1e-9)


def test_solvent_constant_and_positive(content):
    """percent_solvent = 100*(1 - k/Vm), strictly positive, mode-dependent k."""
    for mode in ("P", "D", "C"):
        rv = content.matthewsCoeff(molWt=6000.0, polymerMode=mode)
        for r in rv["results"]:
            expected = 100.0 * (1.0 - SOLVENT_K[mode] / r["matth_coef"])
            assert r["percent_solvent"] == pytest.approx(expected, abs=1e-9)
            assert r["percent_solvent"] > 0.0


def test_vm_and_solvent_monotonic(content):
    rv = content.matthewsCoeff(molWt=6000.0, polymerMode="P")
    vms = [r["matth_coef"] for r in rv["results"]]
    solv = [r["percent_solvent"] for r in rv["results"]]
    assert vms == sorted(vms, reverse=True)      # Vm decreases with nmol
    assert solv == sorted(solv, reverse=True)     # solvent decreases with nmol


def test_probabilities_normalised(content):
    for mode in ("P", "D", "C"):
        rv = content.matthewsCoeff(molWt=8000.0, polymerMode=mode)
        total = sum(r["prob_matth"] for r in rv["results"])
        assert total == pytest.approx(1.0, abs=1e-9)
        assert all(0.0 <= r["prob_matth"] <= 1.0 for r in rv["results"])


def test_unknown_mode_falls_back_to_protein(content):
    assert (content.matthewsCoeff(molWt=8000.0, polymerMode="ZZZ")["results"]
            == content.matthewsCoeff(molWt=8000.0, polymerMode="P")["results"])


def test_nres_estimates_weight(content):
    """nRes path estimates MW = 112.5 * nRes."""
    by_nres = content.matthewsCoeff(nRes=80, polymerMode="P")
    by_mw = content.matthewsCoeff(molWt=112.5 * 80, polymerMode="P")
    assert [r["matth_coef"] for r in by_nres["results"]] == \
           pytest.approx([r["matth_coef"] for r in by_mw["results"]])


def test_zero_weight_raises(content):
    from ccp4i2.core.base_object.error_reporting import CException
    with pytest.raises(CException):
        content.matthewsCoeff(molWt=0.0)
