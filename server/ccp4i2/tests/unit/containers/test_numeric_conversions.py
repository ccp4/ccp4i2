"""
CInt and CFloat must be usable wherever Python expects a number.

Found by C1 (docs/error-handling-remediation.md): chainsaw's
processOutputFiles() computes a sequence identity as

    float(outputData.NO_CONSERVED) / float(outputData.NUMRES_MODEL)

which raised ``TypeError: float() argument must be a string or a real number,
not 'CInt'`` on every run. The job was recorded as Finished, because the
exception was caught and printed to the server console.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.core.base_object.fundamental_types import CFloat, CInt


def test_float_of_a_cint():
    """The chainsaw case, exactly."""
    conserved, total = CInt(), CInt()
    conserved.set(42)
    total.set(84)
    assert float(conserved) / float(total) == pytest.approx(0.5)


def test_int_of_a_cint_is_unchanged():
    n = CInt()
    n.set(7)
    assert int(n) == 7


def test_a_cint_can_index_a_sequence():
    i = CInt()
    i.set(2)
    assert "abcde"[i] == "c"
    assert list(range(i)) == [0, 1]


def test_int_of_a_cfloat_truncates():
    x = CFloat()
    x.set(2.75)
    assert int(x) == 2
    assert float(x) == pytest.approx(2.75)


@pytest.mark.parametrize("factory,setter,expected", [
    (CInt, 3, 3.0),
    (CFloat, 3.5, 3.5),
])
def test_both_types_convert_to_float(factory, setter, expected):
    obj = factory()
    obj.set(setter)
    assert float(obj) == pytest.approx(expected)
