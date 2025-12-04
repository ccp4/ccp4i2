from math import isclose
from pathlib import Path
from os import environ
from ...wrappers.refmacat import refmacat


def test_1rxf():
    hklin = Path(environ["CCP4"], "examples", "data", "1rxf.mtz")
    xyzin = Path(environ["CCP4"], "examples", "data", "1rxf_randomise.pdb")
    result = refmacat(hklin, xyzin)
    assert isclose(result.initial_r_work, 0.3352, abs_tol=0.01)
    assert isclose(result.initial_r_free, 0.3390, abs_tol=0.01)
    assert isclose(result.initial_fsc, 0.9025, abs_tol=0.01)
    assert result.r_work < result.initial_r_work
    assert result.r_free < result.initial_r_free
    assert result.fsc > result.initial_fsc
    assert isclose(result.data_completeness, 95.261, abs_tol=0.01)
    assert isclose(result.resolution_low, 53.106, abs_tol=0.01)
    assert isclose(result.resolution_high, 1.501, abs_tol=0.01)
    result.delete_files()
