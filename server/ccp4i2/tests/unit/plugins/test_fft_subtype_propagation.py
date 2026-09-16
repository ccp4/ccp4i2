"""fft types its output map from the input coefficients (issue #524).

fft's output map subtype can't be fixed in the def.xml -- it depends on the input
coefficients (a normal, difference, or anomalous map). Rather than a user param,
fft propagates the coefficients' subType (1/2/3) straight to the output map, which
is authoritative and needs no UI. Untyped input leaves the output untyped.
"""

import pytest

from ccp4i2.core.tasks import get_plugin_class


def _fft():
    return get_plugin_class("fft")()


@pytest.mark.parametrize("coeff_subtype", [1, 2, 3])
def test_output_map_inherits_coeff_subtype(coeff_subtype):
    p = _fft()
    p.container.inputData.FPHIIN.subType.set(coeff_subtype)
    p._propagate_map_subtype()
    assert int(p.container.outputData.MAPOUT.subType) == coeff_subtype


def test_default_input_gives_normal_output():
    # CMapCoeffsDataFile defaults to subType 1, so a plain coefficients file
    # yields a normal (1) map -- fft never leaves the output bare.
    p = _fft()
    p._propagate_map_subtype()
    assert int(p.container.outputData.MAPOUT.subType) == 1
