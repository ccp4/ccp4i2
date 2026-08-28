"""
aimless's override blocks rest on a coupling nothing enforces.

The four `*_OVERRIDE` toggles gate blocks of parameters, and inside each block
aimless asks `isDefault()` to decide what the user actually touched:

    if not par.OUTLIER_OVERRIDE:
        return
    if not par.OUTLIER_EMAX.isDefault():
        self.appendCommandScript("REJECT EMAX %f" % par.OUTLIER_EMAX)

That design is better than it first looks, and deliberately two-level: the
toggle is a use/don't-use gate, while provenance stays per parameter. Because
the parameters answer for themselves rather than inheriting from the flag,
unticking the box **keeps** a value the user determined --- it round-trips as
EXPLICITLY_SET and re-ticking recovers it. Had the flag carried the provenance,
unticking would have destroyed it.

But `isDefault()` is a state test, and `not isDefault()` is *true* for a
parameter that was never set at all. A parameter with no `<default>` in the
def.xml is NOT_SET at construction, so the guard would pass and aimless would
format an unset value into a command line.

It does not happen, because every parameter these guards protect declares a
default. That is the coupling: not a property of the code, a property of the
def.xml agreeing with it. This test is what makes them agree on purpose.

Pure Python -- no CCP4 binaries needed.
"""
import re
from pathlib import Path

import pytest

AIMLESS = Path("ccp4i2/wrappers/aimless/script")


def _guarded_parameter_names():
    """Every `par.X.isDefault()` aimless gates a command on."""
    source = (AIMLESS / "aimless.py").read_text()
    return sorted(set(re.findall(r"par\.([A-Z0-9_]+)\.isDefault\(\)", source)))


def test_the_guards_are_still_there():
    """If aimless stops using isDefault(), this file should be revisited."""
    names = _guarded_parameter_names()
    assert len(names) >= 6, f"expected aimless to gate on isDefault(); found {names}"


@pytest.mark.parametrize("name", _guarded_parameter_names())
def test_every_isDefault_guarded_parameter_declares_a_default(name):
    """No default means NOT_SET, and `not isDefault()` lets NOT_SET through."""
    def_xml = (AIMLESS / "aimless.def.xml").read_text()
    block = re.search(
        rf'<content id="{name}">.*?</content>', def_xml, re.S
    )
    assert block, f"{name} is guarded in aimless.py but absent from aimless.def.xml"
    assert "<default>" in block.group(0), (
        f"{name} has no <default>, so it is NOT_SET at construction and "
        f"`not {name}.isDefault()` is True --- aimless would format an unset "
        f"value into a command line"
    )
