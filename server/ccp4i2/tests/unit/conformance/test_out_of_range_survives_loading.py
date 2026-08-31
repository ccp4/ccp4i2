"""A value the file records must not be replaced by a default in silence.

`_import_parameter_value` ended with `except Exception: return False`. Two
different failures raise ValueError inside it and only one is a conversion
problem:

  - "that is not a number" --- a malformed file. Returning False is right.
  - a *bounds* violation --- the assignment-time min/max check refusing a value
    somebody chose. Returning False left the parameter at its default, so a
    params.xml holding FRAC=1.5 loaded as FRAC=0.05 and said nothing, and a
    rerun used a value the file does not record.

The value is now kept and `validity()` reports it, which it does correctly once
the value is actually stored.

That numeric types reject at *assignment* while string types accept and report
at *validation* is a real inconsistency in the model, not addressed here. This
only stops the difference costing a silent substitution on load.

Pure Python -- no CCP4 binaries needed.
"""
import os
import tempfile

import pytest

django = pytest.importorskip("django")


def _round_trip(value):
    """Save freerflag, rewrite FRAC to *value* in the file, load it back."""
    from ccp4i2.core.tasks import get_plugin_class

    work = tempfile.mkdtemp()
    path = os.path.join(work, "params.xml")
    saved = get_plugin_class("freerflag")(workDirectory=work, name="f")
    saved.container.controlParameters.FRAC = 0.3
    saved.saveDataToXml(path)

    with open(path) as handle:
        body = handle.read()
    assert "<FRAC>0.3</FRAC>" in body
    with open(path, "w") as handle:
        handle.write(body.replace("<FRAC>0.3</FRAC>", f"<FRAC>{value}</FRAC>"))

    loaded = get_plugin_class("freerflag")(workDirectory=work, name="f")
    loaded.loadDataFromXml(path)
    return loaded.container.controlParameters.FRAC


def test_an_in_range_value_loads():
    """The control: without this the rest proves nothing."""
    assert _round_trip(0.7).value == 0.7


def test_an_out_of_range_value_is_kept_not_replaced():
    """FRAC declares min 0.0 max 1.0. 1.5 used to load as 0.05."""
    assert _round_trip(1.5).value == 1.5


def test_and_validity_says_why_it_is_wrong():
    reported = [str(e.get("details")) for e in _round_trip(1.5).validity().entries()]
    assert any("above maximum" in r for r in reported), \
        f"an out-of-range value loaded without being reported: {reported}"


def test_a_value_that_is_not_a_number_is_still_rejected():
    """The other half of the ValueError: a malformed file is not a bounds
    violation, and must not be stored as though it were."""
    frac = _round_trip("not-a-number")
    assert frac.value == 0.05, "garbage was stored rather than refused"
