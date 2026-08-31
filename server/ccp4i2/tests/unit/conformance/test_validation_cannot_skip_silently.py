"""A check that raised is a check that did not happen.

`CData.validity()` recursed into children inside `try/except Exception: pass`,
commented "Don't let one child's validation failure stop others". The first
half was right; the second was not. A child whose validation *raised* was
reported as **valid**, so a job could run on an assurance nobody had given.

Unlike the swallow in `_apply_metadata_attributes`, this one was firing. Making
it report surfaced three defects that had been invisible for as long as they
have existed:

  1. `maxLength` and `minLength` arrive from a def.xml as *text* and were never
     coerced, so `len(value) > maxLength` raised TypeError. 99 numeric
     qualifiers across the registry were strings, and every bounds check on
     them had therefore never run.
  2. A def.xml writing `<min>None</min>` to mean "no bound" produced the string
     'None', which is truthy --- so the comparison fired and raised, by a
     second route.
  3. `CRangeSelection.validity(self, arg)` required a positional argument that
     the framework does not pass, so *no* CRangeSelection had ever been
     validated.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

django = pytest.importorskip("django")


def test_numeric_qualifiers_are_numbers():
    """Whatever route they arrived by: def.xml, class metadata, CList subitem."""
    from ccp4i2.core.tasks import TASKS, get_plugin_class

    NUMERIC = ("min", "max", "minLength", "maxLength",
               "listMinLength", "listMaxLength")
    strings, checked = [], 0
    for name in TASKS:
        try:
            plugin = get_plugin_class(name)(parent=None, name=name)
        except Exception:
            continue
        seen = set()

        def walk(obj, path, depth=0):
            nonlocal checked
            if id(obj) in seen or depth > 7:
                return
            seen.add(id(obj))
            for key in NUMERIC:
                value = obj.get_qualifier(key) if hasattr(obj, "get_qualifier") else None
                if value is not None:
                    checked += 1
                    if isinstance(value, str):
                        strings.append(f"{name}.{path}: {key}={value!r}")
            for child in obj.children():
                walk(child, f"{path}.{child._name}", depth + 1)

        walk(plugin.container, "container")

    assert checked > 200, f"only {checked} qualifiers seen --- proves little"
    assert not strings, (
        f"{len(strings)} numeric qualifiers are strings, so every bounds check "
        f"on them raises and is skipped: {strings[:6]}"
    )


def test_a_length_bound_is_enforced():
    """LIGAND_CODE declares maxLength 3. Before, that never ran."""
    from ccp4i2.core.tasks import get_plugin_class

    plugin = get_plugin_class("SubstituteLigand")(parent=None, name="s")
    code = plugin.container.controlParameters.LIGAND_CODE
    assert code.get_qualifier("maxLength") == 3

    code.value = "TOOLONG"
    assert list(code.validity().entries()), "a 7-character value passed a max of 3"

    code.value = "DRG"
    assert not list(code.validity().entries()), "a valid value was rejected"


def test_a_range_selection_is_validated():
    """Its signature made it uncallable by the framework, so it never was."""
    from ccp4i2.core.CCP4Data import CRangeSelection

    bad = CRangeSelection()
    bad.value = "1-10,oops"
    assert list(bad.validity().entries()), "a malformed range was accepted"

    good = CRangeSelection()
    good.value = "1-10,15"
    assert not list(good.validity().entries()), "a valid range was rejected"


def test_an_unset_range_selection_is_not_reported_as_malformed():
    """An empty selection is not a malformed one.

    An unset CString holds '' rather than None, so without care the syntax
    check runs on '' and calls it malformed --- reporting every unset
    excludeSelection, a warning that is always wrong and therefore worse than
    none.

    Being *unset* may still be reported where the declaration does not allow
    it; that is a different statement, and the one `allowUndefined` governs.
    """
    from ccp4i2.core.tasks import get_plugin_class

    plugin = get_plugin_class("pointless")(parent=None, name="p")
    selection = plugin.container.inputData.UNMERGEDFILES[0].excludeSelection
    assert not list(selection.validity().entries()), \
        "an unset, optional range selection was reported"


def test_a_child_that_cannot_be_validated_is_reported():
    """The elevation itself: siblings still validate, and the failure shows."""
    from ccp4i2.core.base_object.ccontainer import CContainer
    from ccp4i2.core.base_object.fundamental_types import CInt

    class Exploding(CInt):
        def validity(self):
            raise RuntimeError("boom")

    container = CContainer(name="root")
    container.ok = CInt(name="ok")
    container.bad = Exploding(name="bad")

    reported = [str(e) for e in container.validity().entries()]
    assert any("could not validate" in r and "RuntimeError" in r for r in reported), \
        f"a child whose validation raised was reported as valid: {reported}"
