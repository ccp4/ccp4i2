"""
Construction is not a user choice.

A freshly built container has had no user near it, so no parameter in it may
claim to have been explicitly set. That sounds obvious enough not to need a
test, which is exactly why it went wrong: `CString.__init__` assigned
`self.value`, assignment records `EXPLICITLY_SET`, and `CInt`, `CFloat` and
`CBoolean` declare no `__init__` at all and so never did. The asymmetry made
812 of 5,803 leaf parameters --- 14%, across 65 of 171 tasks --- claim a user
had chosen them before a user had seen them.

Nothing observed it, because `CString.isSet()` overrides the base to test the
value for emptiness rather than trust the state. Every consumer got the right
answer from the wrong reasoning, and the state was free to drift.

That cover disappears the moment anything reads `ValueState` directly, which
is what provenance-driven transmission does. `ParamsXmlHandler._is_explicitly_set`
already reads it and already answers True for these; it governs a different
branch than `controlParameters` takes, so a clean `params.xml` rested on which
of two disagreeing filters happened to run.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.core.base_object.cdata import ValueState
from ccp4i2.core.base_object.fundamental_types import (
    CBoolean, CFloat, CInt, CString,
)

SCALARS = (CString, CInt, CFloat, CBoolean)


@pytest.mark.parametrize("klass", SCALARS, ids=[k.__name__ for k in SCALARS])
def test_bare_construction_is_not_set(klass):
    """No scalar leaf may claim to be set before anything sets it."""
    obj = klass()
    assert obj.getValueState("value") == ValueState.NOT_SET, (
        f"{klass.__name__}() reports {obj.getValueState('value').name} --- "
        "construction is not a user choice"
    )
    assert obj.isSet() is False


def test_the_four_scalars_agree_with_each_other():
    """The bug was an asymmetry, so pin the symmetry, not just the value."""
    states = {k.__name__: k().getValueState("value") for k in SCALARS}
    assert len(set(states.values())) == 1, (
        f"scalar leaves disagree at construction: "
        f"{ {n: s.name for n, s in states.items()} }"
    )


def test_a_real_value_still_records_a_choice():
    """The fix must not cost us the ability to record a genuine one."""
    obj = CString("GEN_NEW")
    assert obj.getValueState("value") == ValueState.EXPLICITLY_SET
    assert obj.isSet() is True


def test_an_applied_default_survives_construction():
    """Seeding the empty value must not overwrite a default already applied."""
    obj = CString()
    obj.setDefault("auto")
    assert obj.value == "auto"
    assert obj.getValueState("value") == ValueState.DEFAULT
    assert obj.isDefault() is True


def test_no_task_has_explicitly_set_parameters_at_construction():
    """The exhaustive form: every task, every leaf, straight after building."""
    django = pytest.importorskip("django")
    from ccp4i2.core.tasks import TASKS, get_plugin_class

    offenders = []
    examined = built = 0
    for name in TASKS:
        try:
            plugin = get_plugin_class(name)(parent=None, name=name)
        except Exception:
            continue          # unbuildable tasks are a separate concern
        built += 1
        for section in ("controlParameters", "inputData", "outputData"):
            container = getattr(plugin.container, section, None)
            if container is None:
                continue
            for child in container.children():
                try:
                    state = child.getValueState("value")
                except Exception:
                    continue
                examined += 1
                if state == ValueState.EXPLICITLY_SET:
                    offenders.append(
                        f"{name}.{section}.{child._name} "
                        f"({type(child).__name__})"
                    )

    assert built > 100, f"only {built} tasks built --- the sweep proves little"
    assert not offenders, (
        f"{len(offenders)} of {examined} leaf parameters claim EXPLICITLY_SET "
        f"at construction, e.g. {offenders[:8]}"
    )
