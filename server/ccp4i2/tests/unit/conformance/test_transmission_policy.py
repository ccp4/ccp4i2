"""
Declaring when a parameter should reach the program.

Wrappers kept inventing this. phaser carries `requiredDefaultList`, a class
attribute naming parameters the program must be given even when they sit on
their default. aimless declares four `*_OVERRIDE` booleans in its def.xml,
each labelled "override default ...", whose real job is to let a user say
*the values below are mine* --- a provenance signal presented in the interface
as though it were science.

They solve opposite halves of one missing distinction: phaser needs to send a
default, aimless needs to withhold one. Both were reaching for a policy that
had nowhere to be declared.

`<sendWhen>` gives it somewhere. The def.xml qualifier parser is open --- any
tag inside `<qualifiers>` is retained --- so nothing needed to change to carry
it; what was missing was a named vocabulary and a reader that refuses a typo.

Nothing consumes this yet. These tests fix the meaning before anything
depends on it.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

from ccp4i2.core.base_object.cdata import CData, ValueState
from ccp4i2.core.base_object.fundamental_types import CFloat


def _param(state, value=0.05):
    """A CFloat parked in one of the three states."""
    obj = CFloat()
    if state is ValueState.DEFAULT:
        obj.setDefault(value)
    elif state is ValueState.EXPLICITLY_SET:
        obj.value = value
    assert obj.getValueState("value") == state
    return obj


# state, policy -> should it be sent
MATRIX = [
    (ValueState.NOT_SET,        CData.SEND_ALWAYS,    False),
    (ValueState.DEFAULT,        CData.SEND_ALWAYS,    True),
    (ValueState.EXPLICITLY_SET, CData.SEND_ALWAYS,    True),
    (ValueState.NOT_SET,        CData.SEND_IF_CHOSEN, False),
    (ValueState.DEFAULT,        CData.SEND_IF_CHOSEN, False),
    (ValueState.EXPLICITLY_SET, CData.SEND_IF_CHOSEN, True),
    (ValueState.NOT_SET,        CData.SEND_IF_SET,    False),
    (ValueState.EXPLICITLY_SET, CData.SEND_IF_SET,    True),
]


@pytest.mark.parametrize("state,policy,expected", MATRIX,
                         ids=[f"{s.name}-{p}" for s, p, _ in MATRIX])
def test_policy_matrix(state, policy, expected):
    assert _param(state).shouldSend(policy) is expected


def test_always_still_declines_an_unset_parameter():
    """`always` means whatever its provenance, not whatever its emptiness.

    A parameter nobody set and that carries no default has no value to send,
    and emitting it would put an empty string on the command line --- which is
    the class of silent nonsense this mechanism exists to stop.
    """
    assert _param(ValueState.NOT_SET).shouldSend(CData.SEND_ALWAYS) is False


def test_the_two_conditional_policies_differ_only_on_a_default():
    """If they agreed here there would be no reason for both to exist."""
    on_default = _param(ValueState.DEFAULT)
    assert on_default.shouldSend(CData.SEND_IF_SET) is True
    assert on_default.shouldSend(CData.SEND_IF_CHOSEN) is False


def test_the_declaration_beats_the_callers_baseline():
    """The point of declaring: the def.xml overrides what the wrapper assumes."""
    obj = _param(ValueState.DEFAULT)
    assert obj.shouldSend(CData.SEND_IF_CHOSEN) is False
    obj.set_qualifier("sendWhen", CData.SEND_ALWAYS)
    assert obj.shouldSend(CData.SEND_IF_CHOSEN) is True


def test_an_unknown_policy_raises_rather_than_falling_back():
    """A policy that silently means its opposite is the bug being removed."""
    obj = _param(ValueState.DEFAULT)
    obj.set_qualifier("sendWhen", "ifchosen")     # plausible casing slip
    with pytest.raises(ValueError, match="ifchosen"):
        obj.shouldSend()


def test_the_error_names_the_parameter():
    """A def.xml typo must say which parameter, or it is a hunt."""
    django = pytest.importorskip("django")
    from ccp4i2.core.tasks import get_plugin_class

    plugin = get_plugin_class("freerflag")(parent=None, name="freerflag")
    frac = plugin.container.controlParameters.FRAC
    frac.set_qualifier("sendWhen", "whenever")
    with pytest.raises(ValueError, match="FRAC"):
        frac.shouldSend()


def test_the_declared_set_is_exactly_what_phaser_used_to_hold_in_code():
    """What replaced `requiredDefaultList`, pinned against what it contained.

    The old lists were `['PART_VARI', 'PART_DEVI']` on phaser_EP_AUTO and
    `['PART_VARI', 'PART_DEVI', 'LLGM']` on phaser_EP_LLG --- so LLGM is
    required by the second and not the first, which is expressible only
    because each declares its own LLGM while PART_VARI and PART_DEVI are
    declared once on EP_AUTO and inherited.

    Any further declaration is a deliberate act and should land here with it.
    """
    django = pytest.importorskip("django")
    from ccp4i2.core.tasks import TASKS, get_plugin_class

    expected = {
        "phaser_EP_AUTO.keywords.PART_VARI",
        "phaser_EP_AUTO.keywords.PART_DEVI",
        "phaser_EP_LLG.keywords.PART_VARI",
        "phaser_EP_LLG.keywords.PART_DEVI",
        "phaser_EP_LLG.keywords.LLGM",
    }

    def walk(container, prefix, found):
        for child in container.children():
            name = f"{prefix}.{child._name}"
            if child.get_qualifier("sendWhen") is not None:
                found.add(name)
            if child.dataOrder():
                walk(child, name, found)

    found = set()
    for name in TASKS:
        try:
            plugin = get_plugin_class(name)(parent=None, name=name)
        except Exception:
            continue
        for section in ("controlParameters", "inputData", "outputData", "keywords"):
            container = getattr(plugin.container, section, None)
            if container is not None:
                walk(container, f"{name}.{section}", found)

    assert found == expected, (
        f"declared but not expected: {sorted(found - expected)}; "
        f"expected but not declared: {sorted(expected - found)}"
    )


def test_phaser_EP_AUTO_does_not_require_LLGM():
    """The asymmetry is the whole reason this had to be per-task."""
    django = pytest.importorskip("django")
    from ccp4i2.core.tasks import get_plugin_class

    auto = get_plugin_class("phaser_EP_AUTO")(parent=None, name="a")
    llg = get_plugin_class("phaser_EP_LLG")(parent=None, name="l")
    assert auto.container.keywords.LLGM.sendWhen(CData.SEND_IF_CHOSEN) == CData.SEND_IF_CHOSEN
    assert llg.container.keywords.LLGM.sendWhen(CData.SEND_IF_CHOSEN) == CData.SEND_ALWAYS
