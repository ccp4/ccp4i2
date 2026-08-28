"""Qualifiers are shared by declaration and private by instance.

The merged `_qualifiers_template` depends only on the class, so it is merged
once per class rather than once per object. That is safe exactly as long as no
instance can reach another's qualifiers --- or the cache itself.

The risk is not theoretical. `MakeLink` sets `allowUndefined` on around twenty
parameters *conditionally, inside a method*, so those instances genuinely
diverge from their declaration; `xia2_xds` rewrites `enumerators` and `default`
at runtime. If a `set_qualifier` on one container could reach another, one
job's state would leak into the next --- and because containers are rebuilt per
API request, "the next" means the next parameter edit by anyone.

These pin the isolation the cache depends on, which is otherwise only visible
by reading `CData.__init__` and believing it.

Pure Python -- no CCP4 binaries needed.
"""
import pytest

django = pytest.importorskip("django")


@pytest.fixture
def two_freerflags():
    from ccp4i2.core.tasks import get_plugin_class
    klass = get_plugin_class("freerflag")
    a = klass(parent=None, name="a")
    b = klass(parent=None, name="b")
    return (a.container.controlParameters.GEN_MODE,
            b.container.controlParameters.GEN_MODE)


def test_two_instances_do_not_share_the_qualifiers_dict(two_freerflags):
    a, b = two_freerflags
    assert a._qualifiers is not b._qualifiers


def test_two_instances_do_not_share_nested_qualifier_values(two_freerflags):
    """A shared list would let `.append()` on one instance reach every other."""
    a, b = two_freerflags
    assert a.get_qualifier("enumerators") is not b.get_qualifier("enumerators")


def test_set_qualifier_does_not_reach_a_sibling(two_freerflags):
    a, b = two_freerflags
    a.set_qualifier("guiLabel", "CHANGED ON A")
    assert a.get_qualifier("guiLabel") == "CHANGED ON A"
    assert b.get_qualifier("guiLabel") != "CHANGED ON A"


def test_set_qualifier_does_not_poison_later_instances():
    """The cache outlives the instance, so this is the one that matters."""
    from ccp4i2.core.tasks import get_plugin_class
    klass = get_plugin_class("freerflag")

    before = klass(parent=None, name="before")
    original = before.container.controlParameters.GEN_MODE.get_qualifier("guiLabel")

    victim = klass(parent=None, name="victim")
    victim.container.controlParameters.GEN_MODE.set_qualifier("guiLabel", "POISON")

    after = klass(parent=None, name="after")
    assert after.container.controlParameters.GEN_MODE.get_qualifier("guiLabel") == original


def test_a_declared_qualifier_still_arrives():
    """Isolation must not have been bought by dropping the declaration."""
    from ccp4i2.core.tasks import get_plugin_class
    p = get_plugin_class("freerflag")(parent=None, name="f")
    frac = p.container.controlParameters.FRAC
    assert frac.get_qualifier("max") == 1.0
    assert frac.get_qualifier("min") == 0.0
