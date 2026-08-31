"""What a container is *declared* to contain bounds what it serialises.

`__setattr__` appends to `_data_order` for any CData assigned to a
non-underscore name, so before this anything dropped on a container joined its
model and was written to `params.xml` under its own `_name` --- not even the
attribute name it was assigned to:

    c.casual_cdata = CString('dropped')
    c.dataOrder()  -> [..., 'casual_cdata']
    getEtree()     -> <CString_4443846000>dropped</CString_4443846000>

`_declared_content` records only what an installation path put there. Plain
assignment does not touch it, so it cannot silently extend the model.

The object remains a perfectly usable attribute --- this bounds the *model*,
not what a plugin may hold on itself.

Gating the *assignment* rather than the serialisation was tried and is not
available: 9,443 of 9,445 CData assignments during construction happen while
the name is still undeclared, because every installation path calls setattr
first and declares afterwards.

Pure Python -- no CCP4 binaries needed.
"""
import xml.etree.ElementTree as ET

import pytest

django = pytest.importorskip("django")


@pytest.fixture
def inputs():
    """Yields the container *and* keeps the plugin alive.

    Returning `plugin.container.inputData` and letting the plugin go is enough
    to lose the children: the plugin's destroy() clears them, and a
    sub-container handed out on its own does not keep its parent alive. That
    is long-standing behaviour, identical before and after the child-structure
    collapse --- but it makes a fixture that discards the plugin quietly test
    an empty container.
    """
    from ccp4i2.core.tasks import get_plugin_class
    plugin = get_plugin_class("freerflag")(parent=None, name="f")
    yield plugin.container.inputData
    del plugin


def test_a_dropped_object_does_not_join_the_model(inputs):
    from ccp4i2.core.base_object.fundamental_types import CString

    before = list(inputs.dataOrder())
    inputs.casual_cdata = CString("dropped")
    assert list(inputs.dataOrder()) == before


def test_a_dropped_object_is_not_serialised(inputs):
    from ccp4i2.core.base_object.fundamental_types import CString

    inputs.casual_cdata = CString("dropped")
    written = ET.tostring(inputs.getEtree("inputData", excludeUnset=True),
                          encoding="unicode")
    assert "casual_cdata" not in written
    assert "CString_" not in written, \
        "serialised under its own _name, which is the older and worse symptom"


def test_it_is_still_an_ordinary_attribute(inputs):
    """Bounding the model must not stop a plugin holding things on itself."""
    from ccp4i2.core.base_object.fundamental_types import CString

    inputs.casual_cdata = CString("dropped")
    assert inputs.casual_cdata.value == "dropped"


def test_declared_content_still_serialises(inputs):
    """The other half: declaring must actually work."""
    inputs.F_SIGF.setFullPath("/tmp/does-not-need-to-exist.mtz")
    written = ET.tostring(inputs.getEtree("inputData", excludeUnset=True),
                          encoding="unicode")
    assert "F_SIGF" in written


def test_every_container_declares_what_it_reports():
    """The equality this rests on, across the registry.

    If the two ever diverge it is an assignment nobody declared, rather than a
    gap in the bookkeeping --- which is the whole point of keeping the list.
    """
    from ccp4i2.core.tasks import TASKS, get_plugin_class
    from ccp4i2.core.base_object.ccontainer import CContainer

    checked, undeclared = 0, []
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
            if isinstance(obj, CContainer):
                checked += 1
                missing = set(obj.dataOrder() or []) - set(obj.declaredContent())
                if missing:
                    undeclared.append(f"{name}.{path}: {sorted(missing)}")
            for child in obj.children():
                walk(child, f"{path}.{child._name}", depth + 1)

        walk(plugin.container, "container")

    assert checked > 500, f"only {checked} containers checked --- proves little"
    assert not undeclared, f"reported but never declared: {undeclared[:6]}"
