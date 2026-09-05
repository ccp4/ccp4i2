"""An unset primitive must be distinguishable from one set to its sentinel.

CInt/CFloat/CString/CBoolean hold 0 / 0.0 / "" / False when nothing has been
set, so the JSON the GUI reads must carry the set-state alongside _value, and
the GUI must be able to send null to return a field to unset.
"""
import json

from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core.base_object.fundamental_types import CBoolean, CFloat, CInt, CString
from ccp4i2.lib.utils.containers.json_encoder import CCP4i2JsonEncoder


def encode(obj):
    return json.loads(json.dumps(obj, cls=CCP4i2JsonEncoder))


def test_unset_float_reports_not_set_with_sentinel_value():
    encoded = encode(CFloat())
    assert encoded["_value"] == 0.0
    assert encoded["_valueState"] == "NOT_SET"


def test_default_and_explicit_states_are_reported():
    x = CFloat()
    x.setDefault(0.05)
    assert encode(x)["_valueState"] == "DEFAULT"
    x.set(0.3)
    assert encode(x)["_valueState"] == "EXPLICITLY_SET"


def test_non_primitives_do_not_carry_a_value_state():
    c = CContainer(name="c")
    c.addContent(CInt, "N")
    encoded = encode(c)
    assert "_valueState" not in encoded
    assert encoded["_value"]["N"]["_valueState"] == "NOT_SET"


def test_set_none_unsets_and_keeps_the_type_sentinel():
    for cls, sentinel, value in (
        (CInt, 0, 7), (CFloat, 0.0, 0.3), (CString, "", "abc"), (CBoolean, False, True),
    ):
        x = cls()
        x.set(value)
        assert x.isSet()
        x.set(None)
        assert not x.isSet()
        assert x._value == sentinel, cls.__name__
        assert str(x) == str(sentinel), cls.__name__
        assert encode(x)["_value"] == sentinel


def test_set_parameter_with_none_returns_primitive_to_unset():
    c = CContainer(name="root")
    c.addContent(CFloat, "FRAC")
    c.set_parameter("root.FRAC", 0.3)
    assert c.FRAC.isSet() and float(c.FRAC) == 0.3
    c.set_parameter("root.FRAC", None)
    assert not c.FRAC.isSet()
    assert float(c.FRAC) == 0.0
    assert encode(c.FRAC)["_valueState"] == "NOT_SET"


DEF_XML = """<?xml version='1.0' encoding='ASCII'?>
<ccp4i2>
  <ccp4i2_header>
    <function>DEF</function>
    <pluginName>default_probe</pluginName>
  </ccp4i2_header>
  <ccp4i2_body id="default_probe">
    <container id="controlParameters">
      <content id="FRAC">
        <className>CFloat</className>
        <qualifiers>
          <default>0.05</default>
          <min>0.0</min>
          <max>1.0</max>
        </qualifiers>
      </content>
    </container>
  </ccp4i2_body>
</ccp4i2>
"""


def test_def_xml_default_survives_clearing_the_field(tmp_path):
    """A def.xml <default> is applied to the value, not kept as a qualifier,
    so it must be retained elsewhere for the GUI to show once cleared."""
    from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser

    def_xml = tmp_path / "default_probe.def.xml"
    def_xml.write_text(DEF_XML)
    container = DefXmlParser().parse_def_xml(def_xml)
    frac = container.controlParameters.FRAC

    encoded = encode(frac)
    assert encoded["_valueState"] == "DEFAULT"
    assert encoded["_value"] == 0.05
    assert encoded["_qualifiers"]["default"] == 0.05

    frac.set(None)
    encoded = encode(frac)
    assert encoded["_valueState"] == "NOT_SET"
    assert encoded["_value"] == 0.0
    assert encoded["_qualifiers"]["default"] == 0.05
