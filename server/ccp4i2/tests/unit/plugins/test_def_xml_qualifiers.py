"""Qualifiers a def.xml declares must survive parsing.

Two of them were being dropped, which mattered because the interface builds
folders and expert-level filtering from exactly these:

* A <container> carries <qualifiers> just as <content> does, and that is where
  a folder's guiLabel lives. They were never read, so every folder in a
  def.xml-derived task rendered as an unlabelled grey bar.
* <guiDefinition> holds elements, not a value. Reading its .text yielded the
  whitespace between the tags and discarded <expertLevel>, so no def.xml task
  has ever had a working expert level — with a PHIL-generated task such as
  xia2_dials, that is the difference between 22 parameters and 281.
"""

import pytest

from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser

DEF_XML = """<?xml version='1.0' encoding='ASCII'?>
<ccp4i2>
  <ccp4i2_header>
    <function>DEF</function>
    <pluginName>qualifier_probe</pluginName>
  </ccp4i2_header>
  <ccp4i2_body id="qualifier_probe">
    <container id="controlParameters">
      <container id="scope">
        <qualifiers>
          <guiLabel>Scope settings</guiLabel>
          <guiDefinition/>
        </qualifiers>
        <content id="scope__basic">
          <className>CBoolean</className>
          <qualifiers>
            <guiLabel>Basic thing</guiLabel>
            <guiDefinition><expertLevel>0</expertLevel></guiDefinition>
            <default>False</default>
          </qualifiers>
        </content>
        <content id="scope__expert">
          <className>CBoolean</className>
          <qualifiers>
            <guiLabel>Expert thing</guiLabel>
            <guiDefinition><expertLevel>2</expertLevel></guiDefinition>
            <default>False</default>
          </qualifiers>
        </content>
      </container>
    </container>
  </ccp4i2_body>
</ccp4i2>
"""


@pytest.fixture(name="parsed")
def parsed_fixture(tmp_path):
    def_xml = tmp_path / "qualifier_probe.def.xml"
    def_xml.write_text(DEF_XML)
    return DefXmlParser().parse_def_xml(def_xml)


def test_a_container_keeps_its_gui_label(parsed):
    scope = parsed.controlParameters.scope
    assert scope.get_qualifier("guiLabel") == "Scope settings"


def test_gui_definition_keeps_its_children(parsed):
    scope = parsed.controlParameters.scope
    basic = scope.scope__basic.get_qualifier("guiDefinition")
    expert = scope.scope__expert.get_qualifier("guiDefinition")

    assert basic == {"expertLevel": 0}
    # An int, not "2" — the interface compares it against a threshold.
    assert expert == {"expertLevel": 2}


def test_an_empty_gui_definition_is_harmless(parsed):
    """<guiDefinition/> has no children and must not become a stray dict."""
    scope = parsed.controlParameters.scope
    assert not scope.get_qualifier("guiDefinition")
