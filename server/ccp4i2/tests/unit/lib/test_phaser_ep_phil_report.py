"""The Phaser EP pipeline's report draws its SHELXC/D section.

It imported the SHELX report from a module that does not exist and called a
method it does not have, inside a try/except, so the section only ever read
"SHELXC/D ran (report unavailable)". No SHELX is needed to check it: the
pipeline files the SHELXC/D record under a ShelxCD node, and the report has
only to reach the class that draws it.
"""
import xml.etree.ElementTree as ET

from lxml import etree

from ccp4i2.pipelines.phaser_ep_phil.script.phaser_ep_phil_report import (
    phaser_ep_phil_report,
)


def test_the_shelx_section_is_drawn_not_excused():
    root = etree.fromstring(
        "<PhaserEpPipeline><ShelxCD><Shelxc/></ShelxCD></PhaserEpPipeline>")
    report = phaser_ep_phil_report(xmlnode=root, jobStatus="Finished", jobInfo={})
    xml = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "Substructure search (SHELXC/D)" in xml
    assert "report unavailable" not in xml
