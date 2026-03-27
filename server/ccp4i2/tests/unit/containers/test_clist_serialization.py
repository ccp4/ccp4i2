"""Test CAsuContentSeqList serialization to verify XML export works correctly."""

import xml.etree.ElementTree as ET

from ccp4i2.core.CCP4ModelData import CAsuContentSeq, CAsuContentSeqList


def test_clist_serialization():
    """Test that CAsuContentSeqList properly serializes populated items."""
    seq_list = CAsuContentSeqList(name="ASU_CONTENT")

    item1 = seq_list.makeItem()
    item1.sequence = "HPETLVKVK"
    item1.nCopies = 1
    item1.name = "BETA"
    item1.polymerType = "PROTEIN"
    item1.description = "Beta-lactamase"
    seq_list.append(item1)

    item2 = seq_list.makeItem()
    item2.sequence = "AGVMTGAKFT"
    item2.nCopies = 1
    item2.name = "BLIP"
    item2.polymerType = "PROTEIN"
    item2.description = "Beta-lactamase inhibitory protein"
    seq_list.append(item2)

    assert len(seq_list) == 2

    etree = seq_list.getEtree()
    items = etree.findall('.//CAsuContentSeq')
    assert len(items) == 2

    for i, item_elem in enumerate(items):
        children = list(item_elem)
        assert len(children) > 0, f"Item {i} has no child elements in XML"

    xml_str = ET.tostring(etree, encoding='unicode')
    assert '<name>BETA</name>' in xml_str
    assert '<name>BLIP</name>' in xml_str
