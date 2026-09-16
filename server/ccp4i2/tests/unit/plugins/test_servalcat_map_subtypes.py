"""servalcat's map slots require the right subtype.

The half-map inputs must require half maps (subType 5) and the mask input a mask
(subType 4). The old "5,1,0" list silently disabled subtype filtering entirely
(a list containing 0 matches everything in get_by_context), so autopopulation
grabbed the first map it saw -- a trimmed map -- for both half-map slots, and
none of them a half map or a mask.
"""

import xml.etree.ElementTree as ET

from ccp4i2.core.tasks import locate_def_xml


def _local(tag):
    return tag.rsplit("}", 1)[-1]


def _required_subtypes(content_id):
    root = ET.parse(locate_def_xml("servalcat")).getroot()
    for el in root.iter():
        if _local(el.tag) == "content" and el.get("id") == content_id:
            for sub in el.iter():
                if _local(sub.tag) == "requiredSubType":
                    return {int(x) for x in sub.text.split(",")}
            return set()   # slot present but declares no requiredSubType
    raise AssertionError(f"no content id={content_id} in servalcat.def.xml")


def test_half_map_slots_require_half_maps():
    assert _required_subtypes("MAPIN1") == {5}
    assert _required_subtypes("MAPIN2") == {5}


def test_mask_slot_requires_a_mask():
    assert _required_subtypes("MAPMASK") == {4}
