"""Difference-map outputs declare subType 2 (issue #524, F2 curation).

subType is unknowable from a real-space map's bytes, so a difference map must be
*declared* as one at its producing parameter, else it is gleaned untyped and
silently read as a normal map. These are the non-default outputs curated per the
#524 plan: servalcat's SPA Fo-Fc map and SubtractNative's isomorphous-difference
map. (Normal/default maps stay bare; the default already types them.)
"""

import xml.etree.ElementTree as ET

import pytest

from ccp4i2.core.tasks import locate_def_xml


def _local(tag):
    return tag.rsplit("}", 1)[-1]


def _output_subtype(task, content_id):
    path = locate_def_xml(task)
    assert path, f"{task} not registered / def.xml not found"
    root = ET.parse(path).getroot()
    for el in root.iter():
        if _local(el.tag) == "content" and el.get("id") == content_id:
            for sub in el.iter():
                if _local(sub.tag) == "subType" and (sub.text or "").strip():
                    return int(sub.text)
            return None  # declared, but no subType
    raise AssertionError(f"no output content id={content_id} in {task}.def.xml")


@pytest.mark.parametrize("task,content_id", [
    ("servalcat", "MAP_FOFC"),
    ("SubtractNative", "MAPOUT"),
])
def test_difference_map_output_declares_subtype_2(task, content_id):
    assert _output_subtype(task, content_id) == 2


def test_servalcat_fo_map_stays_untyped():
    # The plain Fo density map is the default kind -> deliberately bare.
    assert _output_subtype("servalcat", "MAP_FO") is None
