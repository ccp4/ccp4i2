"""Half maps are a first-class ``CMapDataFile`` sub-type.

A cryo-EM half map is a distinct *kind* of map (one of a pair, for FSC
cross-validation), not a convertible representation -- so it is a ``subType``,
alongside normal/difference/anom/mask. These tests pin the constant and the
def.xml wiring that makes producers tag half maps and consumers recognise them.

CCP4-free: only the core data classes and the def.xml parser, no binaries.
"""

from ccp4i2.core.CCP4XtalData import CMapDataFile
from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser
from ccp4i2.core.tasks import locate_def_xml


def test_halfmap_subtype_constant_is_distinct():
    assert CMapDataFile.SUBTYPE_HALFMAP == 5
    values = {
        CMapDataFile.SUBTYPE_NORMAL,
        CMapDataFile.SUBTYPE_DIFFERENCE,
        CMapDataFile.SUBTYPE_ANOM_DIFFERENCE,
        CMapDataFile.SUBTYPE_MASK,
        CMapDataFile.SUBTYPE_HALFMAP,
    }
    assert len(values) == 5  # all distinct


def test_molrep_map_tags_half_map_outputs_and_inputs():
    c = DefXmlParser().parse_def_xml(locate_def_xml("molrep_map"))
    # Outputs: half maps carry the half-map sub-type, masks carry the mask one.
    assert c.outputData.HALFMAPOUT1.subType == CMapDataFile.SUBTYPE_HALFMAP
    assert c.outputData.HALFMAPOUT2.subType == CMapDataFile.SUBTYPE_HALFMAP
    assert c.outputData.ORIGINALMASK.subType == CMapDataFile.SUBTYPE_MASK
    assert c.outputData.INVERTEDMASK.subType == CMapDataFile.SUBTYPE_MASK
    assert c.outputData.ORIGINALTRIMMEDMAP.subType == CMapDataFile.SUBTYPE_NORMAL
    # Inputs accept half maps, and (for back-compat) normal/untagged maps too.
    for name in ("HALFMAP1", "HALFMAP2"):
        rq = getattr(c.inputData, name).qualifiers("requiredSubType")
        assert CMapDataFile.SUBTYPE_HALFMAP in rq
        assert CMapDataFile.SUBTYPE_NORMAL in rq and 0 in rq


def _subtype_set(rq):
    return set(rq) if isinstance(rq, (list, tuple)) else {int(rq)}


def test_servalcat_map_inputs_require_the_right_subtype():
    c = DefXmlParser().parse_def_xml(locate_def_xml("servalcat"))
    # Half-map slots require *strictly* half maps. Admitting 1/0 (as the slots
    # used to) also admits trimmed and normal maps -- a trimmed map is itself
    # subType 1 -- and, because a requiredSubType list containing 0 disables
    # subtype filtering entirely in get_by_context, admits *every* map: that is
    # how autopopulation grabbed a trimmed map for both half-map slots. There is
    # no way to accept an old subType-1 "half map" without also accepting
    # subType-1 non-half-maps, so the slot requires 5 and files are captured with
    # the right subtype on import instead (see upload_param).
    for name in ("MAPIN1", "MAPIN2"):
        rq = getattr(c.inputData, name).qualifiers("requiredSubType")
        assert _subtype_set(rq) == {CMapDataFile.SUBTYPE_HALFMAP}
    # The mask slot requires a mask, not any map.
    rqm = c.inputData.MAPMASK.qualifiers("requiredSubType")
    assert _subtype_set(rqm) == {CMapDataFile.SUBTYPE_MASK}
