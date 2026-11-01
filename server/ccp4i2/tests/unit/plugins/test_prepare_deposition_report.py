"""The sequence coverage track says where the unmodelled residues are.

The progress bar it replaces showed matches over length, which for a chain
missing its tag, an internal loop and two C-terminal residues read as one
block missing at the C-terminus.
"""
import importlib

from ccp4i2.core.base_object.fundamental_types import CList, CString
from ccp4i2.core.tasks import TASKS, locate_def_xml
from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser
from ccp4i2.wrappers.modelASUCheck.script.modelASUCheck_report import (
    coverage_from_alignment,
    coverage_track_svg,
    describe_runs,
)

# Construct of 12; the model lacks 1-2, 6 and 11-12, and has a different
# residue at 8. (gemmi's match string: | match, . mismatch, space gap.)
ALIGN_SEQ = "MGSAKLVADERK"
ALIGN_MOD = "--SAK-VYDE--"
MATCH = "  |||.|.||  "


def test_coverage_finds_every_gap_not_just_the_end():
    length, modelled, mismatches, runs = coverage_from_alignment(ALIGN_SEQ, ALIGN_MOD, MATCH)
    assert length == 12
    assert runs == [(1, 2), (6, 6), (11, 12)]
    assert describe_runs(runs) == "1-2, 6, 11-12"
    assert mismatches == {8}
    assert modelled == {3, 4, 5, 7, 9, 10}


def test_model_residues_the_construct_lacks_are_not_on_the_axis():
    # An insertion in the model (construct gap) must not shift the numbering.
    length, modelled, _m, runs = coverage_from_alignment("MG-SAK", "MGXSAK", "|| |||")
    assert length == 5 and modelled == {1, 2, 3, 4, 5} and runs == []


def test_track_is_drawn_as_runs():
    svg = coverage_track_svg(12, {3, 4, 5, 7, 9, 10}, {8})
    assert svg.startswith("<svg") and svg.endswith("</svg>")
    assert svg.count('fill="#4caf50"') == 3      # 3-5, 7, 9-10
    assert svg.count('fill="#ff9800"') == 1
    assert coverage_track_svg(0, set(), set()) == ""


def test_the_task_registers_its_own_types():
    assert TASKS["prepare_deposition"].dataTypes == (
        "ccp4i2.wrappers.prepare_deposition.script.prepare_deposition_types",)
    assert TASKS["adding_stats_to_mmcif_i2"].successor == "prepare_deposition"
    types = importlib.import_module(TASKS["prepare_deposition"].dataTypes[0])
    container = DefXmlParser().parse_def_xml(locate_def_xml("prepare_deposition"))
    for path, item_class in (("controlParameters.AUTHORS", "CDepositionAuthor"),
                             ("controlParameters.RELATED", "CRelatedEntry"),
                             ("controlParameters.ENTITIES", "CEntityMetadata"),
                             ("controlParameters.SOFTWARE", "CSoftwareEntry")):
        obj = container
        for part in path.split("."):
            obj = getattr(obj, part)
        assert isinstance(obj, CList), path
        sub = obj.get_qualifier("subItem")
        assert sub and sub["class"].__name__ == item_class, (path, sub)
        assert sub["class"] is not CString
    # The nested list of segments resolves too (the silent failure is a CString).
    entity = types.CEntityMetadata()
    assert isinstance(entity.SEGMENTS, CList)
    assert entity.SEGMENTS.get_qualifier("subItem")["class"] is types.CSourceSegment
