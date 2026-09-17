"""ImportUnmerged refuses a merged file offered as unmerged.

Diamond's auto-processing publishes ``scaled.mtz`` beside ``scaled_unmerged.mtz``
and the names differ by one word, so handing over the merged one is the routine
mistake. Without this check it is not caught at the Run dialog but much later,
inside aimless, as a complaint about missing batches.

The refusal is raised by ``runTimeValidity`` rather than ``process()`` so it
reaches the user before the job is submitted.

CCP4-free: gemmi reads the demo MTZ, no CCP4 binaries.
"""

import pytest

from ccp4i2.core import CCP4ErrorHandling
from ccp4i2.core.CCP4Utils import getCCP4I2Dir
from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser
from ccp4i2.core.tasks import locate_def_xml
from ccp4i2.wrappers.ImportUnmerged.script.ImportUnmerged import (
    ImportUnmerged,
    _looks_like_space_group,
)

pytest.importorskip("gemmi", reason="needs gemmi to read the demo MTZ")

DEMO = getCCP4I2Dir() + "/demo_data/gamma/"
MERGED = DEMO + "merged_intensities_native.mtz"
UNMERGED = DEMO + "HKLOUT_unmerged.mtz"


def _plugin(input_path):
    container = DefXmlParser().parse_def_xml(locate_def_xml("ImportUnmerged"))
    container.inputData.UNMERGEDIN.setFullPath(input_path)
    plugin = object.__new__(ImportUnmerged)
    plugin.container = container
    return plugin


def test_merged_mtz_is_refused_naming_the_input():
    """Code 203 against UNMERGEDIN, so the Run dialog marks the right field."""
    plugin = _plugin(MERGED)
    message = plugin.validate_source(MERGED)
    assert message is not None
    assert "merged" in message.lower()
    assert "ImportObs" in message, "the refusal should say where to go instead"

    error = ImportUnmerged.runTimeValidity(plugin)
    entries = [e for e in error.getErrors() if e.get("code") == 203]
    assert entries, [e.get("code") for e in error.getErrors()]
    assert entries[0].get("name", "").endswith("inputData.UNMERGEDIN"), entries[0]
    assert error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR


def test_genuinely_unmerged_mtz_is_accepted():
    """The guard must not reject the files the task exists to import."""
    plugin = _plugin(UNMERGED)
    assert plugin.validate_source(UNMERGED) is None

    error = ImportUnmerged.runTimeValidity(plugin)
    codes = [e.get("code") for e in error.getErrors()]
    assert 203 not in codes, codes


def test_annotation_carries_the_files_own_facts():
    """The project file list should say what the data is, not just its name."""
    container = DefXmlParser().parse_def_xml(locate_def_xml("ImportUnmerged"))
    out = container.outputData.UNMERGEDOUT
    out.setFullPath(UNMERGED)
    out.loadFile()
    out.annotation.set("Imported unmerged data HKLOUT_unmerged.mtz")

    plugin = object.__new__(ImportUnmerged)
    plugin.container = container
    plugin.finalize_output(out)

    annotation = str(out.annotation)
    assert "MTZ" in annotation
    assert "P 21 21 21" in annotation, annotation


def test_a_mis_parsed_space_group_is_left_out_of_the_annotation():
    """The scalepack reader takes the space group from a fixed position without
    checking it, so a malformed file can yield arbitrary text where a symbol
    belongs. Show the format alone rather than presenting that as fact. (The
    XDS-read-as-scalepack case that exposed this is now fixed in the loader;
    this is the backstop for the positional parse itself.)"""
    assert _looks_like_space_group("P 21 21 21")
    assert _looks_like_space_group("C2221")
    assert _looks_like_space_group("P 43 21 2")
    assert not _looks_like_space_group("MERGE=FALSE    FRIEDEL'S_LAW=FALSE")
    assert not _looks_like_space_group("")
