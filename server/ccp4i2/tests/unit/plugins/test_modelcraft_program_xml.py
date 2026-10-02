"""modelcraft writes program.xml from ModelCraft's own report.

It wrote none: what a build achieved (residues, waters, R factors per
cycle, why it stopped) was only in modelcraft/modelcraft.json, which nothing
reading program.xml, an agent judging the job included, could see.
"""
import pytest

pytest.importorskip("lxml")

from ccp4i2.wrappers.modelcraft.script.modelcraft import program_xml  # noqa: E402

# Trimmed from the Thaumatin build in the docs scenario (ModelCraft 6.1.1).
RESULT = {
    "version": "6.1.1",
    "termination_reason": "Normal",
    "args": ["xray", "--contents", "contents.json", "--model", "xyzin.cif"],
    "jobs": [
        {"name": "csheetbend", "seconds": 2.1},
        {"name": "refmacat", "rwork": 0.27, "rfree": 0.30, "initial_rwork": 0.32,
         "initial_rfree": 0.33, "resolution_high": 1.17, "data_completeness": 97.5},
        {"name": "refmacat", "rwork": 0.21, "rfree": 0.24, "initial_rwork": 0.25,
         "initial_rfree": 0.26, "resolution_high": 1.17, "data_completeness": 97.5},
    ],
    "cycles": [
        {"cycle": 1, "residues": 208, "protein": 208, "nucleic": 0, "waters": 120,
         "r_work": 0.21, "r_free": 0.238},
        {"cycle": 2, "residues": 207, "protein": 207, "nucleic": 0, "waters": 145,
         "r_work": 0.204, "r_free": 0.234},
    ],
    "final": {"cycle": 2, "residues": 207, "protein": 207, "nucleic": 0, "waters": 145,
              "r_work": 0.204, "r_free": 0.234},
}


def test_the_build_is_recorded():
    root = program_xml(RESULT)
    assert root.findtext("TerminationReason") == "Normal"
    assert root.findtext("Final/r_free") == "0.234"
    assert root.findtext("Final/residues") == "207"
    assert [c.findtext("r_free") for c in root.findall("Cycles/Cycle")] == ["0.238", "0.234"]
    # The model it started from: its first refinement's starting R factors.
    assert root.findtext("InputModel/r_free") == "0.33"
    assert root.findtext("ResolutionHigh") == "1.17"


def test_an_early_stop_says_why_and_claims_nothing():
    root = program_xml({"termination_reason": "No residues built", "cycles": []})
    assert root.findtext("TerminationReason") == "No residues built"
    assert root.find("Final") is None and root.find("InputModel") is None


def test_from_phases_alone_there_is_no_input_model():
    # The first refinement is then of the first build (GammaXe, from
    # experimental phases): its starting R is not an input model's.
    result = dict(RESULT, args=["xray", "--contents", "contents.json", "--phases", "x.mtz"])
    root = program_xml(result)
    assert root.find("InputModel") is None
    assert root.findtext("ResolutionHigh") == "1.17"
