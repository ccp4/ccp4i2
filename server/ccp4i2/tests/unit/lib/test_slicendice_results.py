"""SliceNDice's results are judged, not just ranked.

The wrapper reported the lowest R-free as "the best MR solution" whatever it
was. On MDM2 with its AlphaFold model, the one placement SliceNDice tried
(the whole trimmed model, one split) refined to R-free 0.555: no solution,
reported as one. And SliceNDice 0.1.3 runs MR on one split of those it makes
(a loop-indentation slip in its command_line/slice.py), which the report now
says rather than hiding.
"""
import xml.etree.ElementTree as ET

import pytest

sd = pytest.importorskip("ccp4i2.wrappers.slicendice.script.slicendice")

RESULTS = {
    "slice": {
        "split_1": {"residues_ranges": {"0": [["M1:A:26", "M1:A:490"]]}},
        "split_2": {"residues_ranges": {"0": [["M1:A:247", "M1:A:490"]],
                                        "1": [["M1:A:26", "M1:A:111"]]}},
    },
    "dice": {"split_1": {"phaser_llg": 14.0, "phaser_tfz": 3.6,
                         "final_r_fact": 0.571, "final_r_free": 0.5552}},
}


def test_split_ranges():
    assert sd.split_ranges(RESULTS) == {"split_1": ["26-490"],
                                        "split_2": ["247-490", "26-111"]}


def test_a_high_r_free_is_not_a_solution():
    assert sd.best_solution(RESULTS) == ("split_1", False)
    solved = {"dice": {"split_2": {"final_r_fact": 0.30, "final_r_free": 0.34},
                       "split_1": {"final_r_fact": 0.57, "final_r_free": 0.55}}}
    assert sd.best_solution(solved) == ("split_2", True)
    assert sd.best_solution({"dice": {}}) == (None, False)


PROGRAM_XML = """<SliceNDice><RunInfo>
<Best><bid>1</bid><R>0.571</R><RFree>0.5552</RFree><Solved>False</Solved></Best>
<Split id="1"><Model>26-490</Model></Split>
<Split id="2"><Model>247-490</Model><Model>26-111</Model></Split>
<Sol><SolID>1</SolID><llg>14.0</llg><tfz>3.6</tfz><srf>0.571</srf><sre>0.5552</sre></Sol>
</RunInfo></SliceNDice>"""


def report_text(xml, tmp_path):
    from ccp4i2.wrappers.slicendice.script.slicendice_report import slicendice_report
    report = slicendice_report(xmlnode=xml, jobInfo={"fileroot": str(tmp_path)},
                               jobStatus="Finished")
    return ET.tostring(report.as_data_etree(), encoding="unicode")


def test_report_says_no_solution_and_what_was_not_tried(tmp_path):
    log = tmp_path / "slicendice_0"
    log.mkdir()
    (log / "slicendice.log").write_text("R-free < 0.45 & done\n")
    text = report_text(ET.fromstring(PROGRAM_XML), tmp_path)
    assert "No solution" in text and "0.555" in text
    assert "Not tried in molecular replacement: 2 split" in text
    assert "26-111" in text
    assert "R-free &lt; 0.45 &amp; done" in text   # the log, escaped


def test_running_report_needs_no_results(tmp_path):
    text = report_text(None, tmp_path)
    assert "Running" in text
