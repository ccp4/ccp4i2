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


# From the Lck kinase run (AlphaFold model, two splits, against 4c3f): the
# C-lobe placed at TFZ 26.4, the N-lobe nowhere (best TFZ 6.0, LLG +18).
# SliceNDice reported TFZ 4.4 and LLG 1738, the values of Phaser's last
# listed solution, not its best.
PHASER_LOG = """
   Solution #1 annotation (history):
   SOLU SET  RFZ=16.5 TFZ=26.4 PAK=1 LLG=482 TFZ==28.3 RFZ=3.1 TFZ=6.0 PAK=19 LLG=500 TFZ==6.1 LLG=1745 TFZ==6.4 PAK=20
   SOLU SPAC P 21 21 21
   SOLU 6DIM ENSE pdb_AF-P06239-F1-model_v6_kinase_cluster_0 EULER  334.5   73.0  133.5 FRAC  0.43 -0.33 -0.15 BFAC
   SOLU 6DIM ENSE pdb_AF-P06239-F1-model_v6_kinase_cluster_1 EULER  141.2    7.3  350.1 FRAC  0.54 -0.12 -0.09 BFAC
   Solution #2 annotation (history):
   SOLU SET  RFZ=16.5 TFZ=26.4 PAK=1 LLG=482 TFZ==28.3 RFZ=3.1 TFZ=4.4 PAK=34 LLG=495 LLG=1738 PAK=33 LLG=1738
   SOLU 6DIM ENSE pdb_AF-P06239-F1-model_v6_kinase_cluster_0 EULER  334.5   73.0  133.5 FRAC  0.43 -0.33 -0.15 BFAC
   SOLU 6DIM ENSE pdb_AF-P06239-F1-model_v6_kinase_cluster_1 EULER  185.1   86.7  306.5 FRAC  0.25 -0.39 -0.04 BFAC
"""


def test_phaser_components_of_the_best_solution():
    assert sd.phaser_components(PHASER_LOG) == [("0", "26.4", "482", "1"), ("1", "6.0", "500", "19")]
    assert sd.phaser_components("no solutions here") == []


def test_report_names_the_piece_not_placed(tmp_path):
    xml = ET.fromstring("""<SliceNDice><RunInfo>
<Best><bid>2</bid><R>0.4062</R><RFree>0.4274</RFree><Solved>True</Solved><Partial>True</Partial></Best>
<Split id="2"><Model>288-302, 320-506</Model><Model>225-287, 303-319</Model></Split>
<Sol><SolID>2</SolID><llg>1738.0</llg><tfz>4.4</tfz><srf>0.4062</srf><sre>0.4274</sre>
<Component cluster="0" tfz="26.4" llg="482" clashes="1"/>
<Component cluster="1" tfz="6.0" llg="500" clashes="19"/></Sol>
</RunInfo></SliceNDice>""")
    text = report_text(xml, tmp_path)
    assert "Partly solved" in text
    assert "residues 288-302, 320-506: TFZ 26.4, LLG +482" in text
    assert "residues 225-287, 303-319: TFZ 6.0, LLG +18, 19 clashes" in text
    assert "residues 225-287, 303-319 probably did not find its place" in text
