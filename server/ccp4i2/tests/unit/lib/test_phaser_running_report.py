"""The running Phaser report reads its best solution instead of quoting it."""
from xml.etree import ElementTree as ET

from ccp4i2.core.tasks import get_report_class

BEST = ("SOLU SET RFZ=4.6 TFZ=6.2 PAK=0 LLG=40 TFZ==8.1 SOLU SPAC P 21 21 21 "
        "SOLU 6DIM ENSE SearchModel EULER 46.9 112.9 213.6 FRAC -1.28 -0.30 0.62 BFAC 33.06 #TFZ==8.1")


def _render(best):
    root = ET.fromstring(
        "<PHASER><CurrentModule>MOLECULAR REPLACEMENT REFINEMENT AND PHASING</CurrentModule>"
        "<CurrentActivity><label>Generating Statistics</label><max>0</max><value>0</value></CurrentActivity>"
        f"<PhaserCurrentBestSolution>{best}</PhaserCurrentBestSolution></PHASER>")
    report = get_report_class("phaser_mr_auto_phil")(xmlnode=root, jobInfo={"fileroot": "/nowhere"},
                                                     jobStatus="Running")
    return ET.tostring(report.as_data_etree(), encoding="unicode")


def test_the_best_solution_is_summarised_and_tabulated():
    xml = _render(BEST)
    assert "Current best solution: space group P 21 21 21, 1 component placed, refined TFZ 8.1" in xml
    table = ET.fromstring(xml).find(".//*[@key='PhaserCurrentBestTable']")
    assert table is not None
    text = " ".join(table.itertext())
    for cell in ("RFZ", "4.6", "TFZ", "6.2", "Clashes", "0", "LLG", "40"):
        assert cell in text
    # the raw annotation is still there, folded away
    assert "SOLU 6DIM ENSE SearchModel" in xml


def test_two_components_make_two_rows():
    xml = _render("SOLU SET RFZ=10.4 TFZ=10.4 PAK=0 LLG=255 TFZ==12.3 RFZ=5.1 TFZ=18.6 PAK=2 LLG=612 TFZ==20.1 "
                  "SOLU SPAC P 32 2 1")
    assert "2 components placed, refined TFZ 20.1" in xml
    table = ET.fromstring(xml).find(".//*[@key='PhaserCurrentBestTable']")
    text = " ".join(table.itertext())
    assert "18.6" in text and "612" in text
