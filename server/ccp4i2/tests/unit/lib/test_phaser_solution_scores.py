"""Phaser's solution tables show only the scores that were calculated.

Phaser keeps every score on a solution and leaves one it did not calculate
at 0, so a translation-function solution was shown with LLG 0.00, R 0.00
and no clashes, as if they had been measured. The solution's placements
(its annotation) say what was calculated. The XML is a real MR_FTF result
(MDM2, Chainsaw model of MDMX).
"""
import xml.etree.ElementTree as ET

import pytest

pytest.importorskip("lxml")
from ccp4i2.wrappers.phaser_mr_ftf_phil.script.phaser_mr_ftf_phil_report import (  # noqa: E402
    phaser_mr_ftf_phil_report,
)

FTF = """<PhaserReport><Solutions><Solution><Number>1</Number>
<Annotation>RFZ=2.8 TFZ=8.0</Annotation><LLG>0.00</LLG><TFZ>7.96</TFZ><TFZeq>0.00</TFZeq>
<R>0.00</R><PAK>0.00</PAK><spaceGroup>P 65 2 2</spaceGroup><History>RF/TF(23/1:1)</History>
<Placements><Placement><RFZ>2.8</RFZ><TFZ>8.0</TFZ></Placement></Placements>
</Solution></Solutions></PhaserReport>"""


def cells(text, table_id):
    start = text.index(table_id)
    body = text[text.index("<tbody>", start):text.index("</tbody>", start)]
    return [c.text for c in ET.fromstring(body + "</tbody>").iter("td")]


def test_translation_function_shows_no_uncalculated_scores(tmp_path):
    report = phaser_mr_ftf_phil_report(xmlnode=ET.fromstring(FTF), jobStatus="Finished",
                                       jobInfo={"fileroot": str(tmp_path)})
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    # #, space group, LLG, TFZ, TFZ-equiv, R, clashes, annotation
    assert cells(text, "PhaserSolutionsTable") == [
        "1", "P 65 2 2", "–", "7.96", "–", "–", "–", "RFZ=2.8 TFZ=8.0"]
    # RFZ, TFZ, TFZ-equiv, clashes, LLG, note
    assert cells(text, "PhaserPlacements0")[:5] == ["2.8", "8.0", "–", "–", "–"]
    assert "No search attempt recorded" not in text
