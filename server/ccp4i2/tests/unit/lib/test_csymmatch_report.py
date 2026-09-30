"""The csymmatch report says what csymmatch did.

Its program XML is <Csymmatch> itself, and the report looked for
".//Csymmatch/..." below it, so it found nothing and showed only the file
lists: no change of origin, no operators, no scores.
"""
import xml.etree.ElementTree as ET

from ccp4i2.pipelines.phaser_pipeline.script.phaser_pipeline_report import (
    csymmatch_report as registered,
)
from ccp4i2.wrappers.csymmatch.script.csymmatch_report import csymmatch_report

PROGRAM_XML = """<Csymmatch>
  <ChangeOfHand> NO</ChangeOfHand>
  <ChangeOfOrigin> uvw = (0, 0.5, 0)</ChangeOfOrigin>
  <Segment>
    <Range>Chain: A  703- 822 </Range>
    <Operator> x+1/2, -y+1/2, -z</Operator>
    <Shift>     uvw = (0, 0, 0)</Shift>
    <Score>  0.482834</Score>
  </Segment>
</Csymmatch>"""


def test_the_registered_report_is_the_wrapper_one():
    assert registered is csymmatch_report


def test_report_describes_the_origin_shift_and_segments():
    report = csymmatch_report(xmlnode=ET.fromstring(PROGRAM_XML), jobInfo={})
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "A change of origin was applied" in text
    assert "grouped into 1 segments" in text
    assert "x+1/2, -y+1/2, -z" in text
    assert "change of hand was applied" not in text
