"""AREAIMOL's report leads with the numbers read from its log.

The log below is an excerpt of a comparison run (the MDM2 model with and
without Nutlin-3a); the report used to show only the log's summaries, with
the residues Nutlin covers buried halfway down.
"""
from ccp4i2.wrappers.areaimol.areaimol import parse_areas

LOG = """
 TOTAL AREA:     5592.2

 TOTAL AREA:     5704.3

 Differences for individual chains.
 ----------------------------------

 ILE A  19    -0.3    LEU A  54   -28.5    GLY A  58   -26.3
 MET A  62   -12.0    TYR A  67    -1.0    GLN A  72    -0.7
 HIS A  96   -13.3    ILE A  99   -15.0    SO4 A 202   -12.3
<B><FONT COLOR="#FF0000"><!--SUMMARY_BEGIN-->

 Total area difference of chain 'A':     -117.8

 TOTAL AREA DIFFERENCE:     -117.8
"""


def test_totals_and_difference():
    areas = parse_areas(LOG)
    assert areas["totals"] == [5592.2, 5704.3]
    assert areas["difference"] == -117.8


def test_residues_largest_change_first():
    residues = parse_areas(LOG)["residues"]
    assert len(residues) == 9
    assert [(r["name"], r["number"]) for r in residues[:3]] == [
        ("LEU", "54"), ("GLY", "58"), ("ILE", "99")]
    assert residues[0]["chain"] == "A" and residues[0]["change"] == -28.5


def test_single_model_has_no_difference():
    areas = parse_areas(" TOTAL AREA:     5704.3\n")
    assert areas == {"totals": [5704.3], "difference": None, "residues": []}


PROGRAM_XML = """<areaimol><Areas><Total>5592.2</Total><Total>5704.3</Total>
<Difference>-117.8</Difference>
<Residue><name>LEU</name><chain>A</chain><number>54</number><change>-28.5</change></Residue>
</Areas><SASValues/></areaimol>"""


def test_report_leads_with_the_areas():
    # It failed to render at first: an HTML entity (&Aring;) in its text,
    # which the report's XML does not define.
    import xml.etree.ElementTree as ET
    from ccp4i2.wrappers.areaimol.areaimol_report import areaimol_report
    report = areaimol_report(xmlnode=ET.fromstring(PROGRAM_XML), jobInfo={},
                             jobStatus="Finished")
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "117.8 Å² less accessible area" in text
    assert "LEU" in text
