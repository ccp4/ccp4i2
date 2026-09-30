"""PAIREF's report shows what paired refinement decided.

It used to show only a link to PAIREF's own HTML page. The files below are
excerpts of a real run (MDM2, data cut at 1.35 A by Aimless, two 0.05 A
shells added).
"""
from ccp4i2.wrappers.pairef.pairef_report import read_pairef_results

R_VALUES = """\
# Shell      Rwork(init) Rwork(fin) Rwork(diff)   Rfree(init) Rfree(fin) Rfree(diff)
1.35A->1.30A      0.2350     0.2317     -0.0033        0.2467     0.2422     -0.0045
1.30A->1.25A      0.2328     0.2327     -0.0001        0.2429     0.2425     -0.0004
"""

LOG = """\
Suggested cutoff:
1.25 A
Calculation ended.
These warning messages appeared during calculation:
WARNING: Intensities were not found in hklin.mtz. CCwork and CCfree values cannot be calculated.
WARNING: There are only 27 < 50 free reflections in the resolution shell 1.30-1.25 A.
Results are listed in logfile PAIREF_project.html
"""


def write_job(tmp_path):
    project = tmp_path / "pairef_project"
    project.mkdir()
    (project / "project_R-values.csv").write_text(R_VALUES)
    (project / "PAIREF_cutoff.txt").write_text("1.25")
    (tmp_path / "log.txt").write_text(LOG)
    return tmp_path


def test_cutoff_and_shells(tmp_path):
    results = read_pairef_results(write_job(tmp_path))
    assert results["cutoff"] == "1.25"
    assert [(s["from"], s["to"]) for s in results["shells"]] == [
        ("1.35", "1.30"), ("1.30", "1.25")]
    assert results["shells"][0]["rfree"] == (0.2467, 0.2422, -0.0045)


def test_warnings(tmp_path):
    warnings = read_pairef_results(write_job(tmp_path))["warnings"]
    assert len(warnings) == 2
    assert warnings[1].startswith("There are only 27 < 50 free reflections")


def test_nothing_yet(tmp_path):
    assert read_pairef_results(tmp_path) == {"cutoff": None, "shells": [], "warnings": []}


def test_report_renders(tmp_path):
    # It failed to render at first: HTML entities (&Aring;, &rarr;) in its
    # text, which the report's XML does not define.
    import xml.etree.ElementTree as ET
    from ccp4i2.wrappers.pairef.pairef_report import pairef_report
    job = write_job(tmp_path)
    report = pairef_report(xmlnode=ET.Element("pairef"), jobStatus="Finished",
                           jobInfo={"fileroot": str(job), "projectid": "p", "jobnumber": "1"})
    text = ET.tostring(report.as_data_etree(), encoding="unicode")
    assert "Suggested high-resolution cutoff: 1.25 Å" in text
    assert "1.35 → 1.30" in text and "-0.0045" in text
    assert "only 27 &lt; 50 free reflections" in text
