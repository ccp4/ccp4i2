"""A report's select="..." path must not compare numbers.

Report tables and text take their rows with ElementTree's findall(), which
does not evaluate comparisons: select=".//Residue[ZDpa>3.0]" returns nothing,
silently. EDSTATS's four outlier tables (waters and ligands in suspicious
density) were empty for that reason in every job, under headings promising a
list. Select such rows in Python and pass them as selectNodes.
"""
import re
import xml.etree.ElementTree as ET
from pathlib import Path

import ccp4i2

ROOT = Path(ccp4i2.__file__).parent
COMPARISON = re.compile(r"""select\s*=\s*["'][^"']*\[[^\]"']*[<>][^\]"']*\]""")


def test_findall_does_not_compare():
    """The premise: if this ever starts matching, the rule can go."""
    root = ET.fromstring("<a><R><Z>5</Z></R><R><Z>1</Z></R></a>")
    assert root.findall("R[Z>3.0]") == []


def test_no_report_selects_with_a_comparison():
    offenders = []
    for path in list(ROOT.glob("wrappers/*/script/*.py")) + list(ROOT.glob("pipelines/*/script/*.py")):
        for n, line in enumerate(path.read_text(errors="replace").splitlines(), 1):
            if COMPARISON.search(line):
                offenders.append(f"{path.relative_to(ROOT)}:{n}: {line.strip()}")
    assert not offenders, "\n".join(offenders)


def test_edstats_outliers_are_selected_in_python():
    from ccp4i2.wrappers.edstats.script.edstats_report import _outliers
    group = ET.fromstring(
        "<ligands><Residue><Name>NUT</Name><ZDpa>30.2</ZDpa><ZDma>-0.4</ZDma></Residue>"
        "<Residue><Name>SO4</Name><ZDpa>1.0</ZDpa><ZDma>-3.5</ZDma></Residue></ligands>")
    assert [r.findtext("Name") for r in _outliers(group, "ZDpa", above=3.0)] == ["NUT"]
    assert [r.findtext("Name") for r in _outliers(group, "ZDma", below=-3.0)] == ["SO4"]
