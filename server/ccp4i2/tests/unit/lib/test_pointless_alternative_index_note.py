"""The Pointless report's note that the data could be indexed another way
and "no reference data has been given" must appear only when no reference
was given. It read `not HKLREF or XYZIN`, so a coordinate reference raised
the note that advises giving one."""
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import pytest

pointless_report = pytest.importorskip(
    "ccp4i2.wrappers.pointless.script.pointless_report").pointless_report

ALTERNATIVES = """
  <NumberPossibleReindexing>2</NumberPossibleReindexing>
  <PossibleReindexing><ReindexOperator>[h,k,l]</ReindexOperator>
    <CellDiff>0.00</CellDiff><DifferentCell>false</DifferentCell></PossibleReindexing>
  <PossibleReindexing><ReindexOperator>[-h,-k,l]</ReindexOperator>
    <CellDiff>0.00</CellDiff><DifferentCell>false</DifferentCell></PossibleReindexing>
"""


def note(reference):
    xml = ET.fromstring("<POINTLESS>" + ALTERNATIVES + reference + "</POINTLESS>")
    parent = []
    pointless_report.AlternativeIndexWarning(SimpleNamespace(xmlnode=xml), parent=parent)
    return "".join(parent)


def test_note_without_a_reference():
    assert "no reference data has been given" in note("")


@pytest.mark.parametrize("reference", [
    "<ReflectionFile stream='XYZIN'><FileName>model.cif</FileName></ReflectionFile>",
    "<ReflectionFile stream='HKLREF'><FileName>ref.mtz</FileName></ReflectionFile>",
    "<BestReindex><XYZIN>model.cif</XYZIN><ReindexOperator>[-h,-k,l]</ReindexOperator></BestReindex>",
])
def test_no_note_with_a_reference(reference):
    assert note(reference) == ""
