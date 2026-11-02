"""Servalcat, as its author found it (alpha tracker, 2026-10-16).

- MetalCoord never ran from the interface: nothing found the metal sites, so
  the list of monomers to restrain stayed empty and the job succeeded with no
  metal restraints. formats/metal_sites finds them, single ions included.
- From Servalcat 0.4.150 the R against intensities is a conventional R,
  reported as Rwork/Rfree/R, where earlier versions reported R1work/R1free/R1.
  The R name then no longer says which data were refined against; the CC
  names (CCI, CCF) do. Both namings must give the right report columns and the
  right KPIs, and an R1 must never be filed as an R.
- Report table cells written with HTML entities or an unclosed <br> were
  shown literally, and a grouped heading's subtitles were dropped.

gemmi and the stdlib only; no CCP4 binaries.
"""
import xml.etree.ElementTree as ET
from types import SimpleNamespace

import gemmi

from ccp4i2.lib.utils.formats.metal_sites import (
    metal_monomer_codes, metal_monomer_codes_in_file,
)
from ccp4i2.report.elements import _set_cell_content
from ccp4i2.wrappers.servalcat.script import servalcat as servalcat_wrapper
from ccp4i2.wrappers.servalcat.script.servalcat_report import (
    _r_and_cc, _r_and_cc_columns, servalcat_report,
)


# --- metal sites ------------------------------------------------------------

def _structure(*residues):
    st = gemmi.Structure()
    model = gemmi.Model("1")
    chain = gemmi.Chain("A")
    for seq, (name, elements) in enumerate(residues, start=1):
        residue = gemmi.Residue()
        residue.name = name
        residue.seqid = gemmi.SeqId(seq, " ")
        for k, element in enumerate(elements):
            atom = gemmi.Atom()
            atom.name = f"{element}{k}"
            atom.element = gemmi.Element(element)
            atom.pos = gemmi.Position(seq, k, 0)
            residue.add_atom(atom)
        chain.add_residue(residue)
    model.add_chain(chain)
    st.add_model(model)
    return st


def test_single_ions_and_metal_containing_monomers_are_found():
    st = _structure(("ALA", ["N", "C", "C", "O", "C"]), ("ZN", ["Zn"]),
                    ("HEM", ["Fe", "N", "N", "N", "N", "C"]), ("ZN", ["Zn"]),
                    ("HOH", ["O"]), ("ATP", ["P", "P", "P", "O", "O"]))
    assert metal_monomer_codes(st) == ["HEM", "ZN"]


def test_a_model_without_metals_has_none():
    assert metal_monomer_codes(_structure(("ALA", ["N", "C", "C", "O"]), ("HOH", ["O"]))) == []


def test_read_from_a_file_and_an_unreadable_file_has_none(tmp_path):
    path = tmp_path / "model.pdb"
    _structure(("GLY", ["N", "C", "C", "O"]), ("CA", ["Ca"])).write_pdb(str(path))
    assert metal_monomer_codes_in_file(str(path)) == ["CA"]
    assert metal_monomer_codes_in_file(str(tmp_path / "missing.pdb")) == []


# --- R and CC under both namings -------------------------------------------

def _program_xml(**final):
    """Two cycles of servalcat summary, with the given final-cycle keys."""
    cycles = "".join(
        "<cycle><Ncyc>%d</Ncyc><data><summary>%s</summary></data></cycle>" % (
            n, "".join(f"<{k}>{v}</{k}>" for k, v in final.items()))
        for n in (0, 1))
    return ET.fromstring(f"<SERVALCAT>{cycles}</SERVALCAT>")


OLD_INTENSITIES = dict(R1work=0.20, R1free=0.25, CCIworkavg=0.90, CCIfreeavg=0.85)
NEW_INTENSITIES = dict(Rwork=0.21, Rfree=0.26, CCIworkavg=0.90, CCIfreeavg=0.85)
AMPLITUDES = dict(Rwork=0.22, Rfree=0.27, CCFworkavg=0.88, CCFfreeavg=0.83)
NO_FREE_SET_NEW = dict(R=0.23, CCIavg=0.89)


def test_which_r_and_cc_were_reported():
    assert _r_and_cc(_program_xml(**OLD_INTENSITIES)) == ("R1work", "R1free", "I")
    assert _r_and_cc(_program_xml(**NEW_INTENSITIES)) == ("Rwork", "Rfree", "I")
    assert _r_and_cc(_program_xml(**AMPLITUDES)) == ("Rwork", "Rfree", "F")
    assert _r_and_cc(_program_xml(**NO_FREE_SET_NEW)) == ("R", None, "I")
    assert _r_and_cc(_program_xml(FSCaverage=0.8)) is None


def _columns(program_xml):
    report = servalcat_report.__new__(servalcat_report)
    return [title for title, _ in _r_and_cc_columns(report.getCycleData(xmlnode=program_xml))]


def test_new_intensity_output_keeps_its_correlations():
    """Under 0.4.150 the old code took Rwork to mean amplitudes and looked for
    CCF, so the correlation columns were all '-'."""
    assert _columns(_program_xml(**NEW_INTENSITIES)) == ["Rwork", "Rfree", "⟨CCIwork⟩", "⟨CCIfree⟩"]


def test_old_intensity_and_amplitude_output_are_unchanged():
    assert _columns(_program_xml(**OLD_INTENSITIES)) == ["R1work", "R1free", "⟨CCIwork⟩", "⟨CCIfree⟩"]
    assert _columns(_program_xml(**AMPLITUDES)) == ["Rwork", "Rfree", "⟨CCFwork⟩", "⟨CCFfree⟩"]
    assert _columns(_program_xml(**NO_FREE_SET_NEW)) == ["R", "⟨CCI⟩"]


class _Value:
    def __init__(self):
        self.value = None

    def set(self, value):
        self.value = value


def _kpis(summary, free_set=True):
    names = ("RFactor", "RFree", "R", "R1Factor", "R1Free", "R1",
             "CCFwork_avg", "CCFfree_avg", "CCF_avg", "CCIwork_avg", "CCIfree_avg", "CCI_avg")
    indicators = SimpleNamespace(**{name: _Value() for name in names})
    plugin = SimpleNamespace(container=SimpleNamespace(
        outputData=SimpleNamespace(PERFORMANCEINDICATOR=indicators),
        inputData=SimpleNamespace(FREERFLAG=SimpleNamespace(isSet=lambda: free_set)),
        controlParameters=SimpleNamespace(DATA_METHOD="xtal"),
    ))
    servalcat_wrapper.servalcat._extractPerformanceIndicators(plugin, [{"data": {"summary": summary}}])
    return {name: getattr(indicators, name).value for name in names
            if getattr(indicators, name).value is not None}


def test_kpis_under_each_naming():
    assert _kpis(OLD_INTENSITIES) == {"R1Factor": 0.20, "R1Free": 0.25,
                                      "CCIwork_avg": 0.90, "CCIfree_avg": 0.85}
    # A conventional R: the same KPIs as amplitude refinement, never R1
    assert _kpis(NEW_INTENSITIES) == {"RFactor": 0.21, "RFree": 0.26,
                                      "CCIwork_avg": 0.90, "CCIfree_avg": 0.85}
    assert _kpis(AMPLITUDES) == {"RFactor": 0.22, "RFree": 0.27,
                                 "CCFwork_avg": 0.88, "CCFfree_avg": 0.83}
    assert _kpis(NO_FREE_SET_NEW, free_set=False) == {"R": 0.23, "CCI_avg": 0.89}


def test_a_missing_free_r_is_not_recorded_as_zero():
    assert "RFree" not in _kpis(dict(Rwork=0.2, CCFworkavg=0.9))


# --- table cells ------------------------------------------------------------

def _cell(content):
    th = ET.Element("th")
    _set_cell_content(th, content)
    return ET.tostring(th, encoding="unicode")


def test_html_entities_and_br_become_markup_not_text():
    assert _cell("Critical<br>distance (&Aring;)") == "<th>Critical<br />distance (Å)</th>"
    assert _cell("B-value (&Aring;<sup>2</sup>)") == "<th>B-value (Å<sup>2</sup>)</th>"
    assert _cell("&nbsp;bond") == "<th> bond</th>"


def test_xml_entities_and_plain_text_are_unchanged():
    assert _cell("a &amp; b &lt; c") == "<th>a &amp; b &lt; c</th>"
    assert _cell("Rwork") == "<th>Rwork</th>"


def test_a_grouped_heading_keeps_its_subtitles():
    """One heading over Bond and Angle, and the subtitles under it in a row of
    their own, so a reader can put every heading over its own column."""
    from ccp4i2.report.elements import Table

    table = Table()
    table.addData(title="Cycle", data=[0, 1])
    table.addData(title="RMSD", subtitle="Bond", data=[0.01, 0.02])
    table.addData(subtitle="Angle", data=[1.1, 1.2])
    table.addData(title="Weight", data=["-", 0.1])
    groups, subtitles = table.data_as_etree().find(".//thead")
    assert [(th.text, th.get("colspan")) for th in groups] == [
        ("Cycle", None), ("RMSD", "2"), ("Weight", None)]
    assert [th.text for th in subtitles] == [None, "Bond", "Angle", None]
