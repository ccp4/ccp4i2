"""
checkMonomeCoverage must distinguish gemmi's fatal warnings from its advisories.

gemmi's prepare_topology reports every unmatched atom as "definition not
found", then appends a parenthesised remark saying what refinement will do
about it.  "linkage should remove this atom" means the atom is dropped when
the link is applied -- refmac and servalcat both do this silently and
correctly.  Typically it is the extra N-terminal hydrogens (H2/H3) left on a
residue that is no longer a chain terminus.

Treating that advisory as SEVERITY_ERROR blocked submission of a model that
refinement would have handled by itself, which is what happened to a real
Moorhen-edited structure.  A genuinely uncovered ligand must still fail.
"""

import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi (CCP4 python)")

import os

from ccp4i2.core.CCP4PluginScript import CPluginScript

pytestmark = pytest.mark.skipif(
    not (os.environ.get("CLIBD_MON")
         or (os.environ.get("CCP4")
             and os.path.isdir(os.path.join(os.environ["CCP4"], "lib", "data", "monomers")))),
    reason="needs the CCP4 monomer library (CLIBD_MON)",
)


class _Stub(CPluginScript):
    """Bare CPluginScript: checkMonomeCoverage needs no plugin machinery."""

    TASKNAME = "stub"

    def __init__(self):
        self.reports = []

    def appendErrorReport(self, code, details, name=None, severity=None):
        self.reports.append((code, details))


def _dipeptide(with_extra_nterm_h):
    """Two glycines; optionally put H2/H3 on the second (non-terminal) one."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(40, 40, 40, 90, 90, 90)
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")

    geom = [
        # (name, element, x, y, z) -- roughly correct glycine backbone
        [("N", "N", 0.0, 0.0, 0.0), ("CA", "C", 1.458, 0.0, 0.0),
         ("C", "C", 2.0, 1.42, 0.0), ("O", "O", 1.25, 2.39, 0.0)],
        [("N", "N", 3.33, 1.55, 0.0), ("CA", "C", 4.0, 2.84, 0.0),
         ("C", "C", 5.5, 2.7, 0.0), ("O", "O", 6.2, 3.7, 0.0)],
    ]
    for idx, atoms in enumerate(geom, start=1):
        res = gemmi.Residue()
        res.name = "GLY"
        res.seqid = gemmi.SeqId(str(idx))
        for nm, el, x, y, z in atoms:
            a = gemmi.Atom()
            a.name = nm
            a.element = gemmi.Element(el)
            a.pos = gemmi.Position(x, y, z)
            a.occ = 1.0
            a.b_iso = 20.0
            res.add_atom(a)
        if idx == 2 and with_extra_nterm_h:
            for nm in ("H2", "H3"):
                a = gemmi.Atom()
                a.name = nm
                a.element = gemmi.Element("H")
                a.pos = gemmi.Position(3.3, 0.6, 0.3)
                a.occ = 1.0
                a.b_iso = 20.0
                res.add_atom(a)
        chain.add_residue(res)

    model.add_chain(chain)
    st.add_model(model)
    return st


def _write(st, path):
    st.make_mmcif_document().write_file(str(path))
    return str(path)


def test_linkage_removable_atoms_do_not_block(tmp_path):
    """H2/H3 on a mid-chain residue: refinement drops them, so we must pass."""
    path = _write(_dipeptide(with_extra_nterm_h=True), tmp_path / "extra_h.cif")

    plugin = _Stub()
    result = plugin.checkMonomeCoverage(path)

    assert result == CPluginScript.SUCCEEDED, \
        f"linkage advisory wrongly blocked: {plugin.reports}"
    assert plugin.reports == []


def test_clean_model_passes(tmp_path):
    path = _write(_dipeptide(with_extra_nterm_h=False), tmp_path / "clean.cif")

    plugin = _Stub()
    assert plugin.checkMonomeCoverage(path) == CPluginScript.SUCCEEDED
    assert plugin.reports == []


def test_unknown_ligand_still_fails(tmp_path):
    """The check must keep catching what it exists to catch."""
    st = _dipeptide(with_extra_nterm_h=False)
    res = gemmi.Residue()
    res.name = "ZZZ"          # no such monomer in the CCP4 library
    res.seqid = gemmi.SeqId("900")
    res.het_flag = "H"
    for nm, el in [("C1", "C"), ("O1", "O")]:
        a = gemmi.Atom()
        a.name = nm
        a.element = gemmi.Element(el)
        a.pos = gemmi.Position(20, 20, 20)
        a.occ = 1.0
        a.b_iso = 30.0
        res.add_atom(a)
    st[0]["A"].add_residue(res)
    path = _write(st, tmp_path / "unknown.cif")

    plugin = _Stub()
    assert plugin.checkMonomeCoverage(path) == CPluginScript.FAILED
    assert plugin.reports, "an uncovered ligand must be reported"
    assert "ZZZ" in plugin.reports[0][1]


def test_blank_atom_name_is_reported_legibly(tmp_path):
    """The old message read '(replace  with N)', which told the user nothing."""
    st = _dipeptide(with_extra_nterm_h=False)
    for atom in st[0]["A"][1]:
        if atom.name == "O":
            atom.name = ""
            break
    path = _write(st, tmp_path / "blank.cif")

    plugin = _Stub()
    assert plugin.checkMonomeCoverage(path) == CPluginScript.FAILED
    detail = plugin.reports[0][1]
    assert "No atom name in the coordinate file" in detail
    assert "replace  with" not in detail


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
