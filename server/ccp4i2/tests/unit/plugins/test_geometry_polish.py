"""geometry_polish: the pieces that decide what gets changed.

Pure-Python where it can be (MolProbity parsers, the density guard, the CSV
writer, the def.xml and the plugin import); gemmi for the NCS transplant
maths on a toy structure; nothing here needs Coot or MolProbity binaries.
"""
import csv
import importlib
import math
from pathlib import Path

import pytest

from ccp4i2.wrappers.geometry_polish.script import polish_lib as P

gemmi = pytest.importorskip("gemmi")

ROTALYZE = """residue:occupancy:score%:chi1:chi2:chi3:chi4:evaluation:rotamer
 A  12  SER:1.00:0.5:62.1::::OUTLIER:OUTLIER
 A  13  LEU:1.00:45.2:-65.0:175.1:::Favored:mt
 A  14  GLY:1.00:::::::
 B 379  GLU:0.50:1.8:-70.0:160.0:20.0::Allowed:mt-10
SUMMARY: 1 outliers out of 3 residues (numbering may differ)
"""

RAMALYZE = """residue:score%:phi:psi:evaluation:type
 A  12  SER:0.02:-150.0:-170.0:OUTLIER:General
 A  13  LEU:55.60:-65.0:-40.0:Favored:General
 A  14  GLY:3.10:80.0:10.0:Allowed:Glycine
SUMMARY: 1 Favored, 1 Allowed, 1 OUTLIER
"""

CLASHSCORE = """Using electron cloud x-H distances and vdW radii
Bad Clashes >= 0.4 Angstrom:
 C 117  PHE  HB3  C 138  SER  HB3 :1.153
 A  12  SER  OG   B  45  LYS  NZ  :0.612
clashscore = 4.26
"""


def test_parse_rotalyze_keeps_outliers_and_skips_glycine():
    rot = P.parse_rotalyze(ROTALYZE)
    assert rot[("A", 12)] == dict(score=0.5, evaluation="OUTLIER", rotamer="OUTLIER", resname="SER")
    assert rot[("A", 13)]["rotamer"] == "mt"
    assert rot[("B", 379)]["evaluation"] == "Allowed"
    assert ("A", 14) not in rot  # no score: glycine has no rotamer


def test_parse_ramalyze():
    ram = P.parse_ramalyze(RAMALYZE)
    assert ram[("A", 12)]["evaluation"] == "OUTLIER"
    assert ram[("A", 14)] == dict(score=3.1, evaluation="Allowed", type="Glycine", resname="GLY")
    assert len(ram) == 3


def test_parse_clashscore_and_pairs():
    score, clashes = P.parse_clashscore(CLASHSCORE)
    assert score == 4.26
    assert len(clashes) == 2
    assert clashes[0]["overlap"] == 1.153 and clashes[0]["atom1"] == "HB3"
    assert P.clash_pair(clashes[1]) == (("A", 12), ("B", 45))


def test_parse_clash_line_rejects_prose():
    assert P.parse_clash_line("Bad Clashes >= 0.4 Angstrom:") is None


@pytest.mark.parametrize("before, after, ok", [
    (4.0, 3.9, True),      # unchanged within slack
    (4.0, 3.0, True),      # 75% of a strong peak is still plenty
    (4.0, 2.0, False),     # lost half of strong density
    (1.0, 0.75, True),     # within the absolute slack
    (1.0, 0.6, False),     # below the minimum, and more than slack lost
    (0.3, 0.1, True),      # nothing to guard: no density to begin with
    (2.0, float("nan"), True),
])
def test_density_guard(before, after, ok):
    assert P.density_ok(before, after, guard_fraction=0.6, min_density=0.8) is ok


def _toy(seq=("ALA", "SER", "ALA")):
    """Two identical chains of a short helix-ish peptide, chain B rotated and shifted."""
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(60, 60, 60, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    template = [
        ("N", "N", (0.0, 0.0, 0.0)), ("CA", "C", (1.46, 0.0, 0.0)), ("C", "C", (2.0, 1.4, 0.0)),
        ("O", "O", (1.3, 2.4, 0.0)), ("CB", "C", (2.0, -0.8, 1.2)), ("OG", "O", (3.4, -0.9, 1.3)),
    ]
    for chain_name, rot, shift in (("A", 0.0, (0, 0, 0)), ("B", 0.7, (20, 5, 3))):
        chain = gemmi.Chain(chain_name)
        for i, name in enumerate(seq):
            res = gemmi.Residue()
            res.name = name
            res.seqid = gemmi.SeqId(i + 1, " ")
            for atom_name, element, (x, y, z) in template:
                if name == "ALA" and atom_name == "OG":
                    continue
                atom = gemmi.Atom()
                atom.name = atom_name
                atom.element = gemmi.Element(element)
                x = x + 3.8 * i
                rx = x * math.cos(rot) - y * math.sin(rot)
                ry = x * math.sin(rot) + y * math.cos(rot)
                atom.pos = gemmi.Position(rx + shift[0], ry + shift[1], z + shift[2])
                atom.occ = 1.0
                atom.b_iso = 20.0
                res.add_atom(atom)
            chain.add_residue(res)
        model.add_chain(chain)
    st.add_model(model)
    st.setup_entities()
    return st


def test_transplant_sidechain_moves_only_the_side_chain_and_adds_missing_atoms():
    st = _toy()
    target = P.find_residue(st, "A", 2)
    target.remove_atom("OG", " ")             # a truncated side chain
    assert target.find_atom("OG", "*") is None
    ca_before = target.find_atom("CA", "*").pos
    rmsd = P.transplant_sidechain(st, "A", 2, "B")
    assert rmsd is not None and rmsd < 0.05  # identical copies superpose exactly
    og = target.find_atom("OG", "*")
    assert og is not None
    assert target.find_atom("CA", "*").pos.dist(ca_before) < 1e-6
    # the transplanted OG lands where chain A's own OG was, within rounding
    expected = gemmi.Position(3.4 + 3.8, -0.9, 1.3)
    assert og.pos.dist(expected) < 0.05
    assert og.b_iso == target.find_atom("CB", "*").b_iso


def test_transplant_sidechain_refuses_a_different_residue_type():
    st = _toy()
    P.find_residue(st, "B", 2).name = "THR"
    assert P.transplant_sidechain(st, "A", 2, "B") is None


def test_transplant_backbone_window_needs_flanks():
    st = _toy()
    # the whole chain: nothing left to superpose on
    assert P.transplant_backbone_window(st, "A", 1, 3, "B") is None


def test_transplant_backbone_window_moves_the_window_only():
    st = _toy()
    before = {a.name: a.pos for a in P.find_residue(st, "A", 2)}
    flank = P.find_residue(st, "A", 1).find_atom("CA", "*").pos
    rmsd = P.transplant_backbone_window(st, "A", 2, 2, "B")
    assert rmsd is not None and rmsd < 0.05
    after = P.find_residue(st, "A", 2)
    for name, pos in before.items():          # identical copies: nothing moves measurably
        assert after.find_atom(name, "*").pos.dist(pos) < 0.05
    assert P.find_residue(st, "A", 1).find_atom("CA", "*").pos.dist(flank) < 1e-6


def test_write_window_keeps_only_the_range(tmp_path):
    st = _toy()
    path = P.write_window(st, "A", 2, 3, str(tmp_path / "w.pdb"))
    mini = gemmi.read_structure(path)
    assert [r.seqid.num for r in mini[0]["A"]] == [2, 3]
    assert len(mini[0]) == 1


def test_write_neighbourhood_includes_partners_and_nearby(tmp_path):
    st = _toy()
    path = P.write_neighbourhood(st, [("A", 1)], str(tmp_path / "n.pdb"), radius=5.0)
    mini = gemmi.read_structure(path)
    names = {(ch.name, r.seqid.num) for ch in mini[0] for r in ch}
    assert ("A", 1) in names and ("A", 2) in names
    assert not any(ch == "B" for ch, _ in names)  # chain B is 20 A away


def test_changes_csv_has_fixed_columns(tmp_path):
    rows = [dict(chain="A", resno=12, resname="SER", max_shift=1.2, action="rotamers:autofit",
                 density_before=1.1, density_after=1.0, rotamer_before="OUTLIER", rotamer_after="p",
                 rama_before="Favored", rama_after="Favored")]
    path = P.write_changes(str(tmp_path / "c.csv"), rows)
    with open(path) as handle:
        read = list(csv.DictReader(handle))
    assert list(read[0].keys()) == list(P.CHANGES_COLUMNS)
    assert read[0]["action"] == "rotamers:autofit"
    assert read[0]["diff_density_before"] == ""


def test_options_reject_unknown_knob():
    with pytest.raises(TypeError):
        P.Options(bogus=1)


def test_plugin_imports_without_coot_or_molprobity():
    module = importlib.import_module("ccp4i2.wrappers.geometry_polish.script.geometry_polish")
    assert module.geometry_polish.TASKNAME == "geometry_polish"
    assert "molprobity.rotalyze" in module.geometry_polish.AUXILIARY_PROGRAMS
    report = importlib.import_module("ccp4i2.wrappers.geometry_polish.script.geometry_polish_report")
    assert report.geometry_polish_report.TASKNAME == "geometry_polish"


def test_chain_selection_parsing():
    from ccp4i2.wrappers.geometry_polish.script.geometry_polish import _chains
    assert _chains("") is None
    assert _chains("A, B;C") == ["A", "B", "C"]


def test_def_xml_declares_the_contract():
    from ccp4i2.core.task_manager.def_xml_handler import DefXmlParser
    here = Path(__file__).resolve().parents[3] / "wrappers" / "geometry_polish" / "script"
    parsed = DefXmlParser().parse_def_xml(here / "geometry_polish.def.xml")
    text = (here / "geometry_polish.def.xml").read_text()
    for name in ("XYZIN", "FPHIOUT", "DIFFPHIOUT", "DICT", "XYZOUT", "CHANGES_CSV",
                 "DO_ROTAMERS", "DO_RAMA", "DO_CLASHES", "DO_PEPTIDE_FLIPS", "USE_NCS",
                 "DENSITY_GUARD", "MIN_DENSITY", "REFINE_CYCLES", "MAX_CLASH_ROTAMERS",
                 "MIN_CLASH_OVERLAP", "CHAIN_SELECTION"):
        assert 'id="{}"'.format(name) in text
    assert parsed is not None


def test_task_is_registered():
    from ccp4i2.core.tasks import TASKS
    from ccp4i2.core.citations import TASK_CITES
    task = TASKS["geometry_polish"]
    assert task.pluginPath.endswith(":geometry_polish")
    assert TASK_CITES["geometry_polish"] == ["coot", "molprobity"]
