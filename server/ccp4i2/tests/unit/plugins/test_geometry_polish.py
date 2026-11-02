"""geometry_polish: the pieces that decide what gets changed.

Pure-Python where it can be (MolProbity parsers, the density guard, the judge,
the CSV writer, the def.xml and the plugin import); gemmi for the NCS
transplant maths and windows on a toy structure; nothing here needs Coot or
MolProbity binaries.
"""
import csv
import importlib
import math
from pathlib import Path

import pytest

from ccp4i2.wrappers.geometry_polish.script import polish_lib as P

gemmi = pytest.importorskip("gemmi")

# Real molprobity output, 2026-07 CCP4: 4-digit numbers, an insertion code,
# and a 2-character chain (from an mmCIF input).
ROTALYZE = """residue:occupancy:score%:chi1:chi2:chi3:chi4:evaluation:rotamer
 A1001  SER:1.00:76.6:60.2::::Favored:p
 A1005  LYS:1.00:58.9:182.6:170.5:189.4:163.2:Favored:tttt
 A1005A LYS:1.00:0.3:285.1:183.0:200.5:181.9:OUTLIER:OUTLIER
 A1004  GLY:1.00:::::::
AB  11  SER:1.00:78.4:60.7::::Favored:p
 B 379  GLU:0.50:1.8:-70.0:160.0:20.0::Allowed:mt-10
SUMMARY: 1 outliers out of 5 residues (numbering may differ)
"""

RAMALYZE = """residue:score%:phi:psi:evaluation:type
 A1002  ASP:4.06:-90.88:-54.23:Favored:General
 A1005  LYS:71.59:-56.27:-38.75:Favored:General
 A1005A LYS:0.02:-150.0:-170.0:OUTLIER:General
AB  14  GLY:6.69:-44.64:-54.28:Allowed:Glycine
SUMMARY: 2 Favored, 1 Allowed, 1 OUTLIER
"""

CLASHSCORE = """Using electron cloud x-H distances and vdW radii
Adding H/D atoms with reduce...
Bad Clashes >= 0.4 Angstrom:
 C1106  PHE  HB3  C1127A SER  HB3 :1.068
 C1106  PHE  CB   C1127A SER  HB3 :1.052
CX 116  PHE  HB3 CX 137A SER  CB  :0.812
 A  12  SER  OG   B  45  LYS  NZ  :0.612
clashscore = 4.26
"""


def test_residue_field_is_positional():
    assert P.parse_residue_field(" A1005A LYS") == (("A", 1005, "A"), "LYS")
    assert P.parse_residue_field("AB  11  SER") == (("AB", 11, ""), "SER")
    assert P.parse_residue_field(" A  12  SER") == (("A", 12, ""), "SER")
    assert P.parse_residue_field("garbage") is None


def test_parse_rotalyze_four_digits_icodes_and_two_char_chains():
    rot = P.parse_rotalyze(ROTALYZE)
    assert rot[("A", 1001, "")]["rotamer"] == "p"
    assert rot[("A", 1005, "")]["evaluation"] == "Favored"
    assert rot[("A", 1005, "A")] == dict(score=0.3, evaluation="OUTLIER", rotamer="OUTLIER", resname="LYS")
    assert rot[("AB", 11, "")]["resname"] == "SER"
    assert rot[("B", 379, "")]["evaluation"] == "Allowed"
    assert ("A", 1004, "") not in rot  # no score: glycine has no rotamer
    assert len(rot) == 5


def test_parse_ramalyze_four_digits_icodes_and_two_char_chains():
    ram = P.parse_ramalyze(RAMALYZE)
    assert ram[("A", 1005, "")]["evaluation"] == "Favored"
    assert ram[("A", 1005, "A")]["evaluation"] == "OUTLIER"
    assert ram[("AB", 14, "")] == dict(score=6.69, evaluation="Allowed", type="Glycine", resname="GLY")
    assert len(ram) == 4


def test_parse_clashscore_four_digits_icodes_and_two_char_chains():
    score, clashes = P.parse_clashscore(CLASHSCORE)
    assert score == 4.26
    assert len(clashes) == 4
    first = clashes[0]
    assert first["key1"] == ("C", 1106, "") and first["atom1"] == "HB3"
    assert first["key2"] == ("C", 1127, "A") and first["resname2"] == "SER" and first["overlap"] == 1.068
    assert clashes[1]["atom1"] == "CB"
    assert clashes[2]["key1"] == ("CX", 116, "") and clashes[2]["key2"] == ("CX", 137, "A")
    assert clashes[2]["atom2"] == "CB"
    assert P.clash_pair(clashes[3]) == (("A", 12, ""), ("B", 45, ""))


def test_parse_clash_line_rejects_prose():
    assert P.parse_clash_line("Bad Clashes >= 0.4 Angstrom:") is None
    assert P.parse_clash_line("Adding H/D atoms with reduce...") is None


def test_new_clashes_diff_is_by_atom_pair_and_threshold():
    _, before = P.parse_clashscore(CLASHSCORE)
    _, after = P.parse_clashscore(CLASHSCORE + " D 200  LEU  CD1  D 230  VAL  CG1 :0.950\n"
                                  " D 201  LEU  CD1  D 231  VAL  CG1 :0.500\n")
    added = P.new_clashes(before, after, 0.8)
    assert [(c["key1"], c["overlap"]) for c in added] == [(("D", 200, ""), 0.95)]
    # the same atom pair written the other way round is not new
    swapped = dict(after[0])
    swapped["key1"], swapped["key2"] = after[0]["key2"], after[0]["key1"]
    swapped["atom1"], swapped["atom2"] = after[0]["atom2"], after[0]["atom1"]
    assert P.new_clashes(before, [swapped], 0.8) == []


@pytest.mark.parametrize("before, after, ok", [
    (4.0, 3.9, True),      # unchanged within slack
    (4.0, 3.0, True),      # 75% of a strong peak is still plenty
    (4.0, 2.0, False),     # lost half of strong density
    (1.0, 0.85, True),     # within the absolute slack and above the floor
    (1.0, 0.6, False),     # below the minimum
    (0.5, 0.2, False),     # a weak side chain may not lose more than half the slack
    (0.3, 0.1, False),     # nor may a weaker one
    (0.66, 0.66, True),    # a weak flank that did not move is not a veto
    (0.26, 0.23, True),    # nor one that jiggled
    (2.0, float("nan"), True),
])
def test_density_guard(before, after, ok):
    assert P.density_ok(before, after, guard_fraction=0.6, min_density=0.8) is ok


def _toy(seq=("ALA", "SER", "ALA", "LEU", "ALA"), icode_at=None):
    """Two identical chains of a short peptide, chain B rotated and shifted.

    ``icode_at`` gives one residue the previous residue's number with
    insertion code A (a Kabat-style pair), in both chains.
    """
    st = gemmi.Structure()
    st.cell = gemmi.UnitCell(60, 60, 60, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")
    template = [
        ("N", "N", (0.0, 0.0, 0.0)), ("CA", "C", (1.46, 0.0, 0.0)), ("C", "C", (2.0, 1.4, 0.0)),
        ("O", "O", (1.3, 2.4, 0.0)), ("CB", "C", (2.0, -0.8, 1.2)), ("OG", "O", (3.4, -0.9, 1.3)),
        ("CG", "C", (3.4, -0.9, 1.3)), ("CD1", "C", (4.2, -0.2, 2.3)), ("CD2", "C", (4.0, -2.2, 1.0)),
    ]
    wanted = {"ALA": {"N", "CA", "C", "O", "CB"}, "SER": {"N", "CA", "C", "O", "CB", "OG"},
              "LEU": {"N", "CA", "C", "O", "CB", "CG", "CD1", "CD2"}}
    for chain_name, rot, shift in (("A", 0.0, (0, 0, 0)), ("B", 0.7, (20, 5, 3))):
        chain = gemmi.Chain(chain_name)
        for i, name in enumerate(seq):
            res = gemmi.Residue()
            res.name = name
            if icode_at is not None and i == icode_at:
                res.seqid = gemmi.SeqId(i, "A")
            else:
                res.seqid = gemmi.SeqId(i + 1, " ")
            for atom_name, element, (x, y, z) in template:
                if atom_name not in wanted[name]:
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


def test_residue_keys_carry_insertion_codes():
    st = _toy(icode_at=2)   # residues 1, 2, 2A, 4, 5
    assert P.find_residue(st, "A", 2).name == "SER"
    assert P.find_residue(st, "A", 2, "A").name == "ALA"
    assert P.locate(st, ("A", 2, "A")) == 2
    assert P.keys_at(st, "A", [1, 2]) == [("A", 2, ""), ("A", 2, "A")]
    assert P.coot_cid(("A", 2, "A")) == "//A/2.A"
    assert P.coot_cid(("A", 2, "")) == "//A/2"


def test_transplant_sidechain_targets_the_right_icode_twin():
    st = _toy(seq=("ALA", "LEU", "LEU", "ALA", "ALA"), icode_at=2)   # 2 LEU, 2A LEU
    twin = P.find_residue(st, "A", 2, "A")
    twin.remove_atom("CD1", " ")
    plain = P.find_residue(st, "A", 2)
    plain_cd1 = plain.find_atom("CD1", "*").pos
    assert P.transplant_sidechain(st, "A", 2, "B", icode="A") is not None
    assert twin.find_atom("CD1", "*") is not None          # the 2A copy was rebuilt
    assert plain.find_atom("CD1", "*").pos.dist(plain_cd1) < 1e-6   # residue 2 untouched


def test_transplant_sidechain_moves_only_the_side_chain_and_keeps_occupancy():
    st = _toy()
    target = P.find_residue(st, "A", 2)
    target.remove_atom("OG", " ")             # a truncated side chain
    target.find_atom("CB", "*").occ = 0.7
    donor = P.find_residue(st, "B", 2)
    donor.find_atom("OG", "*").occ = 0.3
    ca_before = target.find_atom("CA", "*").pos
    rmsd = P.transplant_sidechain(st, "A", 2, "B")
    assert rmsd is not None and rmsd < 0.05  # identical copies superpose exactly
    og = target.find_atom("OG", "*")
    assert og is not None
    assert abs(og.occ - 0.7) < 1e-6           # the target's occupancy, not the donor's
    assert target.find_atom("CA", "*").pos.dist(ca_before) < 1e-6
    expected = gemmi.Position(3.4 + 3.8, -0.9, 1.3)
    assert og.pos.dist(expected) < 0.05
    assert og.b_iso == target.find_atom("CB", "*").b_iso


def test_transplant_refuses_altlocs_and_different_residue_types():
    st = _toy()
    P.find_residue(st, "B", 2).name = "THR"
    assert P.transplant_sidechain(st, "A", 2, "B") is None
    st = _toy()
    P.find_residue(st, "B", 2).find_atom("OG", "*").altloc = "A"
    assert P.transplant_sidechain(st, "A", 2, "B") is None
    assert P.transplant_backbone_window(st, "A", 2, 2, "B") is None


def test_transplant_backbone_window_needs_flanks():
    st = _toy()
    assert P.transplant_backbone_window(st, "A", 1, 5, "B") is None   # the whole chain


def test_transplant_backbone_window_moves_the_window_only():
    st = _toy()
    before = {a.name: a.pos for a in P.find_residue(st, "A", 3)}
    flank = P.find_residue(st, "A", 1).find_atom("CA", "*").pos
    rmsd = P.transplant_backbone_window(st, "A", 3, 3, "B")
    assert rmsd is not None and rmsd < 0.05
    after = P.find_residue(st, "A", 3)
    for name, pos in before.items():          # identical copies: nothing moves measurably
        assert after.find_atom(name, "*").pos.dist(pos) < 0.05
    assert P.find_residue(st, "A", 1).find_atom("CA", "*").pos.dist(flank) < 1e-6


def test_copy_window_brings_every_atom_of_every_residue():
    src, dst = _toy(), _toy()
    for res in P.chain_residues(src, "A"):
        for atom in res:
            atom.pos = gemmi.Position(atom.pos.x, atom.pos.y, atom.pos.z + 1.0)
    P.copy_window(src, dst, "A", [1, 2, 3])
    for index in (1, 2, 3):
        for atom in P.residue_at(dst, "A", index):
            assert abs(atom.pos.z - P.residue_at(src, "A", index).find_atom(atom.name, "*").pos.z) < 1e-6
    assert P.residue_at(dst, "A", 0).find_atom("CA", "*").pos.z == P.residue_at(src, "A", 0).find_atom("CA", "*").pos.z - 1.0


def test_windows_are_positional(tmp_path):
    st = _toy(icode_at=2)   # 1, 2, 2A, 4, 5
    assert P.window(st, "A", 2, 1, 1) == [1, 2, 3]
    assert P.window(st, "A", 0, 2, 2) == [0, 1, 2]
    path = P.write_index_window(st, "A", 1, 3, str(tmp_path / "w.pdb"))
    mini = gemmi.read_structure(path)
    assert [(r.seqid.num, r.seqid.icode) for r in mini[0]["A"]] == [(2, " "), (2, "A"), (4, " ")]
    path = P.write_window(st, "A", 2, 2, str(tmp_path / "n.pdb"))
    assert len(gemmi.read_structure(path)[0]["A"]) == 2   # both residues numbered 2


def test_write_neighbourhood_can_be_pinned_to_a_selection(tmp_path):
    st = _toy()
    keys = P.neighbourhood_keys(st, [("A", 1, "")], radius=5.0)
    assert ("A", 1, "") in keys and ("A", 2, "") in keys
    assert not any(k[0] == "B" for k in keys)  # chain B is 20 A away
    pinned = {("A", 5, "")}
    path = P.write_neighbourhood(st, [("A", 1, "")], str(tmp_path / "n.pdb"), selection=pinned)
    mini = gemmi.read_structure(path)
    assert [(ch.name, r.seqid.num) for ch in mini[0] for r in ch] == [("A", 5)]


def test_sidechain_density_counts_cb_so_alanine_is_guarded():
    class Grid:
        def interpolate_value(self, pos):
            return 1.0
    ala = P.find_residue(_toy(), "A", 1)
    assert ala.name == "ALA"
    assert P.sidechain_density(Grid(), ala) == 1.0
    assert P.is_amino_acid("ALA") and not P.is_amino_acid("CL") and not P.is_amino_acid("HOH")


def _observation(keys, side, back, rot=None, ram=None, clashes=(), contacts=None):
    return dict(keys=keys, lo=0, hi=len(keys) - 1,
                density={k: (side[k], back[k]) for k in keys},
                contacts={k: (contacts or {}).get(k, set()) for k in keys},
                rot=rot or {}, ram=ram or {}, clashes=list(clashes), selection=set(keys))


def test_judge_guards_every_moved_residue_not_just_the_target():
    keys = [("A", 1, ""), ("A", 2, ""), ("A", 3, "")]
    target = keys[1]
    base = _observation(keys, {k: 2.0 for k in keys}, {k: 3.0 for k in keys},
                        rot={target: dict(evaluation="OUTLIER")}, ram={k: dict(evaluation="Favored") for k in keys})
    good = _observation(keys, {k: 1.9 for k in keys}, {k: 3.0 for k in keys},
                        rot={target: dict(evaluation="Favored", rotamer="p")},
                        ram={k: dict(evaluation="Favored") for k in keys})
    assert P.judge(base, good, target, 0.6, 0.8, 0.8).ok
    # a flank lost its backbone density: the peptide was distorted
    bad = _observation(keys, {k: 1.9 for k in keys}, {keys[0]: 0.5, keys[1]: 3.0, keys[2]: 3.0},
                       rot=good["rot"], ram=good["ram"])
    verdict = P.judge(base, bad, target, 0.6, 0.8, 0.8)
    assert not verdict.ok and "backbone density A 1" in verdict.reasons[0]
    # a flank became a Ramachandran outlier, two residues beyond what moved
    far = ("A", 5, "")
    wider = dict(good["ram"])
    wider[far] = dict(evaluation="OUTLIER")
    bad = _observation(keys, {k: 1.9 for k in keys}, {k: 3.0 for k in keys}, rot=good["rot"], ram=wider)
    verdict = P.judge(base, bad, target, 0.6, 0.8, 0.8)
    assert not verdict.ok and "new Ramachandran outlier A 5" in verdict.reasons[0]


def test_judge_vetoes_a_new_molprobity_clash_anywhere_in_the_neighbourhood():
    keys = [("A", 2, "")]
    base = _observation(keys, {keys[0]: 2.0}, {keys[0]: 3.0})
    _, clashes = P.parse_clashscore(CLASHSCORE)
    trial = _observation(keys, {keys[0]: 2.0}, {keys[0]: 3.0}, clashes=clashes[:1])
    verdict = P.judge(base, trial, keys[0], 0.6, 0.8, 0.8)
    assert not verdict.ok and verdict.reasons == ["new clash PHE C 1106 HB3 - SER C 1127A HB3 (1.07)"]
    # below the threshold it is tolerated
    assert P.judge(base, _observation(keys, {keys[0]: 2.0}, {keys[0]: 3.0}, clashes=clashes[3:]),
                   keys[0], 0.6, 0.8, 0.8).ok


def test_review_rows_keep_chain_number_and_name_together():
    _, clashes = P.parse_clashscore(CLASHSCORE)
    clash = clashes[3]      # A 12 SER OG - B 45 LYS NZ
    assert P.key_label(clash["key1"], clash["resname1"]) == "SER A 12"
    assert P.key_label(clash["key2"], clash["resname2"]) == "LYS B 45"
    assert P.key_label(("C", 1127, "A"), "SER") == "SER C 1127A"


def test_changes_csv_has_fixed_columns_with_icode(tmp_path):
    rows = [dict(chain="A", resno=12, icode="A", resname="SER", max_shift=1.2, action="rotamers:autofit",
                 density_before=1.1, density_after=1.0, rotamer_before="OUTLIER", rotamer_after="p",
                 rama_before="Favored", rama_after="Favored")]
    path = P.write_changes(str(tmp_path / "c.csv"), rows)
    with open(path) as handle:
        read = list(csv.DictReader(handle))
    assert list(read[0].keys()) == list(P.CHANGES_COLUMNS)
    assert read[0]["icode"] == "A" and read[0]["action"] == "rotamers:autofit"
    assert read[0]["diff_density_before"] == ""


def test_map_labels_and_missing_columns_are_reported():
    assert P.map_labels(["F", "PHI"]) == ("F", "PHI")
    assert P.map_labels(["FWT", "PHWT"]) == ("FWT", "PHWT")
    assert P.map_labels(["H", "K", "L"]) == (None, None)


def test_options_reject_unknown_knob():
    with pytest.raises(TypeError):
        P.Options(bogus=1)


def test_plugin_imports_without_coot_or_molprobity():
    module = importlib.import_module("ccp4i2.wrappers.geometry_polish.script.geometry_polish")
    assert module.geometry_polish.TASKNAME == "geometry_polish"
    assert "molprobity.rotalyze" in module.geometry_polish.AUXILIARY_PROGRAMS
    report = importlib.import_module("ccp4i2.wrappers.geometry_polish.script.geometry_polish_report")
    assert report.geometry_polish_report.TASKNAME == "geometry_polish"


def test_program_xml_is_written_atomically(tmp_path, monkeypatch):
    module = importlib.import_module("ccp4i2.wrappers.geometry_polish.script.geometry_polish")
    plugin = module.geometry_polish.__new__(module.geometry_polish)
    target = tmp_path / "program.xml"
    plugin.makeFileName = lambda kind: str(target)
    moves = []
    real_move = module.shutil.move
    monkeypatch.setattr(module.shutil, "move", lambda a, b: (moves.append((a, b)), real_move(a, b)))
    plugin._write_xml(dict(before=dict(rot_outliers=3), stages=[dict(name="rotamers", fixed=1)],
                           changes=[dict(chain="A", resno=5, icode="", resname="SER")], review=[]))
    assert moves == [(str(target) + "_tmp", str(target))]
    assert target.exists() and not (tmp_path / "program.xml_tmp").exists()
    text = target.read_text()
    assert 'rot_outliers="3"' in text and 'icode=""' in text


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
