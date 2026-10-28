"""Servalcat's geometry outliers as places to look at in Moorhen.

The Servalcat author asked (alpha tracker, 2026-10-16) for a way to point at
the outliers in a viewer. The report now offers the worst as buttons; each
builds a scene of the refined model and maps centred on the outlier's atoms.
This checks the atom selections, the choice of the worst, the labels, and
that every scene validates against the ccp4i2 scene schema.

The program.xml below is Servalcat's own shape, trimmed from a real 1tf6 run.
CCP4-free.
"""
import json
import os
import xml.etree.ElementTree as ET

import pytest

import ccp4i2
from ccp4i2.pipelines.servalcat_pipe.script import outlier_places as op

PROGRAM_XML = """<SERVALCAT><SERVALCAT_FIRST>
<cycle><Ncyc>0</Ncyc></cycle>
<cycle><Ncyc>1</Ncyc><geom><outliers>
  <bond><atom1>A/HIS 188/NE2</atom1><atom2>A/ZN 206/ZN</atom2><z>16.2</z></bond>
  <chir><atomc>D/VAL 173/CA</atomc><atom1>D/VAL 173/N</atom1><atom2>D/VAL 173/C</atom2><atom3>D/VAL 173/CB</atom3><z>-25.2</z></chir>
  <plane><label>plan-1</label><atom>E/DG 25/C1'</atom><z>-6.2</z></plane>
  <vdw><atom1>G/DC 2012/OP2.A</atom1><atom2>G/DC 2013A/H6</atom2><z>-5.3</z></vdw>
  <vdw><atom1>A/THR 86/CG2</atom1><atom2>A/ASN 89/HB3</atom2><z>-6.9</z></vdw>
</outliers></geom></cycle>
</SERVALCAT_FIRST></SERVALCAT>"""


@pytest.fixture
def root():
    return ET.fromstring(PROGRAM_XML)


def test_atom_labels_become_coot_selections():
    assert op.atom_cid("A/HIS 188/NE2") == "//A/188/NE2"
    assert op.atom_cid("G/DC 2012/OP2.A") == "//G/2012/OP2"      # altloc dropped
    assert op.atom_cid("G/DC 2013A/H6") == "//G/2013.A/H6"       # insertion code kept
    assert op.atom_cid("E/DG 25/C1'") == "//E/25/C1'"
    assert op.atom_cid("A/ALA -3/CA") == "//A/-3/CA"
    assert op.atom_cid("not an atom") is None


def test_every_kind_is_read_with_all_its_atoms(root):
    found = {(o["kind"], o["index"]): o["atoms"] for o in op.outliers(root)}
    assert found[("chir", 0)] == ["D/VAL 173/CA", "D/VAL 173/N", "D/VAL 173/C", "D/VAL 173/CB"]
    assert found[("plane", 0)] == ["E/DG 25/C1'"]
    assert found[("vdw", 1)] == ["A/THR 86/CG2", "A/ASN 89/HB3"]


def test_the_worst_come_first_whatever_their_kind(root):
    assert [(o["kind"], o["index"]) for o in op.worst(root, limit=3)] == [
        ("chir", 0), ("bond", 0), ("vdw", 1)]


def test_the_after_waters_refinement_is_the_one_read():
    root = ET.fromstring(PROGRAM_XML.replace(
        "</SERVALCAT_FIRST>",
        "</SERVALCAT_FIRST><SERVALCAT_WATERS><cycle><Ncyc>0</Ncyc></cycle><cycle><geom><outliers>"
        "<bond><atom1>B/GLY 1/N</atom1><atom2>B/GLY 1/CA</atom2><z>7.0</z></bond>"
        "</outliers></geom></cycle></SERVALCAT_WATERS>"))
    assert [o["atoms"] for o in op.outliers(root)] == [["B/GLY 1/N", "B/GLY 1/CA"]]


def test_labels_name_a_single_residue_once(root):
    by_kind = {o["kind"]: o for o in op.outliers(root)}
    assert op.describe(by_kind["chir"]) == "Chirality: D/VAL 173 (CA, N, C, CB) (Z -25.2)"
    assert op.describe(by_kind["bond"]) == "Bond: A/HIS 188/NE2 – A/ZN 206/ZN (Z +16.2)"


def test_the_scene_centres_on_the_outlier_with_the_maps(root):
    outlier = op.find(root, "bond", 0)
    scene = op.outlier_scene(outlier, {"XYZOUT": 11, "FPHIOUT": 12, "DIFFPHIOUT": 13}, "uuid")
    assert scene["view"]["centre"]["selection"] == "//A/188/NE2||//A/206/ZN"
    assert scene["view"]["zoom"] < 1  # closer than Moorhen's default, which shows a crystal
    # Marked by its residues as fatter sticks and its atoms as fatter balls coloured
    # whatever their element; never spheres, which bury the residue and its density
    reps = scene["elements"][0]["representations"]
    assert [r.get("selection") for r in reps] == [None, "//A/188", "//A/206", "//A/188/NE2", "//A/206/ZN"]
    assert all(r["style"] == "CBs" for r in reps)
    assert reps[1]["colour"] == op.MARKER_COLOUR and reps[1]["geometry"] == op.RESIDUE_GEOMETRY
    assert reps[3]["colour"]["raw"]["applyColourToNonCarbonAtoms"] is True
    assert reps[3]["geometry"]["ballRadius"] > reps[1]["geometry"]["ballRadius"] > 0.2
    assert [m["name"] for m in scene["maps"]] == ["2mFo-DFc", "mFo-DFc"]
    assert scene["maps"][1]["isDifference"] is True
    assert scene["activeMap"] == "2mFo-DFc"  # Moorhen's refinement tools need an active map
    # No maps recorded: the model alone
    bare = op.outlier_scene(outlier, {"XYZOUT": 11}, "uuid")
    assert bare["maps"] == [] and "activeMap" not in bare


CLIENT_SCHEMA = os.path.join(os.path.dirname(ccp4i2.__file__), "..", "..", "client", "renderer",
                             "lib", "scene", "moorhen-scene.ccp4i2.v1.json")


def test_every_scene_validates(root):
    jsonschema = pytest.importorskip("jsonschema", reason="needs jsonschema")
    if not os.path.isfile(CLIENT_SCHEMA):
        pytest.skip("client scene schema not beside this checkout")
    with open(CLIENT_SCHEMA) as handle:
        schema = json.load(handle)
    for outlier in op.outliers(root):
        jsonschema.validate(op.outlier_scene(outlier, {"XYZOUT": 1, "FPHIOUT": 2, "DIFFPHIOUT": 3}, "u"),
                            schema)


# MolProbity's block (validate_protein, nested under Validation) and the ADP
# analysis, in their real shapes
MORE_XML = """<SERVALCAT><Validation><Molprobity>
  <Summary><Ramachandran_outliers>10.92%</Ramachandran_outliers></Summary>
  <Ramachandran_outliers>
    <Outlier chain="A" seqnum="16" name="SER" score="0.046"/>
    <Outlier chain="A" seqnum="19" name="ASP" score="0.017"/>
  </Ramachandran_outliers>
  <Nonplanar_omegas><Outlier chain="A" seqnum="63" name="HIS" score="N/A"/></Nonplanar_omegas>
  <Rotamer_outliers><Outlier chain="B" seqnum=" 100A" name="TYR" score="0.132"/></Rotamer_outliers>
  <CBeta_outliers>
    <Outlier chain="A" seqnum="46" name="LYS" score="0.308"/>
    <Outlier chain="A" seqnum="63" name="HIS" score="0.409"/>
  </CBeta_outliers>
  <Side_chain_flips><Outlier chain="D" seqnum="95" name="ASN" score="N/A"/></Side_chain_flips>
  <Clashes>
    <Outlier first_atom="A 86 THR CG2" second_atom="A 89 ASN HB3" overlap="-0.6"/>
    <Outlier first_atom="A 187 CYS O" second_atom="A 188A HIS ND1" overlap="-1.305"/>
    <Outlier first_atom="garbled" second_atom="A 1 GLY N" overlap="-2.0"/>
  </Clashes>
</Molprobity></Validation>
<ADP_ANALYSIS><outliers>
  <adp_limit_low>5.00</adp_limit_low><adp_limit_high>100.00</adp_limit_high><iqr_factor>2.00</iqr_factor>
  <low><data><adp>4.00</adp><atom>A/GLY 3/CA</atom></data></low>
  <high>
    <data><adp>150.00</adp><atom>A/LYS 46/NZ</atom></data>
    <data><adp>120.00</adp><atom>A/LYS 46/CE.B</atom></data>
  </high>
</outliers></ADP_ANALYSIS></SERVALCAT>"""


@pytest.fixture
def more():
    return ET.fromstring(MORE_XML)


def test_molprobity_residues_come_section_by_section_worst_first(more):
    residue_places = [p for p in op.molprobity_places(more) if p["kind"] != "clash"]
    assert [(p["kind"], p["index"]) for p in residue_places] == [
        ("rama", 1), ("rama", 0),          # lower probability is worse
        ("rotamer", 0),
        ("cbeta", 1), ("cbeta", 0),        # larger deviation is worse
        ("omega", 0), ("flip", 0)]
    assert residue_places[2]["residues"] == ["//B/100.A"]   # insertion code, padding stripped
    assert op.describe(residue_places[0]) == "Ramachandran: A/ASP 19 (p 0.017)"
    assert op.describe(residue_places[3]) == "Cβ deviation: A/HIS 63 (0.41 Å)"
    assert op.describe(residue_places[5]) == "Non-planar ω: A/HIS 63"
    zero = ET.fromstring('<X><Molprobity><Rotamer_outliers><Outlier chain="A" seqnum="5" name="LYS" '
                         'score="0.0"/></Rotamer_outliers></Molprobity></X>')
    assert op.describe(op.molprobity_places(zero)[0]) == "Rotamer: A/LYS 5 (p < 0.001)"


def test_clashes_name_both_atoms_worst_overlap_first(more):
    clashes = [p for p in op.molprobity_places(more) if p["kind"] == "clash"]
    assert [p["index"] for p in clashes] == [1, 0]           # the garbled one is skipped
    assert clashes[0]["atoms_cids"] == ["//A/187/O", "//A/188.A/ND1"]
    assert clashes[0]["residues"] == ["//A/187", "//A/188.A"]
    assert op.describe(clashes[0]) == "Clash: A/CYS 187/O – A/HIS 188A/ND1 (1.30 Å)"


def test_adp_outliers_furthest_outside_first(more):
    found = op.adp_places(more)
    assert [(p["kind"], p["index"]) for p in found] == [("adp_high", 0), ("adp_high", 1), ("adp_low", 0)]
    assert found[1]["atoms_cids"] == ["//A/46/CE"]
    assert op.describe(found[0]) == "B high: A/LYS 46/NZ (B 150.0, limit 100.0)"


def test_a_residue_outlier_centres_on_its_residue(more):
    scene = op.outlier_scene(op.find(more, "rotamer", 0), {"XYZOUT": 1}, "u")
    assert scene["view"]["centre"]["selection"] == "//B/100.A"
    assert [r.get("selection") for r in scene["elements"][0]["representations"]] == [None, "//B/100.A"]


def test_a_clash_centres_on_its_atoms(more):
    scene = op.outlier_scene(op.find(more, "clash", 1), {"XYZOUT": 1}, "u")
    assert scene["view"]["centre"]["selection"] == "//A/187/O||//A/188.A/ND1"


def test_no_molprobity_or_adp_block_means_no_places(root):
    assert op.molprobity_places(root) == [] and op.adp_places(root) == []


def test_every_other_scene_validates(more):
    jsonschema = pytest.importorskip("jsonschema", reason="needs jsonschema")
    if not os.path.isfile(CLIENT_SCHEMA):
        pytest.skip("client scene schema not beside this checkout")
    with open(CLIENT_SCHEMA) as handle:
        schema = json.load(handle)
    for place in op.places(more):
        jsonschema.validate(op.outlier_scene(place, {"XYZOUT": 1, "FPHIOUT": 2}, "u"), schema)
