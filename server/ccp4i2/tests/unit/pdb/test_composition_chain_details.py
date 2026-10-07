"""
Per-chain polymer / water / other counts in a coordinate file's digest.

The file Digest shows one row per chain.  A chain's whole-residue range runs
over its waters and ligands too, so a chain "A" of residues 1-298 with waters
numbered 1-38 after them reported its range as "1-38" or worse.  The digest
therefore carries the polymer residues' own count and range beside the counts
of waters and of everything else (ligands, ions, buffer) -- with or without
entity records in the file (#681).
"""

import pytest

gemmi = pytest.importorskip("gemmi", reason="needs gemmi (CCP4 python)")

from ccp4i2.core.CCP4ModelData import CPdbDataComposition


def _atom(serial, name, resname, chain, seq, el, het=False):
    rec = "HETATM" if het else "ATOM  "
    return (
        f"{rec}{serial:5d} {name:<4} {resname:>3} {chain}{seq:4d}    "
        f"{1.0 + serial:8.3f}{2.0:8.3f}{3.0:8.3f}  1.00 20.00          {el:>2}"
    )


def _pdb():
    lines = ["CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1"]
    n = 0
    for seq in (3, 4, 5):
        for nm, el in (("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")):
            n += 1
            lines.append(_atom(n, nm, "ALA", "A", seq, el))
    lines.append("TER")
    for nm, el in (("S", "S"), ("O1", "O"), ("O2", "O"), ("O3", "O"), ("O4", "O")):
        n += 1
        lines.append(_atom(n, nm, "SO4", "A", 401, el, het=True))
    for seq in (1, 2):
        n += 1
        lines.append(_atom(n, "O", "HOH", "A", seq, "O", het=True))
    lines.append("END")
    return "\n".join(lines) + "\n"


@pytest.mark.parametrize("entities", [False, True], ids=["no-entities", "entities"])
def test_chain_detail_splits_polymer_water_other(entities):
    st = gemmi.read_pdb_string(_pdb())
    if entities:
        st.setup_entities()
    (detail,) = CPdbDataComposition(st).chainDetails

    assert detail["id"] == "A"
    assert detail["type"] == "protein"
    assert detail["nResidues"] == 6
    assert detail["nPolymer"] == 3
    assert (detail["polymerFirstRes"], detail["polymerLastRes"]) == ("3", "5")
    assert detail["nWater"] == 2
    assert detail["nOther"] == 1
    # The whole-chain range runs into the waters: why the polymer range exists.
    assert detail["lastRes"] == "2"


def test_water_only_chain_has_no_polymer_range():
    pdb = "\n".join(
        [
            "CRYST1   50.000   50.000   50.000  90.00  90.00  90.00 P 1",
            _atom(1, "O", "HOH", "W", 1, "O", het=True),
            _atom(2, "O", "HOH", "W", 2, "O", het=True),
            "END",
        ]
    ) + "\n"
    (detail,) = CPdbDataComposition(gemmi.read_pdb_string(pdb)).chainDetails
    assert detail["nPolymer"] == 0
    assert detail["polymerFirstRes"] == "" and detail["polymerLastRes"] == ""
    assert detail["nWater"] == 2
    assert detail["nOther"] == 0
