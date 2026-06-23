"""
Regression test for atom-name whitespace handling in CPdbData.writeSelection().

Some upstream writers (notably older Moorhen builds that serialised coot/mmdb
molecules through gemmi) emit mmCIF atom names carrying PDB column padding, e.g.
the name is literally ' CA ' / ' N  ' rather than the conformant bare 'CA'/'N'.
gemmi preserves those spaces on read; if writeSelection re-emits them verbatim
the mmCIF contains quoted ' CA ' values that fail monomer-dictionary matching in
refinement (every residue, standard amino acids included, reported as having
"atoms without restraints").

writeSelection must therefore strip leading/trailing whitespace from atom names.
It must NOT, however, damage genuinely special-but-valid names such as the
ribose 2'-oxygen O2' (a prime is part of the name, handled by CIF quoting).
"""

import pytest

import gemmi

from ccp4i2.core.CCP4ModelData import CPdbData


def _write_padded_mmcif(path):
    """Build an mmCIF whose atom names carry PDB-style column padding."""
    st = gemmi.Structure()
    st.spacegroup_hm = "P 1"
    st.cell = gemmi.UnitCell(30, 30, 30, 90, 90, 90)
    model = gemmi.Model(1)
    chain = gemmi.Chain("A")

    res = gemmi.Residue()
    res.name = "ALA"
    res.seqid = gemmi.SeqId("1")
    for nm, el in [(" N  ", "N"), (" CA ", "C"), (" CB ", "C"), (" C  ", "C"), (" O  ", "O")]:
        a = gemmi.Atom()
        a.name = nm           # padded, as a misbehaving writer would store it
        a.element = gemmi.Element(el)
        a.pos = gemmi.Position(1, 2, 3)
        a.occ = 1.0
        a.b_iso = 20.0
        res.add_atom(a)
    chain.add_residue(res)

    # A nucleotide-style residue to prove primes survive (must NOT be stripped).
    res2 = gemmi.Residue()
    res2.name = "A"
    res2.seqid = gemmi.SeqId("2")
    for nm, el in [("O2'", "O"), ("C1'", "C")]:
        a = gemmi.Atom()
        a.name = nm
        a.element = gemmi.Element(el)
        a.pos = gemmi.Position(4, 5, 6)
        a.occ = 1.0
        a.b_iso = 20.0
        res2.add_atom(a)
    chain.add_residue(res2)

    model.add_chain(chain)
    st.add_model(model)
    st.make_mmcif_document().write_file(str(path))


def _all_atom_names(path):
    st = gemmi.read_structure(str(path))
    return [at.name for ch in st[0] for res in ch for at in res]


@pytest.mark.parametrize("out_ext", [".cif", ".pdb"])
def test_writeselection_trims_padded_atom_names(tmp_path, out_ext):
    src = tmp_path / "padded.cif"
    _write_padded_mmcif(src)

    # Sanity: the fixture really does carry padded names on read.
    assert any(n != n.strip() for n in _all_atom_names(src)), \
        "fixture should contain padded atom names"

    pdb = CPdbData()
    err = pdb.loadFile(str(src))
    assert err.count() == 0, f"loadFile errored: {err}"

    n_atoms, atoms = pdb.interpretSelection("A/")
    assert n_atoms == 7

    out = tmp_path / f"selected{out_ext}"
    rc = pdb.writeSelection(atoms, str(out))
    assert rc == 0
    assert out.exists()

    names = _all_atom_names(out)
    # No name may retain leading/trailing whitespace.
    assert all(n == n.strip() for n in names), f"untrimmed names survived: {names!r}"
    # The standard backbone atoms are now bare.
    assert {"N", "CA", "CB", "C", "O"}.issubset(set(names))
    # The prime-bearing names are preserved exactly (not stripped/mangled).
    assert "O2'" in names
    assert "C1'" in names


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
