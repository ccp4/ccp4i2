# Copyright (C) 2026 Newcastle University
#
# This file is part of CCP4i2.
#
# CCP4i2 is free software: you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License version 3,
# modified in accordance with the provisions of the license to address
# the requirements of UK law.
#
# See https://www.ccp4.ac.uk/ccp4license.php for details.
"""Tests for the checkMonomeCoverage pre-flight validation.

These tests exercise the CPluginScript.checkMonomeCoverage method directly
using synthetic structures and dictionaries, without running any refinement.
"""
import os
import tempfile
import textwrap

import gemmi
import pytest

from ccp4i2.core.CCP4PluginScript import CPluginScript


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _make_structure(residues):
    """Build a minimal gemmi.Structure and write it to a temp mmCIF file.

    Args:
        residues: list of (chain_id, res_name, seq_num, entity_type, atoms)
            where atoms is a list of (atom_name, element_symbol) tuples.
    Returns:
        path (str) to a temporary mmCIF file (caller must delete).
    """
    st = gemmi.Structure()
    st.name = "test"
    st.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
    st.spacegroup_hm = "P 1"
    model = gemmi.Model("1")

    chains = {}
    for chain_id, res_name, seq_num, etype, atoms in residues:
        if chain_id not in chains:
            chains[chain_id] = gemmi.Chain(chain_id)
        res = gemmi.Residue()
        res.name = res_name
        res.seqid = gemmi.SeqId(str(seq_num))
        res.entity_type = etype
        x = 0.0
        for aname, elem in atoms:
            atom = gemmi.Atom()
            atom.name = aname
            atom.element = gemmi.Element(elem)
            atom.pos = gemmi.Position(x, 0, 0)
            atom.occ = 1.0
            atom.b_iso = 20.0
            x += 1.5
            res.add_atom(atom)
        chains[chain_id].add_residue(res)

    for ch in chains.values():
        model.add_chain(ch)
    st.add_model(model)

    fd, path = tempfile.mkstemp(suffix=".cif")
    os.close(fd)
    st.make_mmcif_document().write_file(path)
    return path


def _make_dict(monomer_code, atoms, group="non-polymer"):
    """Write a minimal CIF restraint dictionary defining one monomer.

    Args:
        monomer_code: three-letter (or longer) code.
        atoms: list of (atom_name, element_symbol) tuples.
        group: chem_comp group string.
    Returns:
        path (str) to a temporary CIF dictionary file (caller must delete).
    """
    atom_lines = "\n".join(
        f"{monomer_code} {name:5s} {elem} {elem} 0.000"
        for name, elem in atoms
    )
    content = textwrap.dedent(f"""\
        data_comp_list
        loop_
        _chem_comp.id
        _chem_comp.three_letter_code
        _chem_comp.name
        _chem_comp.group
        _chem_comp.number_atoms_all
        _chem_comp.number_atoms_nh
        {monomer_code} {monomer_code} "test monomer" {group} {len(atoms)} {len(atoms)}

        data_comp_{monomer_code}
        loop_
        _chem_comp_atom.comp_id
        _chem_comp_atom.atom_id
        _chem_comp_atom.type_symbol
        _chem_comp_atom.type_energy
        _chem_comp_atom.charge
        {atom_lines}
    """)
    fd, path = tempfile.mkstemp(suffix=".cif")
    with os.fdopen(fd, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def plugin():
    """Create a minimal CPluginScript instance for testing."""
    p = CPluginScript.__new__(CPluginScript)
    # Minimal initialisation needed by appendErrorReport / checkMonomeCoverage
    from ccp4i2.core.base_object.error_reporting import CErrorReport
    p.errorReport = CErrorReport()
    return p


# ---------------------------------------------------------------------------
# tests
# ---------------------------------------------------------------------------

class TestCheckMonomeCoverage:

    def test_standard_amino_acids_only(self, plugin):
        """Structure with only standard amino acids — should always pass."""
        path = _make_structure([
            ("A", "ALA", 1, gemmi.EntityType.Polymer,
             [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")]),
            ("A", "GLY", 2, gemmi.EntityType.Polymer,
             [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O")]),
        ])
        try:
            assert plugin.checkMonomeCoverage(path) == CPluginScript.SUCCEEDED
        finally:
            os.unlink(path)

    def test_known_library_ligand(self, plugin):
        """GOL (glycerol) with correct atoms — should pass via CCP4 monomer library."""
        # GOL in the CCP4 library has: C1 O1 C2 O2 C3 O3
        path = _make_structure([
            ("A", "ALA", 1, gemmi.EntityType.Polymer,
             [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")]),
            ("A", "GOL", 2, gemmi.EntityType.NonPolymer,
             [("C1", "C"), ("O1", "O"), ("C2", "C"), ("O2", "O"),
              ("C3", "C"), ("O3", "O")]),
        ])
        try:
            assert plugin.checkMonomeCoverage(path) == CPluginScript.SUCCEEDED
        finally:
            os.unlink(path)

    def test_code_collision_detected(self, plugin):
        """GOL code but with drug-like atoms — should FAIL (code collision)."""
        path = _make_structure([
            ("A", "GOL", 1, gemmi.EntityType.NonPolymer,
             [("C1", "C"), ("C2", "C"), ("C3", "C"), ("C4", "C"),
              ("N1", "N"), ("S1", "S"), ("F1", "F")]),
        ])
        try:
            rv = plugin.checkMonomeCoverage(path)
            assert rv == CPluginScript.FAILED
            # Error report should mention uncovered atoms
            assert plugin.errorReport.count() > 0
        finally:
            os.unlink(path)

    def test_user_dict_resolves_collision(self, plugin):
        """User-provided dict with correct atoms overrides library — should pass."""
        struct_atoms = [("C1", "C"), ("C2", "C"), ("N1", "N"), ("S1", "S")]
        path = _make_structure([
            ("A", "GOL", 1, gemmi.EntityType.NonPolymer, struct_atoms),
        ])
        dict_path = _make_dict("GOL", struct_atoms)
        try:
            rv = plugin.checkMonomeCoverage(path, [dict_path])
            assert rv == CPluginScript.SUCCEEDED
        finally:
            os.unlink(path)
            os.unlink(dict_path)

    def test_completely_unknown_code(self, plugin):
        """Residue code not in library and no user dict — should FAIL."""
        path = _make_structure([
            ("A", "QQQ", 1, gemmi.EntityType.NonPolymer,
             [("C1", "C"), ("N1", "N"), ("O1", "O")]),
        ])
        try:
            rv = plugin.checkMonomeCoverage(path)
            assert rv == CPluginScript.FAILED
        finally:
            os.unlink(path)

    def test_user_dict_covers_unknown_code(self, plugin):
        """User dict provides definition for unknown code — should pass."""
        struct_atoms = [("C1", "C"), ("N1", "N"), ("O1", "O")]
        path = _make_structure([
            ("A", "QQQ", 1, gemmi.EntityType.NonPolymer, struct_atoms),
        ])
        dict_path = _make_dict("QQQ", struct_atoms)
        try:
            rv = plugin.checkMonomeCoverage(path, [dict_path])
            assert rv == CPluginScript.SUCCEEDED
        finally:
            os.unlink(path)
            os.unlink(dict_path)

    def test_partial_dict_coverage_fails(self, plugin):
        """User dict defines some atoms but not all — should FAIL."""
        path = _make_structure([
            ("A", "QQQ", 1, gemmi.EntityType.NonPolymer,
             [("C1", "C"), ("N1", "N"), ("O1", "O"), ("S1", "S")]),
        ])
        # Dict only defines C1 and N1, missing O1 and S1
        dict_path = _make_dict("QQQ", [("C1", "C"), ("N1", "N")])
        try:
            rv = plugin.checkMonomeCoverage(path, [dict_path])
            assert rv == CPluginScript.FAILED
        finally:
            os.unlink(path)
            os.unlink(dict_path)

    def test_water_and_ions_pass(self, plugin):
        """Water molecules and common ions should pass without dictionaries."""
        path = _make_structure([
            ("A", "ALA", 1, gemmi.EntityType.Polymer,
             [("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C")]),
            ("A", "HOH", 100, gemmi.EntityType.Water,
             [("O", "O")]),
        ])
        try:
            assert plugin.checkMonomeCoverage(path) == CPluginScript.SUCCEEDED
        finally:
            os.unlink(path)

    def test_multiple_instances_all_checked(self, plugin):
        """Multiple instances of same ligand in different chains — all must pass."""
        path = _make_structure([
            # Chain A: GOL with correct atoms
            ("A", "GOL", 1, gemmi.EntityType.NonPolymer,
             [("C1", "C"), ("O1", "O"), ("C2", "C"), ("O2", "O"),
              ("C3", "C"), ("O3", "O")]),
            # Chain B: GOL with an extra wrong atom
            ("B", "GOL", 1, gemmi.EntityType.NonPolymer,
             [("C1", "C"), ("O1", "O"), ("C2", "C"), ("O2", "O"),
              ("C3", "C"), ("O3", "O"), ("BR1", "Br")]),
        ])
        try:
            rv = plugin.checkMonomeCoverage(path)
            assert rv == CPluginScript.FAILED
        finally:
            os.unlink(path)

    def test_empty_structure_passes(self, plugin):
        """Empty structure should pass (nothing to check)."""
        st = gemmi.Structure()
        st.name = "empty"
        st.cell = gemmi.UnitCell(50, 50, 50, 90, 90, 90)
        st.spacegroup_hm = "P 1"
        st.add_model(gemmi.Model("1"))
        fd, path = tempfile.mkstemp(suffix=".cif")
        os.close(fd)
        st.make_mmcif_document().write_file(path)
        try:
            assert plugin.checkMonomeCoverage(path) == CPluginScript.SUCCEEDED
        finally:
            os.unlink(path)
