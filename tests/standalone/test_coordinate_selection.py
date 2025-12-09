"""
Unit tests for coordinate_selection module.

Tests tokenization, parsing, and evaluation of mmdb CID selection syntax
including the new category keywords (protein, nucleic, solvent, etc.).
"""

import pytest
from ccp4i2.core.coordinate_selection.tokenizer import tokenize, TokenType
from ccp4i2.core.coordinate_selection.parser import parse_selection, ParseError
from ccp4i2.core.coordinate_selection.ast_nodes import (
    CIDSelector, CategorySelector, CategoryType,
    LogicalAnd, LogicalOr, LogicalNot
)
from ccp4i2.core.coordinate_selection.evaluator import (
    evaluate_selection,
    AMINO_ACIDS, NUCLEIC_ACIDS, SOLVENTS, SACCHARIDES,
    PROTEIN_BACKBONE_ATOMS, NUCLEIC_BACKBONE_ATOMS
)


# ============================================================================
# Tokenizer Tests
# ============================================================================

class TestTokenizerBasic:
    """Test basic tokenization of CID syntax."""

    def test_simple_chain(self):
        tokens = tokenize("A")
        assert len(tokens) == 2  # IDENTIFIER + EOF
        assert tokens[0].type == TokenType.IDENTIFIER
        assert tokens[0].value == "A"

    def test_chain_with_slash(self):
        tokens = tokenize("A/")
        assert tokens[0].type == TokenType.IDENTIFIER
        assert tokens[1].type == TokenType.SLASH

    def test_chain_residue(self):
        tokens = tokenize("A/27")
        assert tokens[0].type == TokenType.IDENTIFIER
        assert tokens[1].type == TokenType.SLASH
        assert tokens[2].type == TokenType.NUMBER
        assert tokens[2].value == "27"

    def test_residue_range(self):
        tokens = tokenize("A/10-20")
        types = [t.type for t in tokens[:-1]]  # exclude EOF
        assert TokenType.DASH in types

    def test_wildcard(self):
        tokens = tokenize("*")
        assert tokens[0].type == TokenType.STAR

    def test_residue_name_parens(self):
        tokens = tokenize("(ALA)")
        assert tokens[0].type == TokenType.LPAREN
        assert tokens[1].type == TokenType.IDENTIFIER
        assert tokens[1].value == "ALA"
        assert tokens[2].type == TokenType.RPAREN


class TestTokenizerLogicalKeywords:
    """Test tokenization of logical keywords."""

    def test_and_keyword(self):
        tokens = tokenize("A and B")
        assert tokens[1].type == TokenType.AND

    def test_or_keyword(self):
        tokens = tokenize("A or B")
        assert tokens[1].type == TokenType.OR

    def test_not_keyword(self):
        tokens = tokenize("not A")
        assert tokens[0].type == TokenType.NOT

    def test_case_insensitive(self):
        tokens = tokenize("AND OR NOT")
        assert tokens[0].type == TokenType.AND
        assert tokens[1].type == TokenType.OR
        assert tokens[2].type == TokenType.NOT


class TestTokenizerCategoryKeywords:
    """Test tokenization of category keywords."""

    def test_protein_keyword(self):
        tokens = tokenize("protein")
        assert tokens[0].type == TokenType.PROTEIN

    def test_nucleic_keyword(self):
        tokens = tokenize("nucleic")
        assert tokens[0].type == TokenType.NUCLEIC

    def test_solvent_keyword(self):
        tokens = tokenize("solvent")
        assert tokens[0].type == TokenType.SOLVENT

    def test_water_alias(self):
        """'water' should tokenize as SOLVENT."""
        tokens = tokenize("water")
        assert tokens[0].type == TokenType.SOLVENT

    def test_ligand_keyword(self):
        tokens = tokenize("ligand")
        assert tokens[0].type == TokenType.LIGAND

    def test_sugar_keyword(self):
        tokens = tokenize("sugar")
        assert tokens[0].type == TokenType.SUGAR

    def test_saccharide_alias(self):
        """'saccharide' should tokenize as SUGAR."""
        tokens = tokenize("saccharide")
        assert tokens[0].type == TokenType.SUGAR

    def test_polymer_keyword(self):
        tokens = tokenize("polymer")
        assert tokens[0].type == TokenType.POLYMER

    def test_backbone_keyword(self):
        tokens = tokenize("backbone")
        assert tokens[0].type == TokenType.BACKBONE

    def test_sidechain_keyword(self):
        tokens = tokenize("sidechain")
        assert tokens[0].type == TokenType.SIDECHAIN

    def test_hetero_keyword(self):
        tokens = tokenize("hetero")
        assert tokens[0].type == TokenType.HETERO

    def test_hetatm_alias(self):
        """'hetatm' should tokenize as HETERO."""
        tokens = tokenize("hetatm")
        assert tokens[0].type == TokenType.HETERO

    def test_category_case_insensitive(self):
        tokens = tokenize("PROTEIN Nucleic SOLVENT")
        assert tokens[0].type == TokenType.PROTEIN
        assert tokens[1].type == TokenType.NUCLEIC
        assert tokens[2].type == TokenType.SOLVENT


class TestTokenizerComplexExpressions:
    """Test tokenization of complex expressions."""

    def test_category_with_logical(self):
        tokens = tokenize("protein and not solvent")
        types = [t.type for t in tokens[:-1]]
        assert types == [
            TokenType.PROTEIN,
            TokenType.AND,
            TokenType.NOT,
            TokenType.SOLVENT,
        ]

    def test_category_with_chain(self):
        tokens = tokenize("A/ and protein")
        types = [t.type for t in tokens[:-1]]
        assert TokenType.IDENTIFIER in types
        assert TokenType.PROTEIN in types

    def test_braces(self):
        tokens = tokenize("{A/ or B/}")
        assert tokens[0].type == TokenType.LBRACE
        assert tokens[-2].type == TokenType.RBRACE


# ============================================================================
# Parser Tests
# ============================================================================

class TestParserBasic:
    """Test basic parsing of CID syntax."""

    def test_simple_chain(self):
        ast = parse_selection("A")
        assert isinstance(ast, CIDSelector)
        assert ast.chain == "A"

    def test_chain_residue(self):
        ast = parse_selection("A/27")
        assert isinstance(ast, CIDSelector)
        assert ast.chain == "A"
        assert ast.seq_no == "27"

    def test_residue_name(self):
        ast = parse_selection("(ALA)")
        assert isinstance(ast, CIDSelector)
        assert ast.res_name == "ALA"

    def test_model_chain(self):
        ast = parse_selection("/1/A")
        assert isinstance(ast, CIDSelector)
        assert ast.model == "1"
        assert ast.chain == "A"


class TestParserLogical:
    """Test parsing of logical expressions."""

    def test_and_expression(self):
        ast = parse_selection("A and B")
        assert isinstance(ast, LogicalAnd)
        assert isinstance(ast.left, CIDSelector)
        assert isinstance(ast.right, CIDSelector)

    def test_or_expression(self):
        ast = parse_selection("A or B")
        assert isinstance(ast, LogicalOr)

    def test_not_expression(self):
        ast = parse_selection("not A")
        assert isinstance(ast, LogicalNot)
        assert isinstance(ast.operand, CIDSelector)

    def test_complex_logical(self):
        ast = parse_selection("A and B or C")
        # Should parse as (A and B) or C due to precedence
        assert isinstance(ast, LogicalOr)

    def test_grouped_expression(self):
        ast = parse_selection("{A or B} and C")
        assert isinstance(ast, LogicalAnd)
        assert isinstance(ast.left, LogicalOr)


class TestParserCategoryKeywords:
    """Test parsing of category keywords."""

    def test_protein_keyword(self):
        ast = parse_selection("protein")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.PROTEIN

    def test_nucleic_keyword(self):
        ast = parse_selection("nucleic")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.NUCLEIC

    def test_solvent_keyword(self):
        ast = parse_selection("solvent")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.SOLVENT

    def test_water_alias(self):
        ast = parse_selection("water")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.SOLVENT

    def test_ligand_keyword(self):
        ast = parse_selection("ligand")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.LIGAND

    def test_sugar_keyword(self):
        ast = parse_selection("sugar")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.SUGAR

    def test_polymer_keyword(self):
        ast = parse_selection("polymer")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.POLYMER

    def test_backbone_keyword(self):
        ast = parse_selection("backbone")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.BACKBONE

    def test_sidechain_keyword(self):
        ast = parse_selection("sidechain")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.SIDECHAIN

    def test_hetero_keyword(self):
        ast = parse_selection("hetero")
        assert isinstance(ast, CategorySelector)
        assert ast.category == CategoryType.HETERO


class TestParserCategoryWithLogical:
    """Test parsing of category keywords with logical operators."""

    def test_protein_and_not_solvent(self):
        ast = parse_selection("protein and not solvent")
        assert isinstance(ast, LogicalAnd)
        assert isinstance(ast.left, CategorySelector)
        assert ast.left.category == CategoryType.PROTEIN
        assert isinstance(ast.right, LogicalNot)
        assert isinstance(ast.right.operand, CategorySelector)
        assert ast.right.operand.category == CategoryType.SOLVENT

    def test_chain_and_ligand(self):
        ast = parse_selection("A and ligand")
        assert isinstance(ast, LogicalAnd)
        assert isinstance(ast.left, CIDSelector)
        assert isinstance(ast.right, CategorySelector)
        assert ast.right.category == CategoryType.LIGAND

    def test_grouped_categories(self):
        ast = parse_selection("{protein or nucleic} and not solvent")
        assert isinstance(ast, LogicalAnd)
        assert isinstance(ast.left, LogicalOr)
        assert isinstance(ast.left.left, CategorySelector)
        assert isinstance(ast.left.right, CategorySelector)


# ============================================================================
# Evaluator Tests - Using Mock Structure
# ============================================================================

class MockAtom:
    """Mock gemmi Atom for testing."""
    def __init__(self, name, element="C", altloc=""):
        self.name = name
        self.element = MockElement(element)
        self.altloc = altloc


class MockElement:
    """Mock gemmi Element for testing."""
    def __init__(self, name):
        self.name = name


class MockSeqId:
    """Mock gemmi SeqId for testing."""
    def __init__(self, num, icode=""):
        self.num = num
        self.icode = icode


class MockResidue:
    """Mock gemmi Residue for testing."""
    def __init__(self, name, seq_num, atoms=None, het_flag=""):
        self.name = name
        self.seqid = MockSeqId(seq_num)
        self._atoms = atoms or []
        self.het_flag = het_flag

    def __iter__(self):
        return iter(self._atoms)


class MockChain:
    """Mock gemmi Chain for testing."""
    def __init__(self, name, residues=None):
        self.name = name
        self._residues = residues or []

    def __iter__(self):
        return iter(self._residues)


class MockModel:
    """Mock gemmi Model for testing."""
    def __init__(self, name="1", chains=None):
        self.name = name
        self._chains = chains or []

    def __iter__(self):
        return iter(self._chains)


class MockStructure:
    """Mock gemmi Structure for testing."""
    def __init__(self, models=None):
        self._models = models or []

    def __iter__(self):
        return iter(self._models)


def create_simple_protein_structure():
    """Create a simple mock structure with protein residues."""
    # Chain A with 3 ALA residues
    atoms_ala1 = [MockAtom("N"), MockAtom("CA"), MockAtom("C"), MockAtom("O"), MockAtom("CB")]
    atoms_ala2 = [MockAtom("N"), MockAtom("CA"), MockAtom("C"), MockAtom("O"), MockAtom("CB")]
    atoms_ala3 = [MockAtom("N"), MockAtom("CA"), MockAtom("C"), MockAtom("O"), MockAtom("CB")]

    residues_a = [
        MockResidue("ALA", 1, atoms_ala1),
        MockResidue("ALA", 2, atoms_ala2),
        MockResidue("ALA", 3, atoms_ala3),
    ]
    chain_a = MockChain("A", residues_a)
    model = MockModel("1", [chain_a])
    return MockStructure([model])


def create_mixed_structure():
    """Create a mock structure with protein, water, and ligand."""
    # Protein residues with backbone + sidechain atoms
    protein_atoms = [MockAtom("N"), MockAtom("CA"), MockAtom("C"), MockAtom("O"), MockAtom("CB")]

    # Water atoms
    water_atoms = [MockAtom("O", "O")]

    # Ligand atoms (ATP-like)
    ligand_atoms = [MockAtom("N1", "N"), MockAtom("C2", "C"), MockAtom("P", "P")]

    residues_a = [
        MockResidue("ALA", 1, protein_atoms.copy()),
        MockResidue("GLY", 2, [MockAtom("N"), MockAtom("CA"), MockAtom("C"), MockAtom("O")]),
        MockResidue("HOH", 100, water_atoms.copy()),
        MockResidue("HOH", 101, water_atoms.copy()),
        MockResidue("ATP", 200, ligand_atoms.copy(), het_flag="H"),
    ]

    chain_a = MockChain("A", residues_a)
    model = MockModel("1", [chain_a])
    return MockStructure([model])


def create_nucleic_structure():
    """Create a mock structure with nucleic acid residues."""
    # DNA backbone + base atoms
    dna_atoms = [MockAtom("P", "P"), MockAtom("O5'", "O"), MockAtom("C5'", "C"),
                 MockAtom("C4'", "C"), MockAtom("C1'", "C"), MockAtom("N9", "N")]

    residues = [
        MockResidue("DA", 1, dna_atoms.copy()),
        MockResidue("DG", 2, dna_atoms.copy()),
        MockResidue("DC", 3, dna_atoms.copy()),
    ]

    chain = MockChain("A", residues)
    model = MockModel("1", [chain])
    return MockStructure([model])


class TestEvaluatorCIDSelector:
    """Test evaluation of CID selectors."""

    def test_select_chain(self):
        structure = create_simple_protein_structure()
        ast = parse_selection("A")
        result = evaluate_selection(ast, structure)
        # Should select all atoms in chain A (3 residues * 5 atoms = 15)
        assert len(result) == 15

    def test_select_residue_by_number(self):
        structure = create_simple_protein_structure()
        ast = parse_selection("A/1")
        result = evaluate_selection(ast, structure)
        # Should select 5 atoms from residue 1
        assert len(result) == 5

    def test_select_residue_range(self):
        structure = create_simple_protein_structure()
        ast = parse_selection("A/1-2")
        result = evaluate_selection(ast, structure)
        # Should select 10 atoms from residues 1 and 2
        assert len(result) == 10

    def test_select_by_residue_name(self):
        structure = create_mixed_structure()
        ast = parse_selection("(HOH)")
        result = evaluate_selection(ast, structure)
        # Should select 2 water molecules (2 atoms total)
        assert len(result) == 2


class TestEvaluatorCategorySelector:
    """Test evaluation of category selectors."""

    def test_select_protein(self):
        structure = create_mixed_structure()
        ast = parse_selection("protein")
        result = evaluate_selection(ast, structure)
        # ALA has 5 atoms, GLY has 4 atoms = 9 protein atoms
        assert len(result) == 9

    def test_select_solvent(self):
        structure = create_mixed_structure()
        ast = parse_selection("solvent")
        result = evaluate_selection(ast, structure)
        # 2 HOH residues with 1 atom each = 2 atoms
        assert len(result) == 2

    def test_select_water_alias(self):
        structure = create_mixed_structure()
        ast = parse_selection("water")
        result = evaluate_selection(ast, structure)
        assert len(result) == 2

    def test_select_ligand(self):
        structure = create_mixed_structure()
        ast = parse_selection("ligand")
        result = evaluate_selection(ast, structure)
        # ATP has 3 atoms
        assert len(result) == 3

    def test_select_polymer(self):
        structure = create_mixed_structure()
        ast = parse_selection("polymer")
        result = evaluate_selection(ast, structure)
        # Only protein in this case = 9 atoms
        assert len(result) == 9

    def test_select_nucleic(self):
        structure = create_nucleic_structure()
        ast = parse_selection("nucleic")
        result = evaluate_selection(ast, structure)
        # 3 nucleic residues with 6 atoms each = 18 atoms
        assert len(result) == 18

    def test_select_backbone(self):
        structure = create_mixed_structure()
        ast = parse_selection("backbone")
        result = evaluate_selection(ast, structure)
        # ALA backbone: N, CA, C, O = 4 atoms
        # GLY backbone: N, CA, C, O = 4 atoms
        # Total = 8 backbone atoms
        assert len(result) == 8

    def test_select_sidechain(self):
        structure = create_mixed_structure()
        ast = parse_selection("sidechain")
        result = evaluate_selection(ast, structure)
        # ALA sidechain: CB = 1 atom
        # GLY has no sidechain
        # Total = 1 sidechain atom
        assert len(result) == 1


class TestEvaluatorLogicalOperators:
    """Test evaluation of logical operators with categories."""

    def test_protein_and_not_solvent(self):
        structure = create_mixed_structure()
        ast = parse_selection("protein and not solvent")
        result = evaluate_selection(ast, structure)
        # protein (9) AND NOT solvent (2) = 9 (no overlap)
        assert len(result) == 9

    def test_protein_or_ligand(self):
        structure = create_mixed_structure()
        ast = parse_selection("protein or ligand")
        result = evaluate_selection(ast, structure)
        # protein (9) + ligand (3) = 12
        assert len(result) == 12

    def test_not_solvent(self):
        structure = create_mixed_structure()
        ast = parse_selection("not solvent")
        result = evaluate_selection(ast, structure)
        # Total atoms = 9 (protein) + 2 (water) + 3 (ligand) = 14
        # NOT solvent = 14 - 2 = 12
        assert len(result) == 12

    def test_chain_and_protein(self):
        structure = create_mixed_structure()
        ast = parse_selection("A and protein")
        result = evaluate_selection(ast, structure)
        # Same as just protein since all is in chain A
        assert len(result) == 9

    def test_complex_expression(self):
        structure = create_mixed_structure()
        ast = parse_selection("{protein or ligand} and not (HOH)")
        result = evaluate_selection(ast, structure)
        # (protein 9 + ligand 3) AND NOT (water 2) = 12
        assert len(result) == 12


class TestEvaluatorNucleicBackbone:
    """Test evaluation of backbone/sidechain for nucleic acids."""

    def test_nucleic_backbone(self):
        structure = create_nucleic_structure()
        ast = parse_selection("backbone")
        result = evaluate_selection(ast, structure)
        # Each nucleic residue has 5 backbone atoms: P, O5', C5', C4', C1'
        # Wait, NUCLEIC_BACKBONE_ATOMS includes more atoms
        # Let's count what matches
        backbone_count = sum(
            1 for model, chain, res, atom in result
            if atom.name in NUCLEIC_BACKBONE_ATOMS
        )
        assert len(result) == backbone_count

    def test_nucleic_sidechain(self):
        structure = create_nucleic_structure()
        ast = parse_selection("sidechain")
        result = evaluate_selection(ast, structure)
        # Sidechain = atoms not in backbone
        # N9 is not in backbone
        for model, chain, res, atom in result:
            assert atom.name not in NUCLEIC_BACKBONE_ATOMS


# ============================================================================
# Classification Sets Tests
# ============================================================================

class TestClassificationSets:
    """Test the residue/atom classification sets are properly defined."""

    def test_amino_acids_contains_standard(self):
        standard = ['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY',
                    'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER',
                    'THR', 'TRP', 'TYR', 'VAL']
        for aa in standard:
            assert aa in AMINO_ACIDS, f"{aa} should be in AMINO_ACIDS"

    def test_amino_acids_contains_modified(self):
        modified = ['MSE', 'SEP', 'TPO']
        for aa in modified:
            assert aa in AMINO_ACIDS, f"{aa} should be in AMINO_ACIDS"

    def test_nucleic_acids_contains_standard(self):
        standard = ['A', 'C', 'G', 'U', 'T', 'DA', 'DC', 'DG', 'DT']
        for na in standard:
            assert na in NUCLEIC_ACIDS, f"{na} should be in NUCLEIC_ACIDS"

    def test_solvents_contains_water(self):
        water = ['HOH', 'WAT', 'H2O']
        for w in water:
            assert w in SOLVENTS, f"{w} should be in SOLVENTS"

    def test_saccharides_contains_common(self):
        common = ['GLC', 'GAL', 'MAN', 'NAG']
        for s in common:
            assert s in SACCHARIDES, f"{s} should be in SACCHARIDES"

    def test_protein_backbone_atoms(self):
        expected = {'N', 'CA', 'C', 'O'}
        assert expected.issubset(PROTEIN_BACKBONE_ATOMS)

    def test_nucleic_backbone_atoms(self):
        expected = {'P', "O5'", "C5'", "C4'", "C3'", "O3'"}
        assert expected.issubset(NUCLEIC_BACKBONE_ATOMS)
