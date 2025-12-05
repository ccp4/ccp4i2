"""
Evaluator for coordinate selection AST using gemmi.

Applies selection criteria to gemmi.Structure objects.
"""

from typing import List, Set, Tuple
import re
from .ast_nodes import (
    CIDSelector, CategorySelector, CategoryType,
    LogicalAnd, LogicalOr, LogicalNot, SelectionNode
)


# ============================================================================
# Residue and Atom Classification Sets
# ============================================================================

# Standard amino acid residue names (including common modified)
AMINO_ACIDS = {
    'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE',
    'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'MSE', 'SEP', 'TPO', 'PYL', 'SEC',  # Modified amino acids
    'HYP', 'PTR', 'MLY', 'M3L', 'CSO', 'OCS', 'CME',  # More modifications
}

# Standard nucleic acid residue names
NUCLEIC_ACIDS = {
    'A', 'C', 'G', 'U', 'T',          # Standard bases
    'DA', 'DC', 'DG', 'DT', 'DU',     # DNA
    '+A', '+C', '+G', '+U', '+T',     # RNA with modifications
    'ADE', 'CYT', 'GUA', 'URA', 'THY',  # Full names
}

# Common solvent residue names
SOLVENTS = {
    'HOH', 'WAT', 'H2O', 'DOD', 'D2O',  # Water
    'SO4', 'PO4', 'CL', 'NA', 'K', 'MG', 'CA', 'ZN',  # Common ions
    'GOL', 'EDO', 'PEG', 'MPD', 'DMS',  # Common cryo/crystallization agents
}

# Saccharide/sugar residue names
SACCHARIDES = {
    'GLC', 'GAL', 'MAN', 'FUC', 'XYL', 'RIB', 'NAG', 'BMA',
    'FUL', 'SIA', 'NDG', 'BGC', 'GLA', 'GCS', 'NGA', 'A2G',
    'RAM', 'ARA', 'LYX', 'ALL', 'ALT', 'GUP', 'IDO', 'TAL',
}

# Protein backbone atom names
PROTEIN_BACKBONE_ATOMS = {'N', 'CA', 'C', 'O', 'OXT'}

# Nucleic acid backbone atom names
NUCLEIC_BACKBONE_ATOMS = {"P", "OP1", "OP2", "O5'", "C5'", "C4'", "O4'", "C3'", "O3'", "C2'", "O2'", "C1'"}


class SelectionEvaluator:
    """Evaluates selection AST against gemmi Structure."""

    def __init__(self, structure):
        """
        Initialize evaluator with a gemmi Structure.

        Args:
            structure: gemmi.Structure object
        """
        self.structure = structure

    def evaluate(self, node: SelectionNode) -> List:
        """
        Evaluate selection AST and return list of selected atoms.

        Args:
            node: Root node of selection AST

        Returns:
            List of (model, chain, residue, atom) tuples for selected atoms

        """
        if isinstance(node, CIDSelector):
            return self._evaluate_cid(node)
        elif isinstance(node, CategorySelector):
            return self._evaluate_category(node)
        elif isinstance(node, LogicalAnd):
            left_atoms = set(self.evaluate(node.left))
            right_atoms = set(self.evaluate(node.right))
            return list(left_atoms & right_atoms)
        elif isinstance(node, LogicalOr):
            left_atoms = set(self.evaluate(node.left))
            right_atoms = set(self.evaluate(node.right))
            return list(left_atoms | right_atoms)
        elif isinstance(node, LogicalNot):
            # Get all atoms
            all_atoms = set(self._get_all_atoms())
            # Subtract selected atoms
            selected_atoms = set(self.evaluate(node.operand))
            return list(all_atoms - selected_atoms)
        else:
            raise ValueError(f"Unknown node type: {type(node)}")

    def _get_all_atoms(self) -> List[Tuple]:
        """Get all atoms in the structure."""
        atoms = []
        for model in self.structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atoms.append((model, chain, residue, atom))
        return atoms

    def _evaluate_cid(self, cid: CIDSelector) -> List[Tuple]:
        """
        Evaluate a CID selector.

        Returns list of (model, chain, residue, atom) tuples.
        """
        selected_atoms = []

        for model in self.structure:
            if not self._matches_model(model, cid):
                continue

            for chain in model:
                if not self._matches_chain(chain, cid):
                    continue

                for residue in chain:
                    if not self._matches_residue(residue, cid):
                        continue

                    for atom in residue:
                        if self._matches_atom(atom, cid):
                            selected_atoms.append((model, chain, residue, atom))

        return selected_atoms

    def _evaluate_category(self, cat: CategorySelector) -> List[Tuple]:
        """
        Evaluate a category selector.

        Returns list of (model, chain, residue, atom) tuples matching the category.
        """
        selected_atoms = []
        category = cat.category

        for model in self.structure:
            for chain in model:
                for residue in chain:
                    res_name = residue.name

                    # Determine residue type
                    is_protein = res_name in AMINO_ACIDS
                    is_nucleic = res_name in NUCLEIC_ACIDS
                    is_solvent = res_name in SOLVENTS
                    is_sugar = res_name in SACCHARIDES
                    is_polymer = is_protein or is_nucleic or is_sugar
                    # Ligand = not polymer, not solvent (i.e., small molecule)
                    is_ligand = not is_polymer and not is_solvent
                    # Hetero = has het flag or is not standard polymer
                    is_hetero = getattr(residue, 'het_flag', '') == 'H' or (not is_protein and not is_nucleic)

                    # Check if residue matches category (for residue-level categories)
                    residue_matches = False
                    if category == CategoryType.PROTEIN:
                        residue_matches = is_protein
                    elif category == CategoryType.NUCLEIC:
                        residue_matches = is_nucleic
                    elif category == CategoryType.SOLVENT:
                        residue_matches = is_solvent
                    elif category == CategoryType.SUGAR:
                        residue_matches = is_sugar
                    elif category == CategoryType.POLYMER:
                        residue_matches = is_polymer
                    elif category == CategoryType.LIGAND:
                        residue_matches = is_ligand
                    elif category == CategoryType.HETERO:
                        residue_matches = is_hetero
                    elif category in (CategoryType.BACKBONE, CategoryType.SIDECHAIN):
                        # For backbone/sidechain, we match residue but filter atoms
                        residue_matches = is_protein or is_nucleic

                    if not residue_matches:
                        continue

                    # Now filter atoms based on category
                    for atom in residue:
                        atom_matches = True

                        if category == CategoryType.BACKBONE:
                            # Backbone atoms only
                            if is_protein:
                                atom_matches = atom.name in PROTEIN_BACKBONE_ATOMS
                            elif is_nucleic:
                                atom_matches = atom.name in NUCLEIC_BACKBONE_ATOMS
                            else:
                                atom_matches = False
                        elif category == CategoryType.SIDECHAIN:
                            # Sidechain atoms only (not backbone)
                            if is_protein:
                                atom_matches = atom.name not in PROTEIN_BACKBONE_ATOMS
                            elif is_nucleic:
                                atom_matches = atom.name not in NUCLEIC_BACKBONE_ATOMS
                            else:
                                atom_matches = False

                        if atom_matches:
                            selected_atoms.append((model, chain, residue, atom))

        return selected_atoms

    def _matches_model(self, model, cid: CIDSelector) -> bool:
        """Check if model matches CID criteria."""
        if cid.model is None:
            return True

        model_name = str(model.name) if hasattr(model, 'name') and model.name else str(model.get_subchain(0).subchain_id()[0]) if model else "1"

        # Try to get model number - gemmi models don't have a direct number attribute
        # but we can use the model's position in the structure
        model_num = None
        for idx, m in enumerate(self.structure, start=1):
            if m == model:
                model_num = str(idx)
                break

        if model_num is None:
            model_num = "1"

        return self._matches_value(model_num, cid.model)

    def _matches_chain(self, chain, cid: CIDSelector) -> bool:
        """Check if chain matches CID criteria."""
        if cid.chain is None:
            return True

        chain_id = chain.name
        return self._matches_value(chain_id, cid.chain)

    def _matches_residue(self, residue, cid: CIDSelector) -> bool:
        """Check if residue matches CID criteria."""
        # Check sequence number
        if cid.seq_no is not None:
            seq_num = str(residue.seqid.num)
            if not self._matches_value(seq_num, cid.seq_no):
                return False

        # Check residue name
        if cid.res_name is not None:
            res_name = residue.name
            if not self._matches_value(res_name, cid.res_name):
                return False

        # Check insertion code
        if cid.ins_code is not None:
            ins_code = residue.seqid.icode if residue.seqid.icode else ''
            if not self._matches_value(ins_code, cid.ins_code):
                return False

        return True

    def _matches_atom(self, atom, cid: CIDSelector) -> bool:
        """Check if atom matches CID criteria."""
        # Check atom name
        if cid.atom_name is not None:
            atom_name = atom.name
            if not self._matches_value(atom_name, cid.atom_name):
                return False

        # Check chemical element
        if cid.chem_elem is not None:
            elem = atom.element.name
            if not self._matches_value(elem, cid.chem_elem):
                return False

        # Check alternate location
        if cid.alt_loc is not None:
            alt_loc = atom.altloc if atom.altloc else ''
            if not self._matches_value(alt_loc, cid.alt_loc):
                return False

        return True

    def _matches_value(self, actual: str, pattern: str) -> bool:
        """
        Check if actual value matches the pattern.

        Pattern can be:
        - Exact match: "A"
        - Wildcard: "*"
        - List: "A,B,C"
        - Range: "10-20" (for numbers only)
        """
        if pattern == '*':
            return True

        # Check for comma-separated list
        if ',' in pattern:
            values = [v.strip() for v in pattern.split(',')]
            return actual in values

        # Check for range (only for numeric values)
        if '-' in pattern and pattern[0].isdigit():
            parts = pattern.split('-')
            if len(parts) == 2:
                try:
                    start = int(parts[0])
                    end = int(parts[1])
                    actual_num = int(actual)
                    return start <= actual_num <= end
                except ValueError:
                    pass

        # Exact match
        return actual == pattern


def evaluate_selection(ast: SelectionNode, structure) -> List:
    """
    Evaluate a selection AST against a gemmi Structure.

    Args:
        ast: Root node of selection AST
        structure: gemmi.Structure object

    Returns:
        List of (model, chain, residue, atom) tuples for selected atoms

    Example:
        >>> import gemmi
        >>> structure = gemmi.read_structure('model.pdb')
        >>> from coordinate_selection import parse_selection, evaluate_selection
        >>> ast = parse_selection("A/27.A")
        >>> selected = evaluate_selection(ast, structure)
        >>> print(f"Selected {len(selected)} atoms")
    """
    evaluator = SelectionEvaluator(structure)
    return evaluator.evaluate(ast)
