"""
Evaluator for coordinate selection AST using gemmi.

Applies selection criteria to gemmi.Structure objects.
"""

from typing import List, Set, Tuple
import re
from .ast_nodes import (
    CIDSelector, LogicalAnd, LogicalOr, LogicalNot, SelectionNode
)


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
