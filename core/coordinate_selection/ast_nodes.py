"""
AST (Abstract Syntax Tree) node classes for coordinate selection.

Represents the parsed structure of mmdb CID selection expressions.
"""

from typing import Optional, Union, List
from dataclasses import dataclass


@dataclass
class CIDSelector:
    """
    Represents a basic mmdb CID selector.

    Format: /model/chain/seqNo(resName).insCode/atomName[chemElem]:altLoc

    Each field can be:
    - None (matches any)
    - A specific value
    - A wildcard '*' (matches any)
    - A range 'start-end' (for seqNo)
    - A comma-separated list 'A,B,C'
    """
    model: Optional[str] = None      # Model number (e.g., "1", "*", "1,2")
    chain: Optional[str] = None      # Chain ID (e.g., "A", "*", "A,B,C")
    seq_no: Optional[str] = None     # Sequence number (e.g., "27", "10-20", "*")
    res_name: Optional[str] = None   # Residue name (e.g., "ALA", "*")
    ins_code: Optional[str] = None   # Insertion code (e.g., "A", "*")
    atom_name: Optional[str] = None  # Atom name (e.g., "CA", "*")
    chem_elem: Optional[str] = None  # Chemical element (e.g., "C", "CA")
    alt_loc: Optional[str] = None    # Alternate location (e.g., "A", "*")

    def __repr__(self):
        parts = []
        if self.model is not None:
            parts.append(f"model={self.model}")
        if self.chain is not None:
            parts.append(f"chain={self.chain}")
        if self.seq_no is not None:
            parts.append(f"seq={self.seq_no}")
        if self.res_name is not None:
            parts.append(f"res={self.res_name}")
        if self.ins_code is not None:
            parts.append(f"ins={self.ins_code}")
        if self.atom_name is not None:
            parts.append(f"atom={self.atom_name}")
        if self.chem_elem is not None:
            parts.append(f"elem={self.chem_elem}")
        if self.alt_loc is not None:
            parts.append(f"alt={self.alt_loc}")
        return f"CIDSelector({', '.join(parts)})"


@dataclass
class LogicalAnd:
    """Logical AND of two selection expressions."""
    left: 'SelectionNode'
    right: 'SelectionNode'

    def __repr__(self):
        return f"({self.left} AND {self.right})"


@dataclass
class LogicalOr:
    """Logical OR of two selection expressions."""
    left: 'SelectionNode'
    right: 'SelectionNode'

    def __repr__(self):
        return f"({self.left} OR {self.right})"


@dataclass
class LogicalNot:
    """Logical NOT of a selection expression."""
    operand: 'SelectionNode'

    def __repr__(self):
        return f"(NOT {self.operand})"


# Type alias for any selection node
SelectionNode = Union[CIDSelector, LogicalAnd, LogicalOr, LogicalNot]
