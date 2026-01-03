"""
Coordinate selection system for mmdb CID syntax.

This module provides parsing and evaluation of mmdb-style coordinate
selection syntax, replacing the legacy mmdb/ccp4mg dependency with
pure Python + gemmi.

Supports both traditional mmdb CID syntax and semantic category keywords:
    - CID syntax: A/27.A, {A/ or B/}, not (HOH), etc.
    - Categories: protein, nucleic, solvent, ligand, sugar, polymer,
                  backbone, sidechain, hetero

Public API:
    parse_selection(selection_string) -> SelectionAST
    evaluate_selection(ast, gemmi_structure) -> list[gemmi.Atom]
"""

from .parser import parse_selection
from .evaluator import evaluate_selection
from .ast_nodes import (
    CategoryType, CategorySelector, CIDSelector,
    LogicalAnd, LogicalOr, LogicalNot, SelectionNode
)

__all__ = [
    'parse_selection',
    'evaluate_selection',
    'CategoryType',
    'CategorySelector',
    'CIDSelector',
    'LogicalAnd',
    'LogicalOr',
    'LogicalNot',
    'SelectionNode',
]
