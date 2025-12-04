"""
Coordinate selection system for mmdb CID syntax.

This module provides parsing and evaluation of mmdb-style coordinate
selection syntax, replacing the legacy mmdb/ccp4mg dependency with
pure Python + gemmi.

Public API:
    parse_selection(selection_string) -> SelectionAST
    evaluate_selection(ast, gemmi_structure) -> list[gemmi.Atom]
"""

from .parser import parse_selection
from .evaluator import evaluate_selection

__all__ = ['parse_selection', 'evaluate_selection']
