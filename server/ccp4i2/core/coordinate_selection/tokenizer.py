"""
Tokenizer for mmdb CID selection syntax.

Converts selection strings into tokens for parsing.
"""

import re
from enum import Enum, auto
from typing import List, Tuple
from dataclasses import dataclass


class TokenType(Enum):
    """Token types for selection syntax."""
    # Literals
    SLASH = auto()           # /
    LPAREN = auto()          # (
    RPAREN = auto()          # )
    LBRACKET = auto()        # [
    RBRACKET = auto()        # ]
    LBRACE = auto()          # {
    RBRACE = auto()          # }
    COLON = auto()           # :
    DOT = auto()             # .
    COMMA = auto()           # ,
    DASH = auto()            # - (for ranges)
    STAR = auto()            # * (wildcard)

    # Logical keywords
    AND = auto()             # and
    OR = auto()              # or
    NOT = auto()             # not

    # Category keywords (semantic selectors)
    PROTEIN = auto()         # protein - amino acid residues
    NUCLEIC = auto()         # nucleic - nucleic acid residues
    SOLVENT = auto()         # solvent/water - HOH, WAT, etc.
    LIGAND = auto()          # ligand - non-polymer, non-solvent
    SUGAR = auto()           # sugar/saccharide - carbohydrate residues
    POLYMER = auto()         # polymer - protein + nucleic + sugar
    BACKBONE = auto()        # backbone - backbone atoms only
    SIDECHAIN = auto()       # sidechain - sidechain atoms only
    HETERO = auto()          # hetero/hetatm - HETATM records

    # Values
    IDENTIFIER = auto()      # Chain IDs, residue names, atom names, etc.
    NUMBER = auto()          # Model numbers, sequence numbers

    # Special
    EOF = auto()             # End of input


@dataclass
class Token:
    """A token from the selection string."""
    type: TokenType
    value: str
    position: int

    def __repr__(self):
        return f"Token({self.type.name}, '{self.value}', pos={self.position})"


class Tokenizer:
    """Tokenizes mmdb CID selection strings."""

    # Keywords (case-insensitive)
    KEYWORDS = {
        # Logical operators
        'and': TokenType.AND,
        'or': TokenType.OR,
        'not': TokenType.NOT,
        # Category selectors
        'protein': TokenType.PROTEIN,
        'nucleic': TokenType.NUCLEIC,
        'solvent': TokenType.SOLVENT,
        'water': TokenType.SOLVENT,      # alias for solvent
        'ligand': TokenType.LIGAND,
        'sugar': TokenType.SUGAR,
        'saccharide': TokenType.SUGAR,   # alias for sugar
        'polymer': TokenType.POLYMER,
        'backbone': TokenType.BACKBONE,
        'sidechain': TokenType.SIDECHAIN,
        'hetero': TokenType.HETERO,
        'hetatm': TokenType.HETERO,      # alias for hetero
    }

    def __init__(self, input_string: str):
        self.input = input_string
        self.position = 0
        self.tokens: List[Token] = []

    def tokenize(self) -> List[Token]:
        """
        Tokenize the input string.

        Returns:
            List of tokens

        Raises:
            ValueError: If tokenization fails
        """
        while self.position < len(self.input):
            # Skip whitespace (but track it for significance)
            if self.input[self.position].isspace():
                self.position += 1
                continue

            # Try to match a token
            if not self._match_next_token():
                raise ValueError(
                    f"Unexpected character '{self.input[self.position]}' "
                    f"at position {self.position}"
                )

        # Add EOF token
        self.tokens.append(Token(TokenType.EOF, '', self.position))
        return self.tokens

    def _match_next_token(self) -> bool:
        """
        Try to match and consume the next token.

        Returns:
            True if a token was matched, False otherwise
        """
        char = self.input[self.position]
        start_pos = self.position

        # Single-character tokens
        single_char_tokens = {
            '/': TokenType.SLASH,
            '(': TokenType.LPAREN,
            ')': TokenType.RPAREN,
            '[': TokenType.LBRACKET,
            ']': TokenType.RBRACKET,
            '{': TokenType.LBRACE,
            '}': TokenType.RBRACE,
            ':': TokenType.COLON,
            '.': TokenType.DOT,
            ',': TokenType.COMMA,
            '-': TokenType.DASH,
            '*': TokenType.STAR,
        }

        if char in single_char_tokens:
            self.tokens.append(Token(single_char_tokens[char], char, start_pos))
            self.position += 1
            return True

        # Numbers
        if char.isdigit():
            return self._match_number()

        # Identifiers (including keywords)
        if char.isalpha() or char == '_':
            return self._match_identifier()

        return False

    def _match_number(self) -> bool:
        """Match a number token."""
        start_pos = self.position
        value = ''

        while self.position < len(self.input) and self.input[self.position].isdigit():
            value += self.input[self.position]
            self.position += 1

        self.tokens.append(Token(TokenType.NUMBER, value, start_pos))
        return True

    def _match_identifier(self) -> bool:
        """Match an identifier (or keyword) token."""
        start_pos = self.position
        value = ''

        # Identifiers can contain letters, digits, underscores, and some special chars
        while self.position < len(self.input):
            char = self.input[self.position]
            if char.isalnum() or char in ['_', "'", '"']:
                value += char
                self.position += 1
            else:
                break

        # Check if it's a keyword (case-insensitive)
        lower_value = value.lower()
        if lower_value in self.KEYWORDS:
            token_type = self.KEYWORDS[lower_value]
        else:
            token_type = TokenType.IDENTIFIER

        self.tokens.append(Token(token_type, value, start_pos))
        return True


def tokenize(selection_string: str) -> List[Token]:
    """
    Tokenize a selection string.

    Args:
        selection_string: The selection string to tokenize

    Returns:
        List of tokens

    Example:
        >>> tokens = tokenize("A/27.A")
        >>> [t.type for t in tokens]
        [TokenType.IDENTIFIER, TokenType.SLASH, TokenType.NUMBER, ...]
    """
    tokenizer = Tokenizer(selection_string)
    return tokenizer.tokenize()
