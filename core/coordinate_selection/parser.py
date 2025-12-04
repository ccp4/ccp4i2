"""
Parser for mmdb CID selection syntax.

Converts tokens into an AST (Abstract Syntax Tree).

Grammar (simplified):
    expression := or_expr
    or_expr := and_expr ('or' and_expr)*
    and_expr := not_expr ('and' not_expr)*
    not_expr := 'not' not_expr | primary
    primary := cid_selector | '(' expression ')' | '{' expression '}'
    cid_selector := ['/' model] ['/' chain] ['/' residue] ['/' atom]
"""

from typing import List, Optional
from .tokenizer import Token, TokenType, tokenize
from .ast_nodes import (
    CIDSelector, LogicalAnd, LogicalOr, LogicalNot, SelectionNode
)


class ParseError(Exception):
    """Raised when parsing fails."""
    pass


class Parser:
    """Recursive descent parser for selection syntax."""

    def __init__(self, tokens: List[Token]):
        self.tokens = tokens
        self.position = 0

    def current_token(self) -> Token:
        """Get the current token."""
        if self.position < len(self.tokens):
            return self.tokens[self.position]
        return self.tokens[-1]  # EOF

    def peek_token(self, offset: int = 1) -> Token:
        """Peek ahead at a token."""
        pos = self.position + offset
        if pos < len(self.tokens):
            return self.tokens[pos]
        return self.tokens[-1]  # EOF

    def advance(self) -> Token:
        """Consume and return the current token."""
        token = self.current_token()
        if token.type != TokenType.EOF:
            self.position += 1
        return token

    def expect(self, token_type: TokenType) -> Token:
        """Expect a specific token type, raise error if not found."""
        token = self.current_token()
        if token.type != token_type:
            raise ParseError(
                f"Expected {token_type.name} but got {token.type.name} "
                f"at position {token.position}"
            )
        return self.advance()

    def parse(self) -> SelectionNode:
        """Parse the entire expression."""
        result = self.parse_or_expr()
        if self.current_token().type != TokenType.EOF:
            raise ParseError(
                f"Unexpected token {self.current_token().type.name} "
                f"at position {self.current_token().position}"
            )
        return result

    def parse_or_expr(self) -> SelectionNode:
        """Parse OR expression (lowest precedence)."""
        left = self.parse_and_expr()

        while self.current_token().type == TokenType.OR:
            self.advance()  # consume 'or'
            right = self.parse_and_expr()
            left = LogicalOr(left, right)

        return left

    def parse_and_expr(self) -> SelectionNode:
        """Parse AND expression (medium precedence)."""
        left = self.parse_not_expr()

        while self.current_token().type == TokenType.AND:
            self.advance()  # consume 'and'
            right = self.parse_not_expr()
            left = LogicalAnd(left, right)

        return left

    def parse_not_expr(self) -> SelectionNode:
        """Parse NOT expression (highest precedence)."""
        if self.current_token().type == TokenType.NOT:
            self.advance()  # consume 'not'
            operand = self.parse_not_expr()  # Right-associative
            return LogicalNot(operand)

        return self.parse_primary()

    def parse_primary(self) -> SelectionNode:
        """Parse primary expression (CID selector or grouped expression)."""
        token = self.current_token()

        # Grouped expression with parentheses OR residue name selector
        if token.type == TokenType.LPAREN:
            # Look ahead to distinguish (resName) from (expr)
            # If it's (identifier) or (identifier,identifier,...) with no operators/slashes,
            # treat it as a residue name selector
            if self._is_residue_name_selector():
                # It's a (resName) selector - parse as CID with just res_name
                selector = CIDSelector()
                self.advance()  # consume '('
                selector.res_name = self._parse_value()
                self.expect(TokenType.RPAREN)
                return selector
            else:
                # It's a grouped expression like (A or B)
                self.advance()  # consume '('
                expr = self.parse_or_expr()
                self.expect(TokenType.RPAREN)
                return expr

        # Grouped expression with braces (CCP4i2 extension)
        if token.type == TokenType.LBRACE:
            self.advance()  # consume '{'
            expr = self.parse_or_expr()
            self.expect(TokenType.RBRACE)
            return expr

        # Otherwise, it must be a CID selector
        return self.parse_cid_selector()

    def _is_residue_name_selector(self) -> bool:
        """
        Check if current position is (resName) selector vs (expression).

        Returns True if pattern matches: ( identifier[,identifier]* )
        with no operators (and/or/not) or slashes.
        """
        if self.current_token().type != TokenType.LPAREN:
            return False

        # Look ahead past the (
        pos = self.position + 1
        if pos >= len(self.tokens):
            return False

        # Must start with identifier or number
        if self.tokens[pos].type not in (TokenType.IDENTIFIER, TokenType.NUMBER):
            return False

        # Scan until we find ) or something that indicates it's an expression
        pos += 1
        while pos < len(self.tokens):
            token = self.tokens[pos]

            if token.type == TokenType.RPAREN:
                # Found closing paren - it's a residue name
                return True
            elif token.type in (TokenType.SLASH, TokenType.AND, TokenType.OR,
                              TokenType.NOT, TokenType.LPAREN, TokenType.LBRACE):
                # These indicate it's an expression, not a simple residue name
                return False
            elif token.type in (TokenType.IDENTIFIER, TokenType.NUMBER,
                              TokenType.COMMA, TokenType.DASH, TokenType.WILDCARD):
                # These are OK in a residue name list
                pos += 1
            else:
                # Anything else ends the check
                return False

        return False

    def parse_cid_selector(self) -> CIDSelector:
        """
        Parse a CID selector.

        Format: /model/chain/seqNo(resName).insCode/atomName[chemElem]:altLoc

        Examples:
            A              -> chain='A'
            /1             -> model='1'
            /1/A           -> model='1', chain='A'
            A/27           -> chain='A', seq_no='27'
            A/27.A         -> chain='A', seq_no='27', ins_code='A'
            A/(ALA)        -> chain='A', res_name='ALA'
            A/CA[C]        -> chain='A', atom_name='CA', chem_elem='C'
            CA[CA]         -> atom_name='CA', chem_elem='CA'
            /*/A/10-20     -> chain='A', seq_no='10-20'
        """
        selector = CIDSelector()

        # Check if it starts with / (explicit model specification)
        starts_with_slash = self.current_token().type == TokenType.SLASH

        if starts_with_slash:
            self.advance()  # consume initial '/'
            # Next token is model number (or * or identifier)
            selector.model = self._parse_value()

            # Now we expect more fields separated by /
            if self.current_token().type == TokenType.SLASH:
                self.advance()
                selector.chain = self._parse_value()

                if self.current_token().type == TokenType.SLASH:
                    self.advance()
                    self._parse_residue(selector)

                    if self.current_token().type == TokenType.SLASH:
                        self.advance()
                        self._parse_atom(selector)
        else:
            # Starts without /, so first value could be chain or atom name
            first_value = self._parse_value()

            # Check if followed by / (indicates it was a chain)
            if self.current_token().type == TokenType.SLASH:
                selector.chain = first_value
                self.advance()
                self._parse_residue(selector)

                if self.current_token().type == TokenType.SLASH:
                    self.advance()
                    self._parse_atom(selector)
            else:
                # No /, so it might be an atom name or chain
                # Check for [element] which indicates atom
                if self.current_token().type == TokenType.LBRACKET:
                    selector.atom_name = first_value
                    self._parse_element(selector)
                    if self.current_token().type == TokenType.COLON:
                        self._parse_altloc(selector)
                else:
                    # Just a bare identifier - treat as chain
                    selector.chain = first_value

        return selector

    def _parse_value(self) -> str:
        """Parse a simple value (identifier, number, wildcard, or list)."""
        values = []

        while True:
            token = self.current_token()

            if token.type == TokenType.IDENTIFIER:
                values.append(token.value)
                self.advance()
            elif token.type == TokenType.NUMBER:
                # Check for range (e.g., 10-20)
                num = token.value
                self.advance()
                if self.current_token().type == TokenType.DASH:
                    self.advance()
                    if self.current_token().type == TokenType.NUMBER:
                        num += '-' + self.current_token().value
                        self.advance()
                values.append(num)
            elif token.type == TokenType.STAR:
                values.append('*')
                self.advance()
            else:
                break

            # Check for comma (list of values)
            if self.current_token().type == TokenType.COMMA:
                self.advance()
            else:
                break

        if not values:
            raise ParseError(f"Expected value at position {self.current_token().position}")

        return ','.join(values)

    def _parse_residue(self, selector: CIDSelector):
        """Parse residue specification: seqNo(resName).insCode"""
        # Sequence number or residue name in parentheses
        if self.current_token().type == TokenType.LPAREN:
            # (resName)
            self.advance()
            selector.res_name = self._parse_value()
            self.expect(TokenType.RPAREN)
        else:
            # seqNo
            selector.seq_no = self._parse_value()

            # Optional (resName)
            if self.current_token().type == TokenType.LPAREN:
                self.advance()
                selector.res_name = self._parse_value()
                self.expect(TokenType.RPAREN)

            # Optional .insCode
            if self.current_token().type == TokenType.DOT:
                self.advance()
                selector.ins_code = self._parse_value()

    def _parse_atom(self, selector: CIDSelector):
        """Parse atom specification: atomName[chemElem]:altLoc"""
        selector.atom_name = self._parse_value()

        # Optional [chemElem]
        if self.current_token().type == TokenType.LBRACKET:
            self._parse_element(selector)

        # Optional :altLoc
        if self.current_token().type == TokenType.COLON:
            self._parse_altloc(selector)

    def _parse_element(self, selector: CIDSelector):
        """Parse element specification: [chemElem]"""
        self.expect(TokenType.LBRACKET)
        selector.chem_elem = self._parse_value()
        self.expect(TokenType.RBRACKET)

    def _parse_altloc(self, selector: CIDSelector):
        """Parse alternate location: :altLoc"""
        self.expect(TokenType.COLON)
        selector.alt_loc = self._parse_value()


def parse_selection(selection_string: str) -> SelectionNode:
    """
    Parse a selection string into an AST.

    Args:
        selection_string: The selection string to parse

    Returns:
        Root node of the AST

    Raises:
        ParseError: If parsing fails

    Examples:
        >>> ast = parse_selection("A/27.A")
        >>> ast = parse_selection("{A/ or B/} and {(ALA)}")
        >>> ast = parse_selection("not (HOH)")
    """
    tokens = tokenize(selection_string)
    parser = Parser(tokens)
    return parser.parse()
