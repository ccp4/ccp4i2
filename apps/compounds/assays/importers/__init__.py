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
"""
ADME Data Importers

This package provides parsers for importing ADME (Absorption, Distribution,
Metabolism, Excretion) assay data from external vendors into the CCP4i2
Assays data model.

Supported Vendors:
- NCU: Outsourced CRO for ADME studies

Usage:
    from apps.compounds.assays.importers import detect_parser, parse_adme_file

    # Auto-detect and parse
    parser = detect_parser(filepath)
    results = parser.parse(filepath)

    # Or use specific parser
    from apps.compounds.assays.importers.ncu import LiverMicrosomeParser
    parser = LiverMicrosomeParser()
    results = parser.parse(filepath)
"""

from .base import ADMEParser, ParsedResult, ParseResult, ValidationError
from .registry import detect_parser, get_parser, get_parser_by_slug, list_parsers, parse_adme_file

__all__ = [
    'ADMEParser',
    'ParsedResult',
    'ParseResult',
    'ValidationError',
    'detect_parser',
    'get_parser',
    'get_parser_by_slug',
    'list_parsers',
    'parse_adme_file',
]
