"""
KPI utility functions.

Provides helpers for parsing and normalizing KPI units from field names.
"""

import re

# Unit patterns to recognize in field names
# Matches patterns in parentheses or square brackets like "(nM)", "[µM]", "(min)", "(%)"
# Supports:
#   - Molar units: nM, µM/uM, mM, pM, M
#   - Time units: min, s, h
#   - Percentage: %
#   - Rate units: µL/min/mg, uL/min/mg
#   - Permeability: 1e-6 cm/s, 10^-6 cm/s
UNIT_PATTERN = re.compile(
    r'[\(\[]('
    # Molar concentrations
    r'[nμµu]M|'                    # nM, μM, µM, uM (nano/micro molar)
    r'mM|'                         # mM (millimolar)
    r'pM|'                         # pM (picomolar)
    r'M|'                          # M (molar)
    # Time units
    r'min|'                        # minutes
    r's|'                          # seconds
    r'h|'                          # hours
    # Percentage
    r'%|'                          # percentage
    # Rate units (e.g., clearance)
    r'[μµu]L/min/mg|'              # µL/min/mg (microsomal clearance)
    r'mL/min/kg|'                  # mL/min/kg (hepatic clearance)
    # Permeability units
    r'1e-6\s*cm/s|'                # 1e-6 cm/s (Caco-2 Papp)
    r'10\^?-6\s*cm/s|'             # 10^-6 cm/s or 10-6 cm/s
    r'cm/s'                        # cm/s
    r')[\)\]]',
    re.IGNORECASE
)

# Normalize unicode mu (μ, µ) to ASCII 'u' for consistent storage
UNIT_NORMALIZATION = {
    # Molar units
    'μm': 'uM',
    'µm': 'uM',
    'μM': 'uM',
    'µM': 'uM',
    'um': 'uM',
    'UM': 'uM',
    'nm': 'nM',
    'NM': 'nM',
    'mm': 'mM',
    'MM': 'mM',
    'pm': 'pM',
    'PM': 'pM',
    'm': 'M',
    # Time units
    'MIN': 'min',
    'Min': 'min',
    'S': 's',
    'H': 'h',
    # Rate units
    'μL/min/mg': 'uL/min/mg',
    'µL/min/mg': 'uL/min/mg',
    'μl/min/mg': 'uL/min/mg',
    'µl/min/mg': 'uL/min/mg',
    'ul/min/mg': 'uL/min/mg',
    'UL/MIN/MG': 'uL/min/mg',
    'ML/MIN/KG': 'mL/min/kg',
    'ml/min/kg': 'mL/min/kg',
    # Permeability
    '1e-6 cm/s': '1e-6 cm/s',
    '1e-6cm/s': '1e-6 cm/s',
    '10^-6 cm/s': '1e-6 cm/s',
    '10^-6cm/s': '1e-6 cm/s',
    '10-6 cm/s': '1e-6 cm/s',
    '10-6cm/s': '1e-6 cm/s',
    'CM/S': 'cm/s',
}


def normalize_unit(unit: str) -> str:
    """
    Normalize unit string to standard ASCII form.

    Converts unicode mu (μ, µ) to ASCII 'u' and ensures consistent casing.

    Args:
        unit: Raw unit string (e.g., 'µM', 'uM', 'nM')

    Returns:
        Normalized unit string (e.g., 'uM', 'nM', 'mM')
    """
    return UNIT_NORMALIZATION.get(unit, unit)


def parse_unit_from_field_name(field_name: str) -> str | None:
    """
    Extract unit from a KPI field name.

    Parses common unit notation patterns in parentheses or square brackets.
    Supports molar concentrations, time units, percentages, rate units, and
    permeability units.

    Examples:
        "PEO4 GI50 (µM)" -> "uM"
        "IC50 (nM)" -> "nM"
        "IC50 [nM]" -> "nM"          # square brackets
        "Hill coefficient" -> None
        "EC50 (uM)" -> "uM"
        "Inhibition (%)" -> "%"
        "Half-life (min)" -> "min"
        "CLint (µL/min/mg)" -> "uL/min/mg"
        "Papp (1e-6 cm/s)" -> "1e-6 cm/s"

    Args:
        field_name: The KPI field name to parse

    Returns:
        Normalized unit string or None if no unit found
    """
    if not field_name:
        return None

    match = UNIT_PATTERN.search(field_name)
    if match:
        return normalize_unit(match.group(1))
    return None
