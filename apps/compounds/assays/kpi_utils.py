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


# Canonical list of valid units (normalized forms)
VALID_UNITS = {
    # Molar concentrations
    'nM',       # nanomolar
    'uM',       # micromolar (normalized from µM, μM)
    'mM',       # millimolar
    'pM',       # picomolar
    'M',        # molar
    # Time units
    'min',      # minutes
    's',        # seconds
    'h',        # hours
    # Percentage
    '%',        # percent
    # Rate units (e.g., clearance)
    'uL/min/mg',    # microsomal clearance (normalized from µL/min/mg)
    'mL/min/kg',    # hepatic clearance
    # Permeability units
    '1e-6 cm/s',    # Caco-2 Papp
    'cm/s',         # permeability
}

# Common invalid unit patterns with suggested corrections
UNIT_SUGGESTIONS = {
    # Percentage variants
    'pc': '%',
    'percent': '%',
    'pct': '%',
    'percentage': '%',
    # Missing or empty
    '': None,
    # Molar without case
    'nm': 'nM',
    'um': 'uM',
    'mm': 'mM',
    'pm': 'pM',
    # Unicode variants (should be caught by normalize_unit but just in case)
    'µM': 'uM',
    'μM': 'uM',
    'µm': 'uM',
    'μm': 'uM',
}


def validate_unit(unit: str | None) -> tuple[bool, str | None, str | None]:
    """
    Validate a unit string and return validation result.

    Args:
        unit: The unit string to validate (can be None or empty)

    Returns:
        Tuple of (is_valid, normalized_unit, error_message)
        - is_valid: True if unit is valid or None/empty
        - normalized_unit: The normalized form of the unit, or None
        - error_message: Human-readable error if invalid, None otherwise

    Examples:
        validate_unit('nM') -> (True, 'nM', None)
        validate_unit('µM') -> (True, 'uM', None)
        validate_unit('pc') -> (False, None, "Invalid unit 'pc'. Did you mean '%' for percent?")
        validate_unit(None) -> (True, None, None)  # None is acceptable (unit optional)
        validate_unit('') -> (True, None, None)    # Empty is acceptable
    """
    # None or empty is acceptable (unit is optional)
    if not unit or unit.strip() == '':
        return (True, None, None)

    unit = unit.strip()

    # First try to normalize the unit
    normalized = normalize_unit(unit)

    # Check if normalized form is valid
    if normalized in VALID_UNITS:
        return (True, normalized, None)

    # Check if original is valid (in case normalization didn't help)
    if unit in VALID_UNITS:
        return (True, unit, None)

    # Check for common mistakes and provide suggestions
    unit_lower = unit.lower()
    if unit_lower in UNIT_SUGGESTIONS:
        suggestion = UNIT_SUGGESTIONS[unit_lower]
        if suggestion:
            return (False, None, f"Invalid unit '{unit}'. Did you mean '{suggestion}'?")
        else:
            return (False, None, f"Invalid unit '{unit}'. Please specify a valid unit.")

    # Unknown unit - provide helpful error with valid options
    common_units = ['nM', 'uM', 'mM', '%', 'min', 's', 'h']
    return (False, None, f"Invalid unit '{unit}'. Valid units include: {', '.join(common_units)}")
