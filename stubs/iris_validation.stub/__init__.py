"""
Minimal stub for iris_validation module.

This provides minimal stubs for the IRIS (Interactive Ramachandran Interactive System)
validation library used by validate_protein.

Most functionality will raise NotImplementedError if actually called.
This stub allows plugins to import iris_validation for registration purposes.
"""

__version__ = "stub-0.0.1"


def validate(*args, **kwargs):
    """Stub for iris validation function."""
    raise NotImplementedError("iris_validation.validate is not implemented in stub")


class IrisValidator:
    """Stub for IRIS validator class."""

    def __init__(self, *args, **kwargs):
        raise NotImplementedError("IrisValidator is not implemented in stub")
