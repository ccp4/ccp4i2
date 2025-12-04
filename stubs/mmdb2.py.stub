"""
Minimal stub for mmdb2 (MMDB2 - Macromolecular Database) module.

This provides minimal stubs for the mmdb2 functionality used by acedrgNew's atomMatching.
Most functionality will raise NotImplementedError if actually called.
"""

# Constants used by atomMatching.py
nAminoacidNames = 20  # Standard amino acids
nNucleotideNames = 8  # Standard nucleotides

# Selection type constants
STYPE_RESIDUE = 1
SKEY_NEW = 0

# Dummy objects for attribute access
class AAProperties:
    """Stub for amino acid properties."""
    pass

class NucleotideName:
    """Stub for nucleotide names."""
    pass


def InitMatType():
    """Initialize material type (stub - does nothing)."""
    pass


def getAAProperty(props, index):
    """Get amino acid property by index (stub)."""
    class AAProperty:
        name = f"AA{index}"
    return AAProperty()


def getPstr(names, index):
    """Get nucleotide name by index (stub)."""
    return f"NUC{index}"


class intp:
    """Integer pointer wrapper (stub)."""
    def __init__(self):
        self._value = 0

    def value(self):
        return self._value


class Manager:
    """MMDB2 Manager class for handling molecular structures (stub)."""

    def ReadCoorFile(self, filename):
        """Read coordinate file (stub - not implemented)."""
        raise NotImplementedError("mmdb2.Manager.ReadCoorFile is not implemented in stub")

    def NewSelection(self):
        """Create new selection (stub - not implemented)."""
        raise NotImplementedError("mmdb2.Manager.NewSelection is not implemented in stub")

    def Select(self, selHnd, stype, selection, skey):
        """Perform selection (stub - not implemented)."""
        raise NotImplementedError("mmdb2.Manager.Select is not implemented in stub")


def GetResidueSelIndex(molHnd, selHnd, selindexp):
    """Get residue selection index (stub - not implemented)."""
    raise NotImplementedError("mmdb2.GetResidueSelIndex is not implemented in stub")


def getPCResidue(selRes, index):
    """Get residue by index (stub - not implemented)."""
    raise NotImplementedError("mmdb2.getPCResidue is not implemented in stub")
