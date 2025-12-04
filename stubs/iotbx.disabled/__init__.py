"""
Minimal stub for iotbx (IOTBX - Input/Output Toolbox) module.

This provides minimal stubs for the iotbx crystallography I/O library from CCTBX.
The real iotbx module is distributed with CCTBX and contains compiled extensions.

Most functionality will raise NotImplementedError if actually called.
This stub allows plugins to import iotbx for registration purposes.
"""

__version__ = "stub-0.0.1"


class pdb:
    """Stub for iotbx.pdb module."""

    class hierarchy:
        """Stub for PDB hierarchy functionality."""

        @staticmethod
        def input(*args, **kwargs):
            raise NotImplementedError("iotbx.pdb.hierarchy.input is not implemented in stub")


class cif:
    """Stub for iotbx.cif module."""

    @staticmethod
    def reader(*args, **kwargs):
        raise NotImplementedError("iotbx.cif.reader is not implemented in stub")
