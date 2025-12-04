"""
Crystallographic data format conversion modules.

This package provides converters for various crystallographic data transformations:

1. **ObsDataConverter** (obs_data_converter.py)
   MTZ observation data conversions:
   - IPAIR (I+/I-) ↔ FPAIR (F+/F-) ↔ IMEAN (I) ↔ FMEAN (F)
   - Fully implemented with ctruncate and gemmi

2. **PhaseDataConverter** (phase_data_converter.py)
   Phase data conversions:
   - HL (Hendrickson-Lattman) ↔ PHIFOM (PHI/FOM) ↔ FPHI (F/PHI)
   - TODO: Requires implementation

3. **ModelConverter** (model_converter.py)
   Macromolecular model format conversions:
   - PDB ↔ mmCIF
   - TODO: Requires gemmi-based implementation

Usage:
    from core.conversions import ObsDataConverter
    output_path = ObsDataConverter.to_fmean(obs_file, work_directory="./work")
"""

from .obs_data_converter import ObsDataConverter
from .phase_data_converter import PhaseDataConverter
from .model_converter import ModelConverter

__all__ = [
    'ObsDataConverter',
    'PhaseDataConverter',
    'ModelConverter',
]
