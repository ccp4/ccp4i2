# Copyright (C) 2025 University of York
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
    from ccp4i2.core.conversions import ObsDataConverter
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
