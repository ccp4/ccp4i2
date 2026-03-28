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
NCU ADME Data Parsers.

This package provides parsers for ADME data files from NCU (outsourced CRO).

Supported Assay Types:
- LM: Liver Microsome Stability
- BS: Blood/Serum Stability
- GSH: GSH Stability (Glutathione)
- Caco-2: Caco-2 Permeability
"""

from .lm import LiverMicrosomeParser
from .bs import BloodSerumStabilityParser
from .gsh import GSHStabilityParser
from .caco2 import Caco2PermeabilityParser

__all__ = [
    'LiverMicrosomeParser',
    'BloodSerumStabilityParser',
    'GSHStabilityParser',
    'Caco2PermeabilityParser',
]
