"""
NCU ADME Data Parsers.

This package provides parsers for ADME data files from NCU (outsourced CRO).

Supported Assay Types:
- LM: Liver Microsome Stability
- BS: Blood/Serum Stability
- GSH: GSH Stability (Glutathione)
- Caco-2: Caco-2 Permeability
- Hepatocytes: Hepatocyte Stability
- LogD: Distribution coefficient at pH 7.4
- PPB: Plasma Protein Binding
- Solubility: Kinetic solubility at pH 7.4
"""

from .lm import LiverMicrosomeParser
from .bs import BloodSerumStabilityParser
from .gsh import GSHStabilityParser
from .caco2 import Caco2PermeabilityParser
from .hepatocytes import HepatocyteStabilityParser
from .logd import LogDParser
from .ppb import PlasmaProteinBindingParser
from .solubility import SolubilityParser

__all__ = [
    'LiverMicrosomeParser',
    'BloodSerumStabilityParser',
    'GSHStabilityParser',
    'Caco2PermeabilityParser',
    'HepatocyteStabilityParser',
    'LogDParser',
    'PlasmaProteinBindingParser',
    'SolubilityParser',
]
