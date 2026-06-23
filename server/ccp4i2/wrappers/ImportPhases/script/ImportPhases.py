"""Import phases (an HL or Phi/FOM mini-MTZ) from a full MTZ file."""
from ccp4i2.core import CCP4XtalData
from ccp4i2.wrappers.import_common import CImportMiniMtzBase


class ImportPhases(CImportMiniMtzBase):

    TASKNAME = 'ImportPhases'
    GROUP_TYPE = 'Phs'
    OUTPUT_PARAM = 'PHSOUT'
    OUTPUT_CLASS = CCP4XtalData.CPhsDataFile
    TYPE_LABEL = 'phase'
