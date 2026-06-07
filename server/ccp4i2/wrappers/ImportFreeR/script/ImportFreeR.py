"""Import a free-R flag set (a mini-MTZ) from a full MTZ file."""
from ccp4i2.core import CCP4XtalData
from ccp4i2.wrappers.import_common import CImportMiniMtzBase


class ImportFreeR(CImportMiniMtzBase):

    TASKNAME = 'ImportFreeR'
    GROUP_TYPE = 'FreeR'
    OUTPUT_PARAM = 'FREEROUT'
    OUTPUT_CLASS = CCP4XtalData.CFreeRDataFile
    TYPE_LABEL = 'free-R flag'
