"""Import observed reflections (a mini-MTZ) from a full MTZ file."""
from ccp4i2.core import CCP4XtalData
from ccp4i2.wrappers.import_common import CImportMiniMtzBase


class ImportObs(CImportMiniMtzBase):

    TASKNAME = 'ImportObs'
    GROUP_TYPE = 'Obs'
    OUTPUT_PARAM = 'OBSOUT'
    OUTPUT_CLASS = CCP4XtalData.CObsDataFile
    TYPE_LABEL = 'observation'
