"""Import map coefficients (an FPHI mini-MTZ) from a full MTZ file."""
from ccp4i2.core import CCP4XtalData
from ccp4i2.wrappers.import_common import (
    CImportMiniMtzBase,
    infer_map_subtype,
)


class ImportMapCoeffs(CImportMiniMtzBase):

    TASKNAME = 'ImportMapCoeffs'
    GROUP_TYPE = 'MapCoeffs'
    OUTPUT_PARAM = 'FPHIOUT'
    OUTPUT_CLASS = CCP4XtalData.CMapCoeffsDataFile
    TYPE_LABEL = 'map coefficient'

    def finish_output(self, out_obj, group, labels):
        # Infer NORMAL / DIFFERENCE / ANOM-DIFFERENCE from the column names
        try:
            out_obj.subType.set(infer_map_subtype(labels))
        except Exception:
            pass
