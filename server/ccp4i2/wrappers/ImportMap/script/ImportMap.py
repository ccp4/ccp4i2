"""Import a CCP4/MRC map file into the project."""
from ccp4i2.wrappers.import_common import CImportFileBase


class ImportMap(CImportFileBase):

    TASKNAME = 'ImportMap'
    INPUT_PARAM = 'MAPIN'
    OUTPUT_PARAM = 'MAPOUT'
    ANNOTATION_PREFIX = 'Imported map'

    def validate_source(self, src_path):
        try:
            import gemmi
            gemmi.read_ccp4_map(src_path)
        except Exception as e:
            return 'Not a readable CCP4/MRC map file: ' + str(e)
        return None
