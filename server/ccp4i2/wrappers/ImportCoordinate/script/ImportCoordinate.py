"""Import a coordinate file (PDB or mmCIF) into the project."""
from ccp4i2.wrappers.import_common import CImportFileBase


class ImportCoordinate(CImportFileBase):

    TASKNAME = 'ImportCoordinate'
    INPUT_PARAM = 'XYZIN'
    OUTPUT_PARAM = 'XYZOUT'
    ANNOTATION_PREFIX = 'Imported coordinates'

    def validate_source(self, src_path):
        try:
            import gemmi
            structure = gemmi.read_structure(src_path)
            if structure is None or len(structure) == 0:
                return 'Coordinate file contains no models.'
        except Exception as e:
            return 'Not a readable coordinate file: ' + str(e)
        return None
