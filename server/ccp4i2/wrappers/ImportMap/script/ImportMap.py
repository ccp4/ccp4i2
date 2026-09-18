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

    def finalize_output(self, out):
        """Tag the imported map with the kind the user chose (CMapDataFile
        SUBTYPE_*): normal / difference / anomalous / mask / half map. Half maps
        are one of a pair for cross-validated refinement (servalcat --halfmaps)."""
        control = self.container.controlParameters
        if not (hasattr(control, 'MAP_SUBTYPE') and control.MAP_SUBTYPE.isSet()):
            return
        try:
            out.subType.set(int(control.MAP_SUBTYPE))
        except Exception:
            out.subType = int(control.MAP_SUBTYPE)
