"""Import a monomer / restraint dictionary (CIF) into the project."""
from ccp4i2.wrappers.import_common import CImportFileBase


class ImportDictionary(CImportFileBase):

    TASKNAME = 'ImportDictionary'
    INPUT_PARAM = 'DICTIN'
    OUTPUT_PARAM = 'DICTOUT'
    ANNOTATION_PREFIX = 'Imported dictionary'

    def validate_source(self, src_path):
        try:
            import gemmi
            doc = gemmi.cif.read(src_path)
        except Exception as e:
            return 'Not a readable CIF dictionary: ' + str(e)
        # A restraint dictionary should carry at least one _chem_comp block
        has_chem_comp = any(
            block.find_loop('_chem_comp.id')
            or block.find_loop('_chem_comp_atom.atom_id')
            or block.find_pair('_chem_comp.id')
            for block in doc
        )
        if not has_chem_comp:
            return ('CIF file has no _chem_comp data - it does not look like a '
                    'monomer/restraint dictionary.')
        return None
