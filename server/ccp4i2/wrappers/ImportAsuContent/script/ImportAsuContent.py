"""Import an existing ASU contents definition (.asu.xml) into the project.

This is a straight pass-through for an ASU contents file produced elsewhere.
To *build* ASU contents from sequences, use the ProvideAsuContents task instead.
"""
from ccp4i2.wrappers.import_common import CImportFileBase


class ImportAsuContent(CImportFileBase):

    TASKNAME = 'ImportAsuContent'
    INPUT_PARAM = 'ASUIN'
    OUTPUT_PARAM = 'ASUOUT'
    ANNOTATION_PREFIX = 'Imported ASU contents'

    def validate_source(self, src_path):
        try:
            from lxml import etree
            etree.parse(src_path)
        except Exception as e:
            return 'Not a readable ASU contents (XML) file: ' + str(e)
        return None
