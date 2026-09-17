"""Import an unmerged reflection file (MTZ, XDS HKL, scalepack, mmCIF) into the project."""
from ccp4i2.wrappers.import_common import CImportFileBase


def _looks_like_space_group(text):
    """Guard the annotation against a mis-parsed header.

    The unmerged loader reads the space group positionally for the text
    formats, so a file it has misidentified can yield a whole header line
    here -- an XDS_ASCII file taken for scalepack hands back
    ``MERGE=FALSE    FRIEDEL'S_LAW=FALSE``. A real Hermann-Mauguin symbol is
    short and made of letters, digits, spaces, slashes and minus signs, so
    anything else is dropped rather than shown to the user as fact.
    """
    text = text.strip()
    if not text or len(text) > 16:
        return False
    return all(char.isalnum() or char in ' /-' for char in text)


class ImportUnmerged(CImportFileBase):

    TASKNAME = 'ImportUnmerged'
    INPUT_PARAM = 'UNMERGEDIN'
    OUTPUT_PARAM = 'UNMERGEDOUT'
    ANNOTATION_PREFIX = 'Imported unmerged data'

    def validate_source(self, src_path):
        """Refuse anything the unmerged loader cannot read, and refuse a merged
        file offered as unmerged.

        The merged check is the one piece of judgement in the task, and it earns
        its place: Diamond's auto-processing publishes scaled.mtz beside
        scaled_unmerged.mtz, and picking the wrong one otherwise fails much
        later inside aimless with a message about batches.

        The verdict comes from CUnmergedDataContent.merged rather than a direct
        column scan, so it covers scalepack and mmCIF too and stays in step with
        the loader every other consumer of these files uses.
        """
        from ccp4i2.core.CCP4XtalData import CUnmergedDataFile

        unmerged = CUnmergedDataFile()
        unmerged.setFullPath(src_path)
        try:
            unmerged.loadFile()
            content = unmerged.getFileContent()
        except Exception as e:
            return 'Not a readable unmerged reflection file: ' + str(e)

        if content is None:
            return 'Could not read the contents of this reflection file.'

        file_format = str(content.format)
        if not file_format or file_format == 'unk':
            return ('Could not determine the format of this file. Expected an '
                    'unmerged MTZ, XDS HKL, scalepack sca, SAINT, SHELX or mmCIF file.')

        if str(content.merged) == 'merged':
            return ('This file contains merged reflections, not unmerged ones. '
                    'Use ImportObs for merged data.')

        return None

    def finalize_output(self, out):
        """Carry the file's own facts into the annotation, so the project file
        list says what the data is rather than only where it came from."""
        content = out.getFileContent()
        if content is None:
            return

        parts = [str(content.format).upper()]
        if content.spaceGroup.isSet() and _looks_like_space_group(str(content.spaceGroup)):
            parts.append(str(content.spaceGroup))
        if content.knowncell and content.cell.isSet():
            cell = content.cell
            parts.append('%.1f %.1f %.1f' % (float(cell.a), float(cell.b), float(cell.c)))
        if content.knownwavelength and content.wavelength.isSet():
            parts.append('%.4f A' % float(content.wavelength))

        out.annotation.set(str(out.annotation) + ' (' + ', '.join(parts) + ')')
