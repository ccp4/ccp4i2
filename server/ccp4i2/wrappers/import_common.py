"""
Shared helpers for the Import<Type> family of tasks.

These tasks give the CLI/API a clean, validated "load a file into a project"
operation for each key CCP4i2 data type. Two flavours of base class live here:

* ``CImportFileBase`` - a pure pass-through copy for file types that have no
  internal column structure (coordinates, dictionaries, maps, sequences, ASU
  contents). The input file is copied into the job directory as a typed output
  object; the gleaner registers it with the project automatically.

* ``CImportMiniMtzBase`` - for the mini-MTZ family (observations, map
  coefficients, free-R flags, phases). A full MTZ is supplied and a single
  typed mini-MTZ is extracted. Column handling is *hybrid*: if exactly one
  column group of the target type is present it is used silently; if it is
  ambiguous (or absent) the user must name the columns via the COLUMNS
  control parameter, and ``runTimeValidity`` blocks submission until they do.
"""
import shutil
from pathlib import Path

from ccp4i2.core.CCP4PluginScript import CPluginScript
from ccp4i2.core.CCP4Utils import split_mtz_file


# ---------------------------------------------------------------------------
# Column-group utilities (mini-MTZ imports)
# ---------------------------------------------------------------------------

def get_groups_of_type(mtz_full_path, group_type):
    """Return the list of CColumnGroup objects of ``group_type`` in an MTZ.

    ``group_type`` is one of 'Obs', 'Phs', 'MapCoeffs', 'FreeR'. Uses gemmi via
    CMtzData.getColumnGroups(), which pattern-matches MTZ column *types* against
    each mini-MTZ class's signature, so it works regardless of the user's
    column labels.
    """
    from ccp4i2.core.CCP4XtalData import CMtzDataFile

    if not mtz_full_path or not Path(mtz_full_path).exists():
        return []
    mtz = CMtzDataFile()
    mtz.setFullPath(str(mtz_full_path))
    mtz.loadFile()
    contents = mtz.getFileContent()
    if contents is None:
        return []
    groups = contents.getColumnGroups() or []
    return [g for g in groups if str(g.columnGroupType) == group_type]


def group_labels(group):
    """Return the source column labels of a column group as a list of str."""
    labels = []
    for col in group.columnList:
        label = col.columnLabel
        if hasattr(label, 'value'):
            label = label.value
        labels.append(str(label))
    return labels


def group_content_flag(group):
    """Return the integer contentFlag of a column group (or None)."""
    flag = group.contentFlag
    if hasattr(flag, 'value'):
        flag = flag.value
    try:
        return int(flag)
    except (TypeError, ValueError):
        return None


def group_dataset(group):
    """Return the dataset name of a column group ('' if unset)."""
    dataset = group.dataset
    if hasattr(dataset, 'value'):
        dataset = dataset.value
    return str(dataset) if dataset else ''


def resolve_group(candidates, requested_labels):
    """Pick the column group to import.

    Returns ``(group, status)`` where ``status`` is one of:
      * 'ok'        - ``group`` is the one to use
      * 'none'      - no column group of this type exists in the file
      * 'ambiguous' - more than one candidate and the user named no columns
      * 'no_match'  - the user named columns but none of the candidates match

    ``requested_labels`` is a list of source column labels (or None/empty).
    """
    if not candidates:
        return None, 'none'

    if requested_labels:
        wanted = [lbl.strip() for lbl in requested_labels if lbl.strip()]
        if wanted:
            for group in candidates:
                if group_labels(group) == wanted:
                    return group, 'ok'
            # Looser match: first requested label appears in the group
            for group in candidates:
                if wanted[0] in group_labels(group):
                    return group, 'ok'
            return None, 'no_match'

    if len(candidates) == 1:
        return candidates[0], 'ok'
    return None, 'ambiguous'


def describe_candidates(candidates):
    """Human-readable summary of candidate groups for error messages."""
    parts = []
    for group in candidates:
        dataset = group_dataset(group)
        labels = ','.join(group_labels(group))
        prefix = (dataset + '/') if dataset else ''
        parts.append(prefix + '[' + labels + ']')
    return '; '.join(parts)


# Map coefficient column-name heuristics for subType inference (refmac / phenix
# / buster conventions). 1=normal (2Fo-Fc), 2=difference (Fo-Fc), 3=anomalous.
_MAP_DIFFERENCE_PATTERNS = [
    'DELFWT', 'PHDELWT', 'DELPHWT', 'FOFC', 'FOFCWT', 'PHFOFCWT',
    'MFODFCWT', 'PHMFODFCWT', 'DIFF',
]
_MAP_ANOMALOUS_PATTERNS = [
    'FAN', 'PHAN', 'APTS', 'PHATS', 'SIGFAN', 'DANO', 'SIGDANO', 'PHIANO', 'ANOM',
]


def infer_map_subtype(column_labels):
    """Infer CMapCoeffsDataFile subType (1/2/3) from column labels."""
    from ccp4i2.core.CCP4XtalData import CMapCoeffsDataFile

    upper = [lbl.upper() for lbl in column_labels]
    for lbl in upper:
        for pattern in _MAP_ANOMALOUS_PATTERNS:
            if pattern in lbl:
                return CMapCoeffsDataFile.SUBTYPE_ANOM_DIFFERENCE
    for lbl in upper:
        for pattern in _MAP_DIFFERENCE_PATTERNS:
            if pattern in lbl:
                # 2Fo-Fc style maps carry FOFC substrings but are *normal* maps
                if '2FOFC' in lbl or '2FO-FC' in lbl:
                    return CMapCoeffsDataFile.SUBTYPE_NORMAL
                return CMapCoeffsDataFile.SUBTYPE_DIFFERENCE
    return CMapCoeffsDataFile.SUBTYPE_NORMAL


def canonical_mapping(output_class, content_flag, labels):
    """Map source labels to the canonical mini-MTZ labels for this content flag.

    Produces a clean mini-MTZ whose columns are named per the class's
    CONTENT_SIGNATURE_LIST (e.g. F,SIGF for mean SFs; FREER for free flags).
    Falls back to identity if the column count doesn't match.
    """
    signature_list = getattr(output_class, 'CONTENT_SIGNATURE_LIST', None)
    if content_flag and signature_list and 1 <= content_flag <= len(signature_list):
        signature = signature_list[content_flag - 1]
        if signature and len(signature) == len(labels):
            return {src: canon for src, canon in zip(labels, signature)}
    return {lbl: lbl for lbl in labels}


# ---------------------------------------------------------------------------
# Base: pure pass-through file import (no columns)
# ---------------------------------------------------------------------------

class CImportFileBase(CPluginScript):
    """Copy one input file into the job directory as a typed output object.

    Subclasses set ``INPUT_PARAM`` / ``OUTPUT_PARAM`` (names in inputData /
    outputData) and ``ANNOTATION_PREFIX``. ContentFlag is auto-introspected
    (PDB vs mmCIF etc.); the output is registered with the project by the
    gleaner with no further action needed.
    """

    INPUT_PARAM = None
    OUTPUT_PARAM = None
    ANNOTATION_PREFIX = 'Imported file'

    ERROR_CODES = {
        200: {'description': 'Input file not set'},
        201: {'description': 'Input file not found'},
        202: {'description': 'Failed to copy imported file'},
        203: {'description': 'Imported file failed validation'},
    }

    def _input_object(self):
        return getattr(self.container.inputData, self.INPUT_PARAM)

    def _output_object(self):
        return getattr(self.container.outputData, self.OUTPUT_PARAM)

    def validate_source(self, src_path):
        """Optional subclass hook: return an error string if invalid, else None."""
        return None

    def runTimeValidity(self):
        from ccp4i2.core import CCP4ErrorHandling

        error = super().runTimeValidity()
        if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
            return error
        inp = self._input_object()
        if not inp.isSet():
            return error
        src = str(inp.fullPath)
        message = self.validate_source(src)
        if message:
            error.append(
                self.TASKNAME, 203, message,
                self.TASKNAME + '.container.inputData.' + self.INPUT_PARAM,
                CCP4ErrorHandling.SEVERITY_ERROR,
            )
        return error

    def startProcess(self):
        inp = self._input_object()
        out = self._output_object()

        if not inp.isSet():
            self.appendErrorReport(200, self.INPUT_PARAM)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        src = str(inp.fullPath)
        if not src or not Path(src).exists():
            self.appendErrorReport(201, src)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        # Preserve the input extension so the output format is unambiguous
        dst = str(out.fullPath)
        src_ext = Path(src).suffix
        if src_ext and Path(dst).suffix.lower() != src_ext.lower():
            dst = str(Path(dst).with_suffix(src_ext))
            out.setFullPath(dst)

        try:
            shutil.copyfile(src, dst)
        except Exception as e:
            self.appendErrorReport(202, str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        # Auto-detect content flag where the data type supports it
        try:
            out.setContentFlag()
        except Exception:
            pass

        annotation = self.ANNOTATION_PREFIX + ' ' + Path(src).name
        try:
            out.annotation.set(annotation)
        except Exception:
            out.annotation = annotation

        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED


# ---------------------------------------------------------------------------
# Base: mini-MTZ import (column-bearing)
# ---------------------------------------------------------------------------

class CImportMiniMtzBase(CPluginScript):
    """Extract a single typed mini-MTZ from a full MTZ.

    Subclasses set ``GROUP_TYPE`` ('Obs'/'Phs'/'MapCoeffs'/'FreeR'),
    ``OUTPUT_PARAM`` (name in outputData), ``OUTPUT_CLASS`` (the CData class),
    and ``TYPE_LABEL`` (for messages). Input is ``HKLIN`` (CMtzDataFile);
    optional control parameter ``COLUMNS`` (comma-separated source labels)
    disambiguates when more than one column group of the target type exists.
    """

    GROUP_TYPE = None
    OUTPUT_PARAM = None
    OUTPUT_CLASS = None
    TYPE_LABEL = 'reflection'

    ERROR_CODES = {
        210: {'description': 'Input MTZ not found'},
        211: {'description': 'No column group of the requested type found'},
        212: {'description': 'Ambiguous column groups - specify columns'},
        213: {'description': 'Requested columns did not match any column group'},
        214: {'description': 'Failed to split MTZ file'},
    }

    def _requested_labels(self):
        control = self.container.controlParameters
        if hasattr(control, 'COLUMNS') and control.COLUMNS.isSet():
            text = str(control.COLUMNS)
            if text:
                return [lbl.strip() for lbl in text.split(',') if lbl.strip()]
        return None

    def _candidate_groups(self):
        hklin = self.container.inputData.HKLIN
        if not hklin.isSet():
            return []
        return get_groups_of_type(str(hklin.fullPath), self.GROUP_TYPE)

    def runTimeValidity(self):
        from ccp4i2.core import CCP4ErrorHandling

        error = super().runTimeValidity()
        if error.maxSeverity() >= CCP4ErrorHandling.SEVERITY_ERROR:
            return error
        if not self.container.inputData.HKLIN.isSet():
            return error

        candidates = self._candidate_groups()
        group, status = resolve_group(candidates, self._requested_labels())
        columns_name = self.TASKNAME + '.container.controlParameters.COLUMNS'
        hklin_name = self.TASKNAME + '.container.inputData.HKLIN'

        if status == 'none':
            error.append(
                self.TASKNAME, 211,
                'No ' + self.TYPE_LABEL + ' columns were found in this MTZ file.',
                hklin_name, CCP4ErrorHandling.SEVERITY_ERROR,
            )
        elif status == 'ambiguous':
            error.append(
                self.TASKNAME, 212,
                'This MTZ has several ' + self.TYPE_LABEL + ' column groups: '
                + describe_candidates(candidates)
                + '. Enter the columns to import in the Columns field.',
                columns_name, CCP4ErrorHandling.SEVERITY_ERROR,
            )
        elif status == 'no_match':
            error.append(
                self.TASKNAME, 213,
                'The columns you specified did not match any '
                + self.TYPE_LABEL + ' group. Available: '
                + describe_candidates(candidates),
                columns_name, CCP4ErrorHandling.SEVERITY_ERROR,
            )
        return error

    def _set_output_metadata(self, out_obj, content_flag, labels, dataset):
        """Set contentFlag / subType / annotation on the output object."""
        if content_flag is not None:
            out_obj.contentFlag.set(content_flag)
        annotation = (self.TYPE_LABEL[:1].upper() + self.TYPE_LABEL[1:]
                      + ' from ' + Path(str(self.container.inputData.HKLIN.baseName)).stem)
        prefix = (dataset + '/') if dataset else ''
        annotation += ' columns ' + prefix + '[' + ','.join(labels) + ']'
        try:
            out_obj.annotation.set(annotation)
        except Exception:
            out_obj.annotation = annotation

    def startProcess(self):
        inp = self.container.inputData
        out = self.container.outputData

        hklin = str(inp.HKLIN.fullPath)
        if not inp.HKLIN.isSet() or not Path(hklin).exists():
            self.appendErrorReport(210, hklin)
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        candidates = self._candidate_groups()
        group, status = resolve_group(candidates, self._requested_labels())
        if group is None:
            code = {'none': 211, 'ambiguous': 212, 'no_match': 213}.get(status, 211)
            self.appendErrorReport(code, describe_candidates(candidates))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        labels = group_labels(group)
        content_flag = group_content_flag(group)
        dataset = group_dataset(group)
        mapping = canonical_mapping(self.OUTPUT_CLASS, content_flag, labels)

        out_obj = getattr(out, self.OUTPUT_PARAM)
        output_path = str(out_obj.fullPath)

        try:
            split_mtz_file(hklin, output_path, mapping)
        except Exception as e:
            self.appendErrorReport(214, str(e))
            self.reportStatus(CPluginScript.FAILED)
            return CPluginScript.FAILED

        self._set_output_metadata(out_obj, content_flag, labels, dataset)
        self.finish_output(out_obj, group, labels)

        self.reportStatus(CPluginScript.SUCCEEDED)
        return CPluginScript.SUCCEEDED

    def finish_output(self, out_obj, group, labels):
        """Subclass hook for type-specific metadata (e.g. map subType)."""
        return None
