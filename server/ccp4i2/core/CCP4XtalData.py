"""
Implementation classes for CCP4XtalData.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""

from typing import Optional, Any

from ccp4i2.core.CCP4Data import CFloatRange
from ccp4i2.core.base_object.error_reporting import CErrorReport
from ccp4i2.core.cdata_stubs.CCP4XtalData import CAltSpaceGroupStub, CAltSpaceGroupListStub, CAnomalousColumnGroupStub, CAnomalousIntensityColumnGroupStub, CAnomalousScatteringElementStub, CAsuComponentStub, CAsuComponentListStub, CCellStub, CCellAngleStub, CCellLengthStub, CColumnGroupStub, CColumnGroupItemStub, CColumnGroupListStub, CColumnTypeStub, CColumnTypeListStub, CCrystalNameStub, CDatasetStub, CDatasetListStub, CDatasetNameStub, CDialsJsonFileStub, CDialsPickleFileStub, CExperimentalDataTypeStub, CFPairColumnGroupStub, CFSigFColumnGroupStub, CFormFactorStub, CFreeRColumnGroupStub, CFreeRDataFileStub, CGenericReflDataFileStub, CHLColumnGroupStub, CIPairColumnGroupStub, CISigIColumnGroupStub, CImageFileStub, CImageFileListStub, CImosflmXmlDataFileStub, CImportUnmergedStub, CImportUnmergedListStub, CMapCoeffsDataFileStub, CMapColumnGroupStub, CMapDataFileStub, CMergeMiniMtzStub, CMergeMiniMtzListStub, CMiniMtzDataFileStub, CMiniMtzDataFileListStub, CMmcifReflDataStub, CMmcifReflDataFileStub, CMtzColumnStub, CMtzColumnGroupStub, CMtzColumnGroupTypeStub, CMtzDataStub, CMtzDataFileStub, CMtzDatasetStub, CObsDataFileStub, CPhaserRFileDataFileStub, CPhaserSolDataFileStub, CPhiFomColumnGroupStub, CPhsDataFileStub, CProgramColumnGroupStub, CProgramColumnGroup0Stub, CRefmacKeywordFileStub, CReindexOperatorStub, CResolutionRangeStub, CRunBatchRangeStub, CRunBatchRangeListStub, CShelxFADataFileStub, CShelxLabelStub, CSpaceGroupStub, CSpaceGroupCellStub, CUnmergedDataContentStub, CUnmergedDataFileStub, CUnmergedDataFileListStub, CUnmergedMtzDataFileStub, CWavelengthStub, CXia2ImageSelectionStub, CXia2ImageSelectionListStub


def _compact_batch_ranges(numbers: list) -> str:
    """
    Convert a list of batch numbers to compact range notation.

    Example: [1, 2, 3, 4, 5, 302, 303, 304, 305] -> "1-5, 302-305"

    Args:
        numbers: Sorted list of batch numbers

    Returns:
        String with compact range notation
    """
    if not numbers:
        return ""

    ranges = []
    start = numbers[0]
    end = numbers[0]

    for num in numbers[1:]:
        if num == end + 1:
            # Extend current range
            end = num
        else:
            # End current range and start new one
            if start == end:
                ranges.append(str(start))
            else:
                ranges.append(f"{start}-{end}")
            start = num
            end = num

    # Add the last range
    if start == end:
        ranges.append(str(start))
    else:
        ranges.append(f"{start}-{end}")

    return ", ".join(ranges)


class CAltSpaceGroup(CAltSpaceGroupStub):
    """
    A string holding the space group
    
    Extends CAltSpaceGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAltSpaceGroupList(CAltSpaceGroupListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CAltSpaceGroupListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAnomalousColumnGroup(CAnomalousColumnGroupStub):
    """
    Selection of F/I and AnomF/I columns from MTZ.
Expected to be part of ab initio phasing dataset ( CDataset)
    
    Extends CAnomalousColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAnomalousIntensityColumnGroup(CAnomalousIntensityColumnGroupStub):
    """
    Selection of I and AnomI columns from MTZ.
Expected to be part of ab initio phasing dataset ( CDataset)
    
    Extends CAnomalousIntensityColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAnomalousScatteringElement(CAnomalousScatteringElementStub):
    """
    Definition of a anomalous scattering element
    
    Extends CAnomalousScatteringElementStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAsuComponent(CAsuComponentStub):
    """
    A component of the asymmetric unit. This is for use in MR, defining
what we are searching for. 
    
    Extends CAsuComponentStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAsuComponentList(CAsuComponentListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CAsuComponentListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CCell(CCellStub):
    """
    A unit cell

    Extends CCellStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def validity(self):
        """Validate the cell parameters.

        If the cell is optional (allowUndefined=True or not specified) and not set,
        skip validation of child parameters (a, b, c, alpha, beta, gamma) entirely.
        This prevents validation errors for unset optional cell parameters.

        Returns:
            CErrorReport containing any validation errors/warnings
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport

        # Check if this cell is optional
        allow_undefined = self.get_qualifier('allowUndefined')
        is_optional = allow_undefined is None or allow_undefined is True

        # Check if the cell has any parameters set
        has_any_set = any(
            getattr(self, param) and getattr(self, param).isSet()
            for param in ['a', 'b', 'c', 'alpha', 'beta', 'gamma']
        )

        # If optional and nothing is set, skip validation entirely
        if is_optional and not has_any_set:
            return CErrorReport()

        # Otherwise, use the standard validation (which validates children)
        return super().validity()

    def guiLabel(self) -> str:
        """
        Return a compact string representation of cell parameters.
        Format: "a=XX.XX b=XX.XX c=XX.XX α=XX.XX β=XX.XX γ=XX.XX"
        Used by legacy code for logging and display.

        Returns:
            Formatted string with cell parameters (2 decimal places)
        """
        parts = []
        if self.a and self.a.isSet():
            parts.append(f"a={self.a.value:.2f}")
        if self.b and self.b.isSet():
            parts.append(f"b={self.b.value:.2f}")
        if self.c and self.c.isSet():
            parts.append(f"c={self.c.value:.2f}")
        if self.alpha and self.alpha.isSet():
            parts.append(f"α={self.alpha.value:.2f}")
        if self.beta and self.beta.isSet():
            parts.append(f"β={self.beta.value:.2f}")
        if self.gamma and self.gamma.isSet():
            parts.append(f"γ={self.gamma.value:.2f}")

        return " ".join(parts) if parts else "Cell parameters not set"


class CCellAngle(CCellAngleStub):
    """
    A cell angle
    
    Extends CCellAngleStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CCellLength(CCellLengthStub):
    """
    A cell length
    
    Extends CCellLengthStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CColumnGroup(CColumnGroupStub):
    """
    Groups of columns in MTZ - probably from analysis by hklfile

    Extends CColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)
        # Set subItem qualifier on columnList so CList.set() knows to convert
        # dict items to CMtzColumn objects. Accessing columnList will auto-create
        # it if not yet instantiated.
        self.columnList.set_qualifier('subItem', {'class': CMtzColumn})

    def columnListStr(self, withTypes: bool = True, splitter: str = ',') -> str:
        """
        Return a string representation of column labels in this group.

        Args:
            withTypes: If True, append column type in parentheses (e.g., "F(F)")
            splitter: Separator between columns (default: ',')

        Returns:
            String like "F,SIGF" or "F(F),SIGF(Q)" depending on withTypes
        """
        column_list = self.columnList
        if len(column_list) == 0:
            return ''

        parts = []
        for col in column_list:
            label = str(col.columnLabel) if col.columnLabel else ''
            if withTypes and col.columnType:
                label = label + '(' + str(col.columnType) + ')'
            parts.append(label)

        return splitter.join(parts)


class CColumnGroupItem(CColumnGroupItemStub):
    """
    Definition of set of columns that form a 'group'
    
    Extends CColumnGroupItemStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CColumnGroupList(CColumnGroupListStub):
    """
    A list with all items of one CData sub-class

    Extends CColumnGroupListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    SUBITEM = {'class': CColumnGroup}


class CColumnType(CColumnTypeStub):
    """
    A list of recognised MTZ column types
    
    Extends CColumnTypeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CColumnTypeList(CColumnTypeListStub):
    """
    A list of acceptable MTZ column types
    
    Extends CColumnTypeListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CCrystalName(CCrystalNameStub):
    """
    A string
    
    Extends CCrystalNameStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDataset(CDatasetStub):
    """
    The experimental data model for ab initio phasing
    
    Extends CDatasetStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDatasetList(CDatasetListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CDatasetListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDatasetName(CDatasetNameStub):
    """
    A string
    
    Extends CDatasetNameStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDialsJsonFile(CDialsJsonFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CDialsJsonFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDialsPickleFile(CDialsPickleFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CDialsPickleFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CExperimentalDataType(CExperimentalDataTypeStub):
    """
    Experimental data type e.g. native or peak
    
    Extends CExperimentalDataTypeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CFPairColumnGroup(CFPairColumnGroupStub):
    """
    A group of MTZ columns required for program input
    
    Extends CFPairColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CFSigFColumnGroup(CFSigFColumnGroupStub):
    """
    A group of MTZ columns required for program input
    
    Extends CFSigFColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CFormFactor(CFormFactorStub):
    """
    The for factor (Fp and Fpp) for a giving element and wavelength
    
    Extends CFormFactorStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CFreeRColumnGroup(CFreeRColumnGroupStub):
    """
    A group of MTZ columns required for program input
    
    Extends CFreeRColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMiniMtzDataFile(CMiniMtzDataFileStub):
    """
    An MTZ experimental data file

    Extends CMiniMtzDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        super().__init__(file_path=file_path, parent=parent, name=name, **kwargs)
        # Note: fileContentClassName='CMtzData' is already set in CMiniMtzDataFileStub decorator

    # Note: _get_conversion_output_path() is now in CDataFile base class

    def columnNames(self, asString=False, ifString=None):
        """
        Get column names from the file content, or expected names based on contentFlag.

        Legacy compatibility method that inspects fileContent to return column names.
        If fileContent is not available, falls back to CONTENT_SIGNATURE_LIST based
        on the contentFlag.

        Args:
            asString: If True, return comma-separated string. If False, return list.
            ifString: Legacy parameter name for asString (for backwards compatibility)

        Returns:
            list or str: Column names as list or comma-separated string

        Example:
            >>> mtz_file.columnNames(False)
            ['H', 'K', 'L', 'F', 'SIGF']
            >>> mtz_file.columnNames(True)
            'H,K,L,F,SIGF'
        """
        # Support legacy ifString parameter
        if ifString is not None:
            asString = ifString

        col_names = []

        # Try to get from fileContent.listOfColumns
        content = self.getFileContent()
        if content and hasattr(content, 'listOfColumns') and content.listOfColumns:
            # Extract column labels from CMtzColumn objects
            for col in content.listOfColumns:
                if hasattr(col, 'columnLabel'):
                    label = col.columnLabel
                    # Handle both CString with .value and plain strings
                    if hasattr(label, 'value'):
                        col_names.append(str(label.value))
                    else:
                        col_names.append(str(label))

        # Fall back to CONTENT_SIGNATURE_LIST if no fileContent columns
        if not col_names and hasattr(self.__class__, 'CONTENT_SIGNATURE_LIST'):
            content_flag = None
            if hasattr(self, 'contentFlag') and self.contentFlag:
                cf = self.contentFlag
                content_flag = cf.value if hasattr(cf, 'value') else int(cf)

            signature_list = self.__class__.CONTENT_SIGNATURE_LIST
            if content_flag is not None and 1 <= content_flag <= len(signature_list):
                # contentFlag is 1-indexed
                col_names = list(signature_list[content_flag - 1])

        if asString:
            return ','.join(col_names)
        else:
            return col_names

    def _introspect_content_flag(self) -> Optional[int]:
        """
        Auto-detect content flag by reading MTZ columns and matching against
        CONTENT_SIGNATURE_LIST.

        This method uses gemmi to read the MTZ file and extract column labels,
        then compares them against the class's CONTENT_SIGNATURE_LIST to
        determine the appropriate contentFlag value.

        Returns:
            Detected content flag (1-indexed), or None if:
            - File doesn't exist
            - File cannot be read
            - No signature matches the columns
            - Class has no CONTENT_SIGNATURE_LIST

        Example:
            >>> obs_file = CObsDataFile('/path/to/file.mtz')
            >>> obs_file.setContentFlag()  # Auto-detects from file
            >>> print(obs_file.contentFlag)  # 4 (if columns are F, SIGF)
        """
        from pathlib import Path

        # Check if file exists
        file_path = self.getFullPath()
        if not file_path or not Path(file_path).exists():
            return None

        # Check if class has CONTENT_SIGNATURE_LIST
        if not hasattr(self.__class__, 'CONTENT_SIGNATURE_LIST'):
            return None

        signature_list = self.__class__.CONTENT_SIGNATURE_LIST

        try:
            import gemmi

            # Read MTZ file
            mtz = gemmi.read_mtz_file(file_path)

            # Extract column labels (just the names, not types)
            column_labels = [col.label for col in mtz.columns]
            column_set = set(column_labels)

            # Match against signatures
            for idx, required_columns in enumerate(signature_list):
                required_set = set(required_columns)

                # Check if all required columns are present
                if required_set.issubset(column_set):
                    # Match found: return contentFlag (1-indexed)
                    return idx + 1

            # No match found
            return None

        except Exception:
            # File reading error or gemmi not available
            return None

    def datasetName(self) -> str:
        """
        Get the first dataset name from the MTZ file (excluding 'HKL_base').

        Returns:
            Dataset name, or empty string if:
            - File is not set
            - File doesn't exist
            - No datasets found (excluding 'HKL_base')
            - File cannot be read
        """
        if not self.isSet():
            return ''

        file_path = self.getFullPath()
        if not file_path:
            return ''

        try:
            import gemmi
            from pathlib import Path

            if not Path(file_path).exists():
                return ''

            # Read MTZ file
            mtz = gemmi.read_mtz_file(file_path)

            # Get dataset names (excluding HKL_base)
            for dataset in mtz.datasets:
                dataset_name = dataset.dataset_name
                if dataset_name and dataset_name != 'HKL_base':
                    return dataset_name

            return ''

        except Exception:
            # File reading error or gemmi not available
            return ''

    def fileExtensions(self):
        """
        Return file extension for MTZ files.

        All MTZ files (CMiniMtzDataFile and derivatives) use the .mtz extension.

        Returns:
            list: ['mtz']
        """
        return ['mtz']

class CFreeRDataFile(CFreeRDataFileStub, CMiniMtzDataFile):
    """
    An MTZ experimental data file for free-R flags.

    Inherits from CMiniMtzDataFile to gain MTZ-specific methods:
    - columnNames(): Get column names from file content
    - _introspect_content_flag(): Auto-detect content flag from MTZ columns
    - datasetName(): Get dataset name from MTZ
    - fileExtensions(): Return ['mtz']

    Extends CFreeRDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Standard column signature for FreeR files
    # splitHklout() will automatically relabel non-standard names (e.g., 'FreeR_flag') to 'FREER'
    CONTENT_SIGNATURE_LIST = [['FREER']]

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        super().__init__(file_path=file_path, parent=parent, name=name, **kwargs)
        # Note: fileContentClassName='CMtzData' is already set in CFreeRDataFileStub decorator


class CGenericReflDataFile(CGenericReflDataFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    Extends CGenericReflDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def getFormat(self):
        """
        Detect file format from file extension.

        Legacy API compatibility method for import_merged and aimless plugins.

        Returns:
            str: File format ('mtz', 'mmcif', 'cif', or 'unknown')
        """
        # Get the file path
        file_path = self.getFullPath()
        if not file_path:
            # Try baseName if getFullPath() returns None
            if hasattr(self, 'baseName') and self.baseName is not None:
                file_path = str(self.baseName.value if hasattr(self.baseName, 'value') else self.baseName)
            else:
                return 'unknown'

        # Detect format from extension
        file_path_lower = str(file_path).lower()
        if file_path_lower.endswith('.mtz'):
            return 'mtz'
        elif file_path_lower.endswith('.cif') or file_path_lower.endswith('.mmcif') or file_path_lower.endswith('.ent'):
            # .ent files are PDB entry files (modern format is mmCIF)
            return 'mmcif'
        else:
            return 'unknown'

    def getMerged(self):
        """
        Check if reflection data is merged.

        Legacy API compatibility method for import_merged plugin.

        Returns:
            bool: True if data is merged, False if unmerged
        """
        # Default to True (merged) unless we can determine otherwise
        # In a full implementation, this would inspect the MTZ/mmCIF file
        # to check for unmerged data indicators (e.g., I+/I- columns)
        return True


class CHLColumnGroup(CHLColumnGroupStub):
    """
    A group of MTZ columns required for program input
    
    Extends CHLColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CIPairColumnGroup(CIPairColumnGroupStub):
    """
    A group of MTZ columns required for program input
    
    Extends CIPairColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CISigIColumnGroup(CISigIColumnGroupStub):
    """
    A group of MTZ columns required for program input
    
    Extends CISigIColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CImageFile(CImageFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CImageFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CImageFileList(CImageFileListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CImageFileListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CImosflmXmlDataFile(CImosflmXmlDataFileStub):
    """
    An iMosflm data file
    
    Extends CImosflmXmlDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CImportUnmerged(CImportUnmergedStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None

    Extends CImportUnmergedStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def populateFromFile(self):
        """
        Auto-populate dataset, crystalName, cell, and wavelength from MTZ file metadata.

        This should be called after the file is set to populate required fields
        that can be derived from the MTZ header.

        Returns:
            bool: True if fields were populated successfully, False otherwise
        """
        if not hasattr(self, 'file') or self.file is None or not self.file.isSet():
            return False

        file_path = self.file.getFullPath()
        if not file_path:
            return False

        try:
            import gemmi
            from pathlib import Path

            if not Path(file_path).exists():
                return False

            mtz = gemmi.read_mtz_file(file_path)

            # Find the first non-HKL_base dataset
            for dataset in mtz.datasets:
                if dataset.dataset_name and dataset.dataset_name != 'HKL_base':
                    # Set dataset name if not already set
                    if hasattr(self, 'dataset') and self.dataset is not None:
                        if not self.dataset.isSet():
                            self.dataset.set(dataset.dataset_name)

                    # Set crystal name if not already set
                    if hasattr(self, 'crystalName') and self.crystalName is not None:
                        if not self.crystalName.isSet():
                            crystal_name = dataset.crystal_name or dataset.dataset_name
                            self.crystalName.set(crystal_name)

                    # Set cell parameters if not already set
                    if hasattr(self, 'cell') and self.cell is not None:
                        if not self.cell.isSet():
                            self.cell.a.set(mtz.cell.a)
                            self.cell.b.set(mtz.cell.b)
                            self.cell.c.set(mtz.cell.c)
                            self.cell.alpha.set(mtz.cell.alpha)
                            self.cell.beta.set(mtz.cell.beta)
                            self.cell.gamma.set(mtz.cell.gamma)

                    # Set wavelength if available and not already set
                    if hasattr(self, 'wavelength') and self.wavelength is not None:
                        if not self.wavelength.isSet() and dataset.wavelength > 0:
                            self.wavelength.set(dataset.wavelength)

                    return True

            # Fallback: use cell from MTZ even if no named dataset found
            if hasattr(self, 'cell') and self.cell is not None and not self.cell.isSet():
                self.cell.a.set(mtz.cell.a)
                self.cell.b.set(mtz.cell.b)
                self.cell.c.set(mtz.cell.c)
                self.cell.alpha.set(mtz.cell.alpha)
                self.cell.beta.set(mtz.cell.beta)
                self.cell.gamma.set(mtz.cell.gamma)

            return False

        except Exception as e:
            import logging
            logger = logging.getLogger(__name__)
            logger.debug(f"Failed to populate CImportUnmerged from file: {e}")
            return False


class CImportUnmergedList(CImportUnmergedListStub):
    """
    A list with all items of one CData sub-class

    Extends CImportUnmergedListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)
        # Set the subItem qualifier that makeItem() expects
        # This allows i2run to create CImportUnmerged items when parsing arguments
        self.set_qualifier('subItem', {'class': CImportUnmerged, 'qualifiers': {}})


class CMapCoeffsDataFile(CMapCoeffsDataFileStub, CMiniMtzDataFile):
    """
    An MTZ map coefficients data file (FPHI format).

    Inherits from CMiniMtzDataFile to gain MTZ-specific methods:
    - columnNames(): Get column names from file content
    - _introspect_content_flag(): Auto-detect content flag from MTZ columns
    - datasetName(): Get dataset name from MTZ
    - fileExtensions(): Return ['mtz']

    Extends CMapCoeffsDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def as_FPHI(self, work_directory: Optional[Any] = None) -> str:
        """
        Return path to this file as FPHI format.

        Since CMapCoeffsDataFile only has one content type (FPHI),
        this method simply returns the current file path without conversion.

        FPHI format: F, PHI

        This is a thin wrapper around PhaseDataConverter.to_fphi().
        See core.conversions.phase_data_converter for implementation details.

        Args:
            work_directory: Ignored for this class

        Returns:
            Full path to this file (no conversion needed)
        """
        from ccp4i2.core.conversions import PhaseDataConverter
        return PhaseDataConverter.to_fphi(self, work_directory=work_directory)


class CPhsDataFile(CPhsDataFileStub, CMiniMtzDataFile):
    """
    An MTZ phase data file.

    Handles phase data in different formats:
    - HL (1): Hendrickson-Lattman coefficients (HLA, HLB, HLC, HLD)
    - PHIFOM (2): Phase + Figure of Merit (PHI, FOM)

    Inherits from CMiniMtzDataFile to gain MTZ-specific methods:
    - columnNames(): Get column names from file content
    - _introspect_content_flag(): Auto-detect content flag from MTZ columns
    - datasetName(): Get dataset name from MTZ
    - fileExtensions(): Return ['mtz']

    Extends CPhsDataFileStub with conversion methods for transforming
    between different phase data representations.
    """

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        super().__init__(file_path=file_path, parent=parent, name=name, **kwargs)
        # Note: fileContentClassName='CMtzData' is already set in CPhsDataFileStub decorator

    # Note: setContentFlag() inherited from CDataFile base class
    # It uses _introspect_content_flag() from CMiniMtzDataFile, which reads
    # CONTENT_SIGNATURE_LIST = [['HLA', 'HLB', 'HLC', 'HLD'], ['PHI', 'FOM']]
    # defined in CPhsDataFileStub

    def as_HL(self, work_directory: Optional[Any] = None) -> str:
        """
        Convert this file to HL format (Hendrickson-Lattman coefficients).

        HL format: HLA, HLB, HLC, HLD

        This is a thin wrapper around PhaseDataConverter.to_hl().
        See core.conversions.phase_data_converter for implementation details.

        Args:
            work_directory: Directory for output if input dir not writable

        Returns:
            Full path to converted file

        Raises:
            ValueError: If conversion not possible from current format
        """
        from ccp4i2.core.conversions import PhaseDataConverter
        return PhaseDataConverter.to_hl(self, work_directory=work_directory)

    def as_PHIFOM(self, work_directory: Optional[Any] = None) -> str:
        """
        Convert this file to PHIFOM format (Phase + Figure of Merit).

        PHIFOM format: PHI, FOM

        This is a thin wrapper around PhaseDataConverter.to_phifom().
        See core.conversions.phase_data_converter for implementation details.

        Args:
            work_directory: Directory for output if input dir not writable

        Returns:
            Full path to converted file

        Raises:
            ValueError: If conversion not possible from current format
        """
        from ccp4i2.core.conversions import PhaseDataConverter
        return PhaseDataConverter.to_phifom(self, work_directory=work_directory)


class CMapColumnGroup(CMapColumnGroupStub):
    """
    A group of MTZ columns required for program input
    
    Extends CMapColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMapDataFile(CMapDataFileStub):
    """
    A CCP4 Map file
    
    Extends CMapDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMergeMiniMtz(CMergeMiniMtzStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CMergeMiniMtzStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMergeMiniMtzList(CMergeMiniMtzListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CMergeMiniMtzListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass



class CMiniMtzDataFileList(CMiniMtzDataFileListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CMiniMtzDataFileListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMmcifReflData(CMmcifReflDataStub):
    """
    Generic mmCIF data.
This is intended to be a base class for other classes
specific to coordinates, reflections or geometry data.

    Extends CMmcifReflDataStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def loadFile(self, file_path: str = None):
        """
        Load mmCIF reflection file using gemmi library.

        This method:
        1. Reads mmCIF file using gemmi.cif.read()
        2. Extracts cell, space group, wavelength
        3. Detects which column types are present
        4. Uses proper CData setters (NO __dict__ manipulation)

        Args:
            file_path: Optional path to mmCIF file. If None, gets path from parent CDataFile.

        Returns:
            CErrorReport with any errors encountered

        Example:
            # Load from parent file's path
            >>> cif_file = CMmcifReflDataFile()
            >>> cif_file.setFullPath('/path/to/reflections.cif')
            >>> cif_file.fileContent.loadFile()

            # Load from explicit path (legacy pattern)
            >>> cif_data = CMmcifReflData()
            >>> error = cif_data.loadFile('/path/to/reflections.cif')
            >>> if error.count() == 0:
            ...     print(f"Has F columns: {cif_data.haveFobsColumn.value}")
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        from pathlib import Path

        error = CErrorReport()

        # If no path provided, get from parent CDataFile
        if file_path is None:
            parent = self.get_parent()
            if parent is not None and hasattr(parent, 'getFullPath'):
                file_path = parent.getFullPath()

        # Validate file path
        if not file_path:
            return error

        path_obj = Path(file_path)
        if not path_obj.exists() or not path_obj.is_file():
            error.append(
                klass=self.__class__.__name__,
                code=101,
                details=f"mmCIF file does not exist or is not a file: '{file_path}'",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        try:
            import gemmi
        except ImportError as e:
            error.append(
                klass=self.__class__.__name__,
                code=102,
                details=f"Failed to import gemmi library: {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        try:
            # Read mmCIF file using gemmi
            cif_doc = gemmi.cif.read(str(file_path))

            # Get the first data block
            if len(cif_doc) == 0:
                error.append(
                    klass=self.__class__.__name__,
                    code=104,
                    details=f"mmCIF file contains no data blocks: '{file_path}'",
                    name=self.object_name() if hasattr(self, 'object_name') else ''
                )
                return error

            block = cif_doc[0]

            # Extract cell parameters using smart assignment
            if hasattr(self, 'cell') and self.cell is not None:
                try:
                    self.cell.a = float(block.find_value('_cell.length_a'))
                    self.cell.b = float(block.find_value('_cell.length_b'))
                    self.cell.c = float(block.find_value('_cell.length_c'))
                    self.cell.alpha = float(block.find_value('_cell.angle_alpha'))
                    self.cell.beta = float(block.find_value('_cell.angle_beta'))
                    self.cell.gamma = float(block.find_value('_cell.angle_gamma'))
                except (ValueError, RuntimeError):
                    pass  # Cell parameters not available

            # Extract space group using smart assignment
            if hasattr(self, 'spaceGroup') and self.spaceGroup is not None:
                try:
                    sg_name = block.find_value('_symmetry.space_group_name_H-M')
                    if sg_name:
                        self.spaceGroup = sg_name.strip('"').strip("'").strip()
                except RuntimeError:
                    # Try alternative tag
                    try:
                        sg_num = block.find_value('_symmetry.Int_Tables_number')
                        if sg_num:
                            self.spaceGroup = f"Space group {sg_num}"
                    except RuntimeError:
                        pass  # Space group not available

            # Extract wavelength using smart assignment
            if hasattr(self, 'wavelength') and self.wavelength is not None:
                try:
                    wavelength_val = block.find_value('_diffrn_radiation_wavelength.wavelength')
                    if wavelength_val:
                        self.wavelength = float(wavelength_val)
                except (ValueError, RuntimeError):
                    pass  # Wavelength not available

            # Detect which column types are present
            # Check for FreeR column
            if hasattr(self, 'haveFreeRColumn') and self.haveFreeRColumn is not None:
                self.haveFreeRColumn = block.find_loop('_refln.status') is not None

            # Check for F_meas columns
            if hasattr(self, 'haveFobsColumn') and self.haveFobsColumn is not None:
                self.haveFobsColumn = (
                    block.find_loop('_refln.F_meas') is not None or
                    block.find_loop('_refln.F_meas_au') is not None
                )

            # Check for F+/F- columns
            if hasattr(self, 'haveFpmObsColumn') and self.haveFpmObsColumn is not None:
                self.haveFpmObsColumn = block.find_loop('_refln.pdbx_F_plus') is not None

            # Check for intensity columns
            if hasattr(self, 'haveIobsColumn') and self.haveIobsColumn is not None:
                self.haveIobsColumn = (
                    block.find_loop('_refln.intensity_meas') is not None or
                    block.find_loop('_refln.F_squared_meas') is not None
                )

            # Check for I+/I- columns
            if hasattr(self, 'haveIpmObsColumn') and self.haveIpmObsColumn is not None:
                self.haveIpmObsColumn = block.find_loop('_refln.pdbx_I_plus') is not None

            # Store gemmi CIF document for advanced queries
            # Use object.__setattr__ to bypass smart assignment
            object.__setattr__(self, '_gemmi_cif_doc', cif_doc)

            # Emit signal if available
            if hasattr(self, 'dataLoaded'):
                try:
                    self.dataLoaded.emit()
                except:
                    pass  # Signal system may not be available

        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=103,
                details=f"Error reading mmCIF file '{file_path}': {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )

        return error


class CMmcifReflDataFile(CMmcifReflDataFileStub):
    """
    A reflection file in mmCIF format

    Extends CMmcifReflDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        super().__init__(file_path=file_path, parent=parent, name=name, **kwargs)
        # Note: fileContentClassName='CMmcifReflData' is already set in CMmcifReflDataFileStub decorator


class CMtzColumn(CMtzColumnStub):
    """
    An MTZ column with column label and column type
    
    Extends CMtzColumnStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMtzColumnGroup(CMtzColumnGroupStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CMtzColumnGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMtzColumnGroupType(CMtzColumnGroupTypeStub, CColumnType):
    """
    
    Inherits from:
    - CMtzColumnGroupTypeStub: Metadata and structure
    - CColumnType: Shared full-fat methods
    A list of recognised MTZ column types
    
    Extends CMtzColumnGroupTypeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CMtzData(CMtzDataStub):
    """
    Base class for classes holding file contents

    Extends CMtzDataStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)
        # Initialize sub-objects required for loadFile()
        if self.cell is None:
            from ccp4i2.core.cdata_stubs.CCP4XtalData import CCellStub
            self.cell = CCellStub(parent=self, name='cell')
        if self.spaceGroup is None:
            from ccp4i2.core.cdata_stubs.CCP4XtalData import CSpaceGroupStub
            self.spaceGroup = CSpaceGroupStub(parent=self, name='spaceGroup')
        if self.resolutionRange is None:
            from ccp4i2.core.cdata_stubs.CCP4XtalData import CResolutionRangeStub
            self.resolutionRange = CResolutionRangeStub(parent=self, name='resolutionRange')

    def to_dict(self):
        """
        Convert CMtzData to a dictionary for serialization.

        Includes all standard children plus type-specific attributes like
        datasets, wavelengths, crystalNames, listOfColumns, and datasetCells.
        """
        from ccp4i2.lib.utils.parameters.value_dict import value_dict_for_object

        # Start with base serialization from children
        result = {}
        try:
            for child in self.children():
                if hasattr(child, 'objectName') and callable(child.objectName):
                    name = child.objectName()
                elif hasattr(child, 'object_name') and callable(child.object_name):
                    name = child.object_name()
                elif hasattr(child, 'name'):
                    name = child.name
                else:
                    continue
                result[name] = value_dict_for_object(child)
        except Exception:
            pass

        # Add type-specific attributes that may be plain Python lists
        type_specific_attrs = [
            'datasets', 'wavelengths', 'crystalNames', 'datasetCells', 'listOfColumns'
        ]
        for attr_name in type_specific_attrs:
            if attr_name not in result and hasattr(self, attr_name):
                try:
                    attr_value = getattr(self, attr_name)
                    if attr_value is not None:
                        if isinstance(attr_value, list) and len(attr_value) > 0:
                            result[attr_name] = value_dict_for_object(attr_value)
                        elif hasattr(attr_value, 'value') and attr_value.value:
                            result[attr_name] = value_dict_for_object(attr_value)
                        elif hasattr(attr_value, '__iter__'):
                            # Handle CList or other iterables
                            converted = value_dict_for_object(attr_value)
                            if converted:
                                result[attr_name] = converted
                except Exception:
                    pass

        return result

    def loadFile(self, file_path: str = None):
        """
        Load MTZ file using gemmi library.

        This method:
        1. Reads MTZ file using gemmi.read_mtz_file()
        2. Extracts metadata (cell, spacegroup, resolution, datasets)
        3. Populates CData attributes using proper setters (NO __dict__ manipulation)
        4. Stores gemmi Mtz object for advanced queries

        Args:
            file_path: Optional path to MTZ file. If None, gets path from parent CDataFile.

        Returns:
            CErrorReport with any errors encountered

        Example:
            # Load from parent file's path
            >>> mtz_file = CMtzDataFile()
            >>> mtz_file.setFullPath('/path/to/data.mtz')
            >>> mtz_file.fileContent.loadFile()

            # Load from explicit path (legacy pattern)
            >>> mtz_data = CMtzData()
            >>> error = mtz_data.loadFile('/path/to/data.mtz')
            >>> if error.count() == 0:
            ...     print(f"Space group: {mtz_data.spaceGroup}")
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        from pathlib import Path

        error = CErrorReport()

        # If no path provided, get from parent CDataFile
        if file_path is None:
            parent = self.get_parent()
            if parent is not None and hasattr(parent, 'getFullPath'):
                file_path = parent.getFullPath()

        # Validate file path
        if not file_path:
            return error

        path_obj = Path(file_path)
        if not path_obj.exists() or not path_obj.is_file():
            error.append(
                klass=self.__class__.__name__,
                code=101,
                details=f"MTZ file does not exist or is not a file: '{file_path}'",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        try:
            import gemmi
        except ImportError as e:
            error.append(
                klass=self.__class__.__name__,
                code=102,
                details=f"Failed to import gemmi library: {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        try:
            # Read MTZ file using gemmi
            mtz = gemmi.read_mtz_file(str(file_path))

            # Extract cell parameters using smart assignment (NO __dict__)
            if hasattr(self, 'cell') and self.cell is not None:
                self.cell.a = mtz.cell.a
                self.cell.b = mtz.cell.b
                self.cell.c = mtz.cell.c
                self.cell.alpha = mtz.cell.alpha
                self.cell.beta = mtz.cell.beta
                self.cell.gamma = mtz.cell.gamma

            # Extract space group using smart assignment
            if hasattr(self, 'spaceGroup') and self.spaceGroup is not None:
                self.spaceGroup = mtz.spacegroup.hm

            # Extract resolution range using smart assignment
            if hasattr(self, 'resolutionRange') and self.resolutionRange is not None:
                mtz.update_reso()
                self.resolutionRange.low = mtz.resolution_low()
                self.resolutionRange.high = mtz.resolution_high()

            # Extract column information as CMtzColumn objects
            # Skip H-type columns (Miller indices H, K, L)
            if hasattr(self, 'listOfColumns') and self.listOfColumns is not None:
                # Clear existing columns before loading - prevents duplicates when
                # loadFile() is called on content that already has columns from a copy
                self.listOfColumns.clear()
                # Create CMtzColumn objects with columnLabel and columnType
                mtz_columns = []
                for col in mtz.columns:
                    # Skip Miller indices (type 'H')
                    if col.type == 'H':
                        continue

                    mtz_col = CMtzColumn(name=col.label)
                    mtz_col.columnLabel = col.label
                    mtz_col.columnType = col.type
                    # Get dataset name and use dataset ID as groupIndex
                    if col.dataset:
                        mtz_col.dataset = col.dataset.dataset_name
                        mtz_col.groupIndex = col.dataset.id
                    mtz_columns.append(mtz_col)
                # Populate the CList using append (direct assignment to CList is silently ignored)
                for mtz_col in mtz_columns:
                    self.listOfColumns.append(mtz_col)

            # Extract dataset information
            if hasattr(self, 'datasets') and self.datasets is not None:
                dataset_names = [ds.dataset_name for ds in mtz.datasets]
                self.datasets = dataset_names

            # Extract crystal names
            if hasattr(self, 'crystalNames') and self.crystalNames is not None:
                crystal_names = [ds.crystal_name for ds in mtz.datasets]
                self.crystalNames = crystal_names

            # Extract wavelengths
            if hasattr(self, 'wavelengths') and self.wavelengths is not None:
                wavelength_list = [ds.wavelength for ds in mtz.datasets]
                self.wavelengths = wavelength_list

            # Extract dataset cells
            if hasattr(self, 'datasetCells') and self.datasetCells is not None:
                cells = []
                for ds in mtz.datasets:
                    cells.append({
                        'a': ds.cell.a,
                        'b': ds.cell.b,
                        'c': ds.cell.c,
                        'alpha': ds.cell.alpha,
                        'beta': ds.cell.beta,
                        'gamma': ds.cell.gamma
                    })
                self.datasetCells = cells

            # Determine if merged (MTZ files are typically merged)
            if hasattr(self, 'merged') and self.merged is not None:
                self.merged = True

            # Store gemmi Mtz object for advanced queries
            # Use object.__setattr__ to bypass smart assignment
            object.__setattr__(self, '_gemmi_mtz', mtz)

        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=103,
                details=f"Error reading MTZ file '{file_path}': {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )

        return error

    def getListOfWavelengths(self):
        """
        Get list of wavelengths from MTZ datasets.

        Legacy API compatibility method for refmac and other plugins.
        Returns wavelengths extracted from the MTZ file.

        Returns:
            list: List of wavelengths from datasets (floats)

        Example:
            >>> wavelength = mtz_data.getListOfWavelengths()[-1]  # Get last wavelength
        """
        # First try to use stored gemmi Mtz object (most reliable)
        if hasattr(self, '_gemmi_mtz') and self._gemmi_mtz is not None:
            try:
                wavelengths = [ds.wavelength for ds in self._gemmi_mtz.datasets]
                # Filter out zero/invalid wavelengths
                valid_wavelengths = [w for w in wavelengths if w > 0.0]
                return valid_wavelengths if valid_wavelengths else wavelengths
            except Exception:
                pass

        # Fallback: use wavelengths attribute if set
        if hasattr(self, 'wavelengths') and self.wavelengths is not None:
            if hasattr(self.wavelengths, 'value'):
                wl_list = self.wavelengths.value
            else:
                wl_list = self.wavelengths

            # Handle both list and single value
            if isinstance(wl_list, (list, tuple)):
                return list(wl_list)
            else:
                return [wl_list] if wl_list else []

        # Final fallback: return empty list
        return []

    def getListOfColumns(self):
        """
        Get list of column labels from MTZ file.

        Legacy API compatibility method for TestObsConversions and other plugins.
        Returns list of column label strings.

        Returns:
            list: List of column labels (strings)

        Example:
            >>> columns = mtz_data.getListOfColumns()
            >>> print(columns)  # ['H', 'K', 'L', 'F', 'SIGF', ...]
        """
        # First try to use stored gemmi Mtz object (most reliable)
        if hasattr(self, '_gemmi_mtz') and self._gemmi_mtz is not None:
            try:
                return [col.label for col in self._gemmi_mtz.columns]
            except Exception:
                pass

        # Fallback: use listOfColumns attribute if set
        if hasattr(self, 'listOfColumns') and self.listOfColumns is not None:
            # Handle CList or list of CMtzColumn objects
            columns = []
            col_list = self.listOfColumns.value if hasattr(self.listOfColumns, 'value') else self.listOfColumns

            if isinstance(col_list, (list, tuple)):
                for col in col_list:
                    if hasattr(col, 'columnLabel'):
                        # CMtzColumn object
                        label = col.columnLabel.value if hasattr(col.columnLabel, 'value') else col.columnLabel
                        columns.append(str(label))
                    elif isinstance(col, str):
                        columns.append(col)
                return columns

        # Final fallback: return empty list
        return []

    def getColumnGroups(self):
        """
        Build list of column groups from MTZ column info by pattern matching.

        This method scans through the columns and identifies sequential patterns
        that match known content signatures for CMiniMtzData subtypes:
        - CObsDataFile: FQ (F,sigF), JQ (I,sigI), GLGL (F+,sigF+,F-,sigF-), KMKM (I+,sigI+,I-,sigI-)
        - CMapCoeffsDataFile: FP (F,Phi), FQP (F,sigF,Phi)
        - CPhsDataFile: AAAA (HLA,HLB,HLC,HLD), PW (Phi,FOM)
        - CFreeRDataFile: I (integer flag)

        Returns:
            list: List of CColumnGroup objects with columnList, columnGroupType, contentFlag

        Example:
            >>> groups = mtz_data.getColumnGroups()
            >>> for group in groups:
            ...     print(f"{group.columnGroupType}: {group.columnList}")
        """
        groupList = []

        # Get listOfColumns - handle both CList and plain list
        if not hasattr(self, 'listOfColumns') or self.listOfColumns is None:
            return groupList

        columns = self.listOfColumns.value if hasattr(self.listOfColumns, 'value') else self.listOfColumns
        if not columns:
            return groupList

        # Helper to extract column info
        def get_col_info(col):
            label = col.columnLabel.value if hasattr(col.columnLabel, 'value') else col.columnLabel
            ctype = col.columnType.value if hasattr(col.columnType, 'value') else col.columnType
            dataset = col.dataset.value if hasattr(col.dataset, 'value') else col.dataset
            group_idx = col.groupIndex.value if hasattr(col.groupIndex, 'value') else col.groupIndex
            return str(label), str(ctype), str(dataset) if dataset else '', group_idx

        # Build type string for pattern matching
        col_infos = [get_col_info(col) for col in columns]
        type_string = ''.join(info[1] for info in col_infos)

        # Define signature patterns with their data types and content flags
        # Order matters: check longer patterns first to avoid partial matches
        # Format: (pattern, group_type, content_flag)
        from ccp4i2.core.base_object.class_metadata import get_class_metadata_by_type

        # Build pattern list from class metadata
        patterns = []
        for cls, label in [(CObsDataFile, 'Obs'), (CPhsDataFile, 'Phs'),
                          (CMapCoeffsDataFile, 'MapCoeffs'), (CFreeRDataFile, 'FreeR')]:
            meta = get_class_metadata_by_type(cls)
            correct_columns = meta.qualifiers.get('correctColumns') if meta and meta.qualifiers else None
            if correct_columns:
                for idx, sig in enumerate(correct_columns):
                    patterns.append((sig, label, idx + 1))

        # Sort by pattern length (longest first) to match longer patterns before shorter ones
        patterns.sort(key=lambda x: len(x[0]), reverse=True)

        # Track which columns have been assigned to a group
        used = [False] * len(col_infos)

        # Scan for patterns
        i = 0
        while i < len(col_infos):
            if used[i]:
                i += 1
                continue

            matched = False
            for pattern, group_type, content_flag in patterns:
                pattern_len = len(pattern)
                if i + pattern_len > len(col_infos):
                    continue

                # Check if pattern matches at this position
                candidate = type_string[i:i + pattern_len]
                if candidate == pattern:
                    # Special check for FreeR: label must contain 'free'
                    if group_type == 'FreeR':
                        label = col_infos[i][0]
                        if 'free' not in label.lower():
                            continue

                    # For Obs data, columns must come from the same dataset
                    # (when crystal/dataset names are present in the MTZ)
                    if group_type == 'Obs':
                        first_dataset = col_infos[i][2]
                        if first_dataset:  # Only check if dataset info is present
                            same_dataset = all(
                                col_infos[i + j][2] == first_dataset
                                for j in range(pattern_len)
                            )
                            if not same_dataset:
                                continue

                    # Create a new column group
                    dataset = col_infos[i][2]
                    group = CColumnGroup(dataset=dataset)
                    group.columnGroupType = group_type
                    group.contentFlag = content_flag

                    # Add the matching columns
                    for j in range(pattern_len):
                        label, ctype, ds, grp_idx = col_infos[i + j]
                        group.columnList.append(CMtzColumn(
                            columnLabel=label,
                            columnType=ctype,
                            groupIndex=grp_idx,
                            dataset=ds
                        ))
                        used[i + j] = True

                    groupList.append(group)
                    i += pattern_len
                    matched = True
                    break

            if not matched:
                i += 1

        return groupList

    def clipperSameCell(self, other_content, tolerance=None):
        """
        Compare unit cells using Clipper's reciprocal space algorithm.

        Two cells disagree if the difference in their orthogonalisation matrices
        is sufficient to map a reflection from one cell onto a different reflection
        in the other cell at the given tolerance (resolution in Angstroms).

        This implements the Clipper Cell::equals() algorithm which considers
        reciprocal space vectors and finds the resolution at which reflections
        would be mis-indexed by more than 0.5 reciprocal lattice units.

        Args:
            other_content: Another CMtzData or CUnmergedDataContent instance
            tolerance: Resolution tolerance in Angstroms (default 1.0)
                      Cells are compatible if mis-indexing doesn't occur
                      until this resolution or higher.

        Returns:
            dict with keys:
                'validity': bool - True if cells are compatible at tolerance
                'tolerance': float - the resolution tolerance value
                'difference': float - resolution where mis-indexing starts (Å)
                'maximumResolution1': float - max resolution for cell1
                'maximumResolution2': float - max resolution for cell2
        """
        import numpy as np

        # Use default tolerance if None (Clipper default is 1.0 Angstrom)
        if tolerance is None:
            tolerance = 1.0

        # Get cell parameters
        cell1 = self.cell
        cell2 = other_content.cell

        if not cell1 or not cell2:
            return {
                'validity': False,
                'tolerance': tolerance,
                'difference': float('inf'),
                'maximumResolution1': 0.0,
                'maximumResolution2': 0.0
            }

        # Build orthogonalization matrices for both cells
        # Orthogonalization matrix converts fractional to Cartesian coordinates
        def build_orth_matrix(cell):
            """Build orthogonalization matrix from cell parameters."""
            import math
            a, b, c = float(cell.a), float(cell.b), float(cell.c)
            alpha = math.radians(float(cell.alpha))
            beta = math.radians(float(cell.beta))
            gamma = math.radians(float(cell.gamma))

            # Volume calculation
            cos_alpha = math.cos(alpha)
            cos_beta = math.cos(beta)
            cos_gamma = math.cos(gamma)
            sin_gamma = math.sin(gamma)

            volume = a * b * c * math.sqrt(
                1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2
                + 2 * cos_alpha * cos_beta * cos_gamma
            )

            # Orthogonalization matrix (fractional to Cartesian)
            orth = np.array([
                [a, b * cos_gamma, c * cos_beta],
                [0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
                [0, 0, volume / (a * b * sin_gamma)]
            ])
            return orth

        orth1 = build_orth_matrix(cell1)
        orth2 = build_orth_matrix(cell2)

        # Calculate the difference matrix
        diff_matrix = orth1 - orth2

        # Find the resolution at which a reflection would be mis-indexed
        # by 0.5 reciprocal lattice units or more
        #
        # We test reflections along the crystallographic axes (h00, 0k0, 00l)
        # at increasing resolution until we find where the difference exceeds
        # 0.5 reciprocal lattice units

        max_resolution = 1000.0  # Start with very high resolution (small d)
        axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # Test along a*, b*, c*

        for axis in axes:
            h, k, l = axis
            # Compute reciprocal space vector for this reflection
            hkl = np.array([h, k, l])

            # Transform to Cartesian reciprocal space
            # (reciprocal of orth matrix)
            orth1_inv = np.linalg.inv(orth1).T
            orth2_inv = np.linalg.inv(orth2).T

            # Reciprocal space vectors
            s1 = orth1_inv @ hkl
            s2 = orth2_inv @ hkl

            # Difference in reciprocal space (in Å⁻¹)
            diff = np.linalg.norm(s1 - s2)

            # Resolution of this reflection
            d_spacing = 1.0 / np.linalg.norm(s1)

            # At what resolution would this axis cause 0.5 r.l.u. mis-indexing?
            # diff * d = 0.5 / |hkl|
            # So critical resolution = 0.5 / (diff * |hkl|)
            if diff > 1e-10:  # Avoid division by zero
                critical_res = 0.5 / (diff * np.linalg.norm(hkl))
                max_resolution = min(max_resolution, critical_res)

        # Calculate maximum resolutions based on cell dimensions
        a1, b1, c1 = float(cell1.a), float(cell1.b), float(cell1.c)
        a2, b2, c2 = float(cell2.a), float(cell2.b), float(cell2.c)
        max_res1 = min(a1, b1, c1) / 2.0
        max_res2 = min(a2, b2, c2) / 2.0

        # Cells are compatible if mis-indexing doesn't occur until
        # resolution better than (smaller than) tolerance
        validity = max_resolution <= tolerance

        return {
            'validity': validity,
            'tolerance': tolerance,
            'difference': max_resolution,  # Resolution where mis-indexing starts
            'maximumResolution1': max_res1,
            'maximumResolution2': max_res2
        }

    def matthewsCoeff(self, seqDataFile=None, nRes=None, molWt=None, polymerMode=""):
        """
        Calculate Matthews coefficient using the CCP4 matthews_coef program.

        This provides estimates of the number of molecules in the asymmetric unit
        based on cell volume, space group, and molecular weight.

        Args:
            seqDataFile: Optional sequence data file to get molecular weight from
            nRes: Optional number of residues (will estimate MW as 112.5 * nRes)
            molWt: Optional molecular weight in Daltons (preferred if known)
            polymerMode: Optional polymer mode string for matthews_coef

        Returns:
            dict: Results containing:
                - cell_volume: Unit cell volume
                - results: List of dicts with nmol_in_asu, matth_coef, percent_solvent, prob_matth

        Raises:
            CException: If molecular weight cannot be determined or program fails

        Example:
            >>> mtz = CMtzDataFile()
            >>> mtz.setFullPath('/path/to/data.mtz')
            >>> mtz.loadFile()
            >>> result = mtz.fileContent.matthewsCoeff(molWt=25000)
            >>> for r in result['results']:
            ...     print(f"{r['nmol_in_asu']} copies: {r['percent_solvent']:.1f}% solvent")
        """
        import os
        import tempfile
        import math
        from ccp4i2.core.base_object.error_reporting import CException

        # Determine molecular weight
        if seqDataFile is not None:
            try:
                molWt = seqDataFile.fileContent.getAnalysis('molecularWeight')
            except Exception:
                molWt = 0.0
        elif nRes is not None:
            # Estimated residue weight as per ccp4 matthews_coeff documentation
            molWt = 112.5 * float(nRes)

        if molWt is None or molWt < 0.01:
            raise CException(self.__class__, 410, str(seqDataFile))

        # Create temporary files for output
        f1 = tempfile.mkstemp()
        os.close(f1[0])
        f2 = tempfile.mkstemp()
        os.close(f2[0])

        try:
            # Build command text for matthews_coef
            com_text = f'MOLWEIGHT {molWt}\nCELL'

            # Get cell parameters
            for p in ['a', 'b', 'c']:
                cell_val = getattr(self.cell, p, None)
                if cell_val is not None:
                    if hasattr(cell_val, 'value'):
                        com_text += f' {cell_val.value}'
                    else:
                        com_text += f' {cell_val}'

            for p in ['alpha', 'beta', 'gamma']:
                angle = getattr(self.cell, p, None)
                if angle is not None:
                    if hasattr(angle, 'value'):
                        a = float(angle.value)
                    else:
                        a = float(angle)
                    # Convert from radians if needed
                    if a < 3.0:
                        a = a * 180.0 / math.pi
                    com_text += f' {a}'

            # Get space group number
            sg_number = 1
            if hasattr(self, 'spaceGroup') and self.spaceGroup is not None:
                if hasattr(self.spaceGroup, 'number'):
                    sg_number = self.spaceGroup.number()
                elif hasattr(self.spaceGroup, 'value'):
                    # Try to get number from space group name
                    try:
                        import gemmi
                        sg = gemmi.SpaceGroup(str(self.spaceGroup.value))
                        sg_number = sg.number
                    except Exception:
                        sg_number = 1

            com_text += f'\nSYMM {sg_number}\n'
            com_text += 'XMLO\nAUTO\n'

            if polymerMode:
                com_text += f'\nMODE {polymerMode}'

            # Run matthews_coef using subprocess
            import subprocess
            import shutil

            # Find matthews_coef executable
            ccp4_bin = os.environ.get('CCP4', os.environ.get('CBIN', ''))
            if ccp4_bin:
                matthews_exe = os.path.join(ccp4_bin, 'bin', 'matthews_coef')
            else:
                # Fall back to PATH
                matthews_exe = shutil.which('matthews_coef')

            if not matthews_exe or not os.path.exists(matthews_exe):
                raise CException(self.__class__, 411, 'matthews_coef executable not found')

            # Run the command
            cmd = [matthews_exe, 'XMLFILE', f1[1]]
            with open(f2[1], 'w') as log_file:
                result = subprocess.run(
                    cmd,
                    input=com_text,
                    capture_output=True,
                    text=True,
                    timeout=30
                )
                log_file.write(result.stdout)
                if result.stderr:
                    log_file.write('\n--- STDERR ---\n')
                    log_file.write(result.stderr)

            if result.returncode != 0 or not os.path.exists(f1[1]):
                raise CException(self.__class__, 411, f'matthews_coef failed: {result.stderr}')

            # Parse results from XML file
            rv = {'results': []}

            from ccp4i2.core.CCP4Utils import openFileToEtree
            x_tree = openFileToEtree(fileName=f1[1])

            try:
                rv['cell_volume'] = float(x_tree.xpath('cell')[0].get('volume'))
            except Exception:
                pass

            x_result_list = x_tree.xpath('result')
            for x_result in x_result_list:
                rv['results'].append({
                    'nmol_in_asu': int(x_result.get('nmol_in_asu'))
                })
                for item in ['matth_coef', 'percent_solvent', 'prob_matth']:
                    rv['results'][-1][item] = float(x_result.get(item))

            return rv

        finally:
            # Clean up temporary files
            try:
                os.unlink(f1[1])
            except Exception:
                pass
            try:
                os.unlink(f2[1])
            except Exception:
                pass


class CMtzDataFile(CMtzDataFileStub):
    """
    An MTZ experimental data file

    Extends CMtzDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        super().__init__(file_path=file_path, parent=parent, name=name, **kwargs)
        # Note: fileContentClassName='CMtzData' is already set in CMtzDataFileStub decorator
        # Note: MIME type now comes from ccp4i2_static_data.py via get_file_type_from_class()
        self.set_qualifier('guiLabel', 'Experimental data')
        self.set_qualifier('toolTip', 'MTZ format reflection data file')
        self.set_qualifier('fileExtensions', ['mtz'])

    def fileExtensions(self):
        """
        Return file extension for MTZ files.

        All MTZ files (CMtzDataFile and derivatives) use the .mtz extension.

        Returns:
            list: ['mtz']
        """
        return ['mtz']

    def columnNames(self, asString=False, ifString=None):
        """
        Get column names from the file content, or expected names based on contentFlag.

        Legacy compatibility method that inspects fileContent to return column names.
        If fileContent is not available, falls back to CONTENT_SIGNATURE_LIST based
        on the contentFlag.

        Args:
            asString: If True, return comma-separated string. If False, return list.
            ifString: Legacy parameter name for asString (for backwards compatibility)

        Returns:
            list or str: Column names as list or comma-separated string

        Example:
            >>> mtz_file.columnNames(False)
            ['H', 'K', 'L', 'F', 'SIGF']
            >>> mtz_file.columnNames(True)
            'H,K,L,F,SIGF'
        """
        # Support legacy ifString parameter
        if ifString is not None:
            asString = ifString

        col_names = []

        # Try to get from fileContent.listOfColumns
        content = self.getFileContent()
        if content and hasattr(content, 'listOfColumns') and content.listOfColumns:
            # Extract column labels from CMtzColumn objects
            for col in content.listOfColumns:
                if hasattr(col, 'columnLabel'):
                    label = col.columnLabel
                    # Handle both CString with .value and plain strings
                    if hasattr(label, 'value'):
                        col_names.append(str(label.value))
                    else:
                        col_names.append(str(label))

        # Fall back to CONTENT_SIGNATURE_LIST if no fileContent columns
        if not col_names and hasattr(self.__class__, 'CONTENT_SIGNATURE_LIST'):
            content_flag = None
            if hasattr(self, 'contentFlag') and self.contentFlag:
                cf = self.contentFlag
                content_flag = cf.value if hasattr(cf, 'value') else int(cf)

            signature_list = self.__class__.CONTENT_SIGNATURE_LIST
            if content_flag is not None and 1 <= content_flag <= len(signature_list):
                # contentFlag is 1-indexed
                col_names = list(signature_list[content_flag - 1])

        if asString:
            return ','.join(col_names)
        else:
            return col_names


class CMtzDataset(CMtzDatasetStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CMtzDatasetStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CObsDataFile(CObsDataFileStub, CMiniMtzDataFile):
    """
    An MTZ experimental data file

    Inherits from:
    - CObsDataFileStub: Metadata and structure
    - CMiniMtzDataFile: Shared full-fat methods

    Extends CObsDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        super().__init__(file_path=file_path, parent=parent, name=name, **kwargs)
        # Note: MIME type now comes from ccp4i2_static_data.py via get_file_type_from_class()
        self.set_qualifier('guiLabel', 'Observed data')
        self.set_qualifier('toolTip', 'MTZ format observed data file')
        self.set_qualifier('fileExtensions', ['mtz'])

    def as_IPAIR(self, work_directory: Optional[Any] = None) -> str:
        """
        Convert this file to IPAIR format (Anomalous Intensities).

        IPAIR format: Iplus, SIGIplus, Iminus, SIGIminus

        Args:
            work_directory: Directory for output if input dir not writable

        Returns:
            Full path to converted file

        Raises:
            NotImplementedError: Conversion logic not yet implemented
        """
        output_path = self._get_conversion_output_path('IPAIR', work_directory=work_directory)

        # TODO: Implement actual conversion logic using gemmi
        # This will involve:
        # 1. Reading current MTZ file
        # 2. Converting data to IPAIR format
        # 3. Writing to output_path

        raise NotImplementedError(
            f"Conversion to IPAIR format not yet implemented. "
            f"Would output to: {output_path}"
        )

    def _ensure_container_child(self, container, name: str, child_class):
        """
        Ensure a container has a child of specified type.

        Helper method for conversion routines that need to initialize container children.

        Args:
            container: CContainer instance
            name: Name of the child attribute
            child_class: Class to instantiate if child doesn't exist

        Returns:
            The child instance (existing or newly created)
        """
        if not hasattr(container, name) or getattr(container, name) is None:
            child = child_class(name=name)
            setattr(container, name, child)
        return getattr(container, name)

    def as_FPAIR(self, work_directory: Optional[Any] = None) -> str:
        """
        Convert this file to FPAIR format (Anomalous Structure Factors).

        FPAIR format: Fplus, SIGFplus, Fminus, SIGFminus

        Uses ctruncate to convert IPAIR (anomalous intensities) to FPAIR
        (anomalous structure factor amplitudes) via French-Wilson conversion.

        This is a thin wrapper around ObsDataConverter.to_fpair().
        See core.conversions.obs_data_converter for implementation details.

        Args:
            work_directory: Directory for ctruncate working files

        Returns:
            Full path to converted file

        Raises:
            ValueError: If contentFlag cannot be determined or conversion not possible
            RuntimeError: If ctruncate plugin is not available or conversion fails
        """
        from ccp4i2.core.conversions import ObsDataConverter
        return ObsDataConverter.to_fpair(self, work_directory=work_directory)


    def as_IMEAN(self, work_directory: Optional[Any] = None) -> str:
        """
        Convert this file to IMEAN format (Mean Intensities).

        IMEAN format: I, SIGI

        Uses ctruncate to convert IPAIR (anomalous intensities) to IMEAN
        (mean intensities) by averaging I+ and I-.

        This is a thin wrapper around ObsDataConverter.to_imean().
        See core.conversions.obs_data_converter for implementation details.

        Args:
            work_directory: Directory for ctruncate working files

        Returns:
            Full path to converted file

        Raises:
            ValueError: If contentFlag cannot be determined or conversion not possible
            RuntimeError: If ctruncate plugin is not available or conversion fails
        """
        from ccp4i2.core.conversions import ObsDataConverter
        return ObsDataConverter.to_imean(self, work_directory=work_directory)

    def as_FPAIR(self, work_directory: Optional[Any] = None) -> str:
        """
        Convert this file to FPAIR format (Anomalous Structure Factors).

        FPAIR format: F+, SIGF+, F-, SIGF-

        Conversion path:
        - IPAIR → FPAIR: French-Wilson via servalcat fw

        This is a thin wrapper around ObsDataConverter.to_fpair().
        See core.conversions.obs_data_converter for implementation details.

        Args:
            work_directory: Directory for working files

        Returns:
            Full path to converted file

        Raises:
            ValueError: If contentFlag cannot be determined
            RuntimeError: If conversion fails
        """
        from ccp4i2.core.conversions import ObsDataConverter
        return ObsDataConverter.to_fpair(self, work_directory=work_directory)

    def as_FMEAN(self, work_directory: Optional[Any] = None) -> str:
        """
        Convert this file to FMEAN format (Mean Structure Factors).

        FMEAN format: F, SIGF

        Handles multiple input formats:
        - IPAIR → FMEAN: French-Wilson via servalcat fw
        - IMEAN → FMEAN: French-Wilson via servalcat fw
        - FPAIR → FMEAN: Inverse-variance weighted mean via gemmi

        This is a thin wrapper around ObsDataConverter.to_fmean().
        See core.conversions.obs_data_converter for implementation details.

        Args:
            work_directory: Directory for working files

        Returns:
            Full path to converted file

        Raises:
            ValueError: If contentFlag cannot be determined
            RuntimeError: If conversion fails
        """
        from ccp4i2.core.conversions import ObsDataConverter
        return ObsDataConverter.to_fmean(self, work_directory=work_directory)


class CPhaserSolDataFile(CPhaserSolDataFileStub):
    """
    Phaser solution data file (pickle format).

    Extends CPhaserSolDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CPhaserRFileDataFile(CPhaserRFileDataFileStub):
    """
    Phaser R-list data file (pickle format).

    Extends CPhaserRFileDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CProgramColumnGroup(CProgramColumnGroupStub):
    """
    A group of MTZ columns required for program input.

    This class maps between generic column names (e.g., 'Ip', 'SIGIp') and
    actual MTZ column labels (e.g., 'Iplus', 'SIGIplus'). It behaves like
    a dictionary that can be accessed via attributes.

    Example:
        >>> colgroup = CProgramColumnGroup()
        >>> colgroup.set({'Ip': 'Iplus', 'SIGIp': 'SIGIplus'})
        >>> print(colgroup.Ip)  # Returns 'Iplus'
        >>> print(colgroup.isSet())  # Returns True

    This provides compatibility with ccp4i2 wrapper code that expects
    attribute access to column mappings.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        """Initialize with an internal column mapping dict."""
        super().__init__(parent=parent, name=name, **kwargs)
        # Internal storage for column name mappings
        self._column_mapping = {}
        self._is_set = False

    def set(self, mapping: dict):
        """
        Set the column name mappings.

        Args:
            mapping: Dict mapping generic column names to actual MTZ column labels
                    e.g., {'Ip': 'Iplus', 'SIGIp': 'SIGIplus', 'Im': 'Iminus', 'SIGIm': 'SIGIminus'}
        """
        if isinstance(mapping, dict):
            self._column_mapping = mapping.copy()
            self._is_set = True
        else:
            # If it's a CData object, try to extract its value
            if hasattr(mapping, '__dict__') and '_column_mapping' in mapping.__dict__:
                self._column_mapping = mapping._column_mapping.copy()
                self._is_set = mapping._is_set
            else:
                raise TypeError(f"Expected dict, got {type(mapping)}")

    def isSet(self, field_name: str = None, allowUndefined: bool = False,
              allowDefault: bool = False, allSet: bool = True):
        """
        Check if column mappings have been set.

        Args:
            field_name: Ignored for CProgramColumnGroup (provided for API compatibility)
            allowUndefined: Ignored for CProgramColumnGroup (provided for API compatibility)
            allowDefault: Ignored for CProgramColumnGroup (provided for API compatibility)
            allSet: Ignored for CProgramColumnGroup (provided for API compatibility)

        Returns:
            bool: True if set() has been called with mappings AND mappings contain actual column data
        """
        # Check both the flag AND that we have actual column mappings
        # This prevents reporting True for empty/uninitialized column groups
        # Ignore 'name' key as it's just metadata, not a column mapping
        if not self._is_set:
            return False

        # Check if we have any column mappings besides 'name'
        actual_columns = {k: v for k, v in self._column_mapping.items() if k != 'name'}
        return len(actual_columns) > 0

    def __getattr__(self, name):
        """
        Allow attribute access to column mappings.

        Args:
            name: Generic column name (e.g., 'Ip', 'SIGIp')

        Returns:
            The mapped MTZ column label, or None if not found
        """
        # Avoid recursion for internal attributes
        if name.startswith('_'):
            raise AttributeError(f"'{type(self).__name__}' object has no attribute '{name}'")

        # Check if we have this column mapping
        if '_column_mapping' in self.__dict__ and name in self._column_mapping:
            return self._column_mapping[name]

        # For missing column names, return None (matches ccp4i2 behavior)
        # This allows wrapper code to check `if inp.ISIGI.I:` without errors
        return None

    def __setattr__(self, name, value):
        """
        Allow attribute assignment to update column mappings.

        Args:
            name: Generic column name
            value: Actual MTZ column label
        """
        # Internal attributes go through normal path
        if name.startswith('_'):
            super().__setattr__(name, value)
        else:
            # For column names, store in mapping
            if '_column_mapping' not in self.__dict__:
                self.__dict__['_column_mapping'] = {}
            if '_is_set' not in self.__dict__:
                self.__dict__['_is_set'] = False

            self._column_mapping[name] = value
            self._is_set = True

    def get(self, name, default=None):
        """
        Get a mapped column name with optional default.

        Args:
            name: Generic column name
            default: Value to return if name not found

        Returns:
            Mapped column label or default
        """
        return self._column_mapping.get(name, default)

    def keys(self):
        """Return all generic column names."""
        return self._column_mapping.keys()

    def values(self):
        """Return all mapped column labels."""
        return self._column_mapping.values()

    def items(self):
        """Return all (generic_name, mapped_label) pairs."""
        return self._column_mapping.items()


class CProgramColumnGroup0(CProgramColumnGroup0Stub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CProgramColumnGroup0Stub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CRefmacKeywordFile(CRefmacKeywordFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CRefmacKeywordFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CReindexOperator(CReindexOperatorStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CReindexOperatorStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CResolutionRange(CFloatRange, CResolutionRangeStub):
    """
    Resolution range for crystallographic data.

    Resolution range convention:
    - low: larger d-spacing number (lower resolution, e.g., 50.0 Å)
    - high: smaller d-spacing number (higher resolution, e.g., 2.0 Å)

    Properties:
    - .start and .end: Actual CFloat attributes (inherited from CFloatRange)
    - .low and .high: Properties that map to .start and .end respectively

    The .low/.high properties return CFloat objects (not primitive values)
    so that code can call .isSet() on them.

    When only .high (or .end) is set, aimless will use data from infinite
    resolution down to the high resolution cutoff. This is the desired behavior
    for workflows like SubstituteLigand.

    NOTE: The __init__ and _smart_assign_from_cdata methods are inherited from
    CFloatRange, which ensures that .start and .end are NOT_SET by default and
    that smart assignment only copies explicitly set fields.
    """

    @property
    def low(self):
        """
        Low resolution limit (larger d-spacing, e.g., 50.0 Å).

        Returns CFloat object (not primitive value) so .isSet() can be called.
        Maps to .start attribute.
        """
        return self.start

    @low.setter
    def low(self, value):
        """Set low resolution limit via .low property."""
        if isinstance(value, CFloat):
            self.start = value
        else:
            self.start.value = value

    @property
    def high(self):
        """
        High resolution limit (smaller d-spacing, e.g., 2.0 Å).

        Returns CFloat object (not primitive value) so .isSet() can be called.
        Maps to .end attribute.
        """
        return self.end

    @high.setter
    def high(self, value):
        """Set high resolution limit via .high property."""
        if isinstance(value, CFloat):
            self.end = value
        else:
            self.end.value = value

    # Legacy methods for backward compatibility
    def setHigh(self, value: Optional[float]) -> None:
        """Set the high resolution limit (smaller d-spacing value)."""
        self.high = value

    def setLow(self, value: Optional[float]) -> None:
        """Set the low resolution limit (larger d-spacing value)."""
        self.low = value



class CRunBatchRange(CRunBatchRangeStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CRunBatchRangeStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CRunBatchRangeList(CRunBatchRangeListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CRunBatchRangeListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CShelxFADataFile(CShelxFADataFileStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CShelxFADataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CShelxLabel(CShelxLabelStub):
    """
    A string
    
    Extends CShelxLabelStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSpaceGroup(CSpaceGroupStub):
    """
    A string holding the space group

    Extends CSpaceGroupStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def fix(self, value: str) -> str:
        """
        Normalize a space group string.

        Takes a space group string (e.g., from program output) and returns
        a normalized version suitable for storage.

        Args:
            value: Space group string to normalize

        Returns:
            Normalized space group string
        """
        if value is None:
            return ""
        # Strip whitespace and return
        # More sophisticated normalization (e.g., P 21 21 21 -> P212121)
        # can be added here if needed
        return str(value).strip()


class CSpaceGroupCell(CSpaceGroupCellStub):
    """
    Cell space group and parameters
    
    Extends CSpaceGroupCellStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CUnmergedDataContent(CUnmergedDataContentStub):
    """
    Base class for classes holding file contents

    Extends CUnmergedDataContentStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)
        # Initialize attributes not in stub but needed by legacy code
        from ccp4i2.core.base_object.fundamental_types import CList
        if not hasattr(self, 'listOfColumns') or self.listOfColumns is None:
            self.listOfColumns = CList(name='listOfColumns', parent=self)
        if not hasattr(self, 'datasets') or self.datasets is None:
            self.datasets = CList(name='datasets', parent=self)
        if not hasattr(self, 'crystalNames') or self.crystalNames is None:
            self.crystalNames = CList(name='crystalNames', parent=self)
        if not hasattr(self, 'wavelengths') or self.wavelengths is None:
            self.wavelengths = CList(name='wavelengths', parent=self)

    def to_dict(self):
        """
        Convert CUnmergedDataContent to a dictionary for serialization.

        Includes all standard children plus type-specific attributes like
        datasets, wavelengths, crystalNames, listOfColumns, datasetCells, and batchs.
        """
        from ccp4i2.lib.utils.parameters.value_dict import value_dict_for_object

        # Start with base serialization from children
        result = {}
        try:
            for child in self.children():
                if hasattr(child, 'objectName') and callable(child.objectName):
                    name = child.objectName()
                elif hasattr(child, 'object_name') and callable(child.object_name):
                    name = child.object_name()
                elif hasattr(child, 'name'):
                    name = child.name
                else:
                    continue
                result[name] = value_dict_for_object(child)
        except Exception:
            pass

        # Add type-specific attributes that may be plain Python lists or strings
        type_specific_attrs = [
            'datasets', 'wavelengths', 'crystalNames', 'datasetCells', 'listOfColumns', 'batchs'
        ]
        for attr_name in type_specific_attrs:
            if attr_name not in result and hasattr(self, attr_name):
                try:
                    attr_value = getattr(self, attr_name)
                    if attr_value is not None:
                        if isinstance(attr_value, str) and attr_value:
                            # Handle string attributes like batchs
                            result[attr_name] = attr_value
                        elif isinstance(attr_value, list) and len(attr_value) > 0:
                            result[attr_name] = value_dict_for_object(attr_value)
                        elif hasattr(attr_value, 'value') and attr_value.value:
                            result[attr_name] = value_dict_for_object(attr_value)
                        elif hasattr(attr_value, '__iter__'):
                            # Handle CList or other iterables
                            converted = value_dict_for_object(attr_value)
                            if converted:
                                result[attr_name] = converted
                except Exception:
                    pass

        return result

    def loadFile(self, file_path: str = None):
        """
        Load unmerged reflection file using gemmi library.

        Supports multiple formats:
        - MTZ (using gemmi.read_mtz_file)
        - mmCIF (using gemmi.cif.read)
        - Other formats detected by extension

        Args:
            file_path: Optional path to reflection file. If None, gets path from parent CDataFile.

        Returns:
            CErrorReport with any errors encountered

        Example:
            # Load from parent file's path
            >>> unmerged_file = CUnmergedDataFile()
            >>> unmerged_file.setFullPath('/path/to/unmerged.mtz')
            >>> unmerged_file.fileContent.loadFile()

            # Load from explicit path (legacy pattern)
            >>> unmerged = CUnmergedDataContent()
            >>> error = unmerged.loadFile('/path/to/unmerged.mtz')
            >>> if error.count() == 0:
            ...     print(f"Format: {unmerged.format.value}")
            ...     print(f"Merged: {unmerged.merged.value}")
        """
        from ccp4i2.core.base_object.error_reporting import CErrorReport
        from pathlib import Path

        error = CErrorReport()

        # If no path provided, get from parent CDataFile
        if file_path is None:
            parent = self.get_parent()
            if parent is not None and hasattr(parent, 'getFullPath'):
                file_path = parent.getFullPath()

        # Validate file path
        if not file_path:
            return error

        path_obj = Path(file_path)
        if not path_obj.exists() or not path_obj.is_file():
            error.append(
                klass=self.__class__.__name__,
                code=101,
                details=f"File does not exist or is not a file: '{file_path}'",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        try:
            import gemmi
        except ImportError as e:
            error.append(
                klass=self.__class__.__name__,
                code=102,
                details=f"Failed to import gemmi library: {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return error

        # Determine format from extension
        suffix = path_obj.suffix.lower()

        try:
            # Handle MTZ files
            if suffix == '.mtz':
                self._load_mtz_file(file_path, gemmi, error)

            # Handle mmCIF files
            elif suffix in ['.cif', '.mmcif', '.ent']:
                self._load_mmcif_file(file_path, gemmi, error)

            # Handle Scalepack format (.sca, .hkl)
            elif suffix in ['.sca', '.hkl']:
                self._load_scalepack_file(file_path, error)

            # Handle XDS files (INTEGRATE.HKL, XDS_ASCII.HKL)
            elif 'INTEGRATE' in path_obj.name or 'XDS_ASCII' in path_obj.name or suffix == '.hkl':
                # Try XDS format first
                try:
                    self._load_xds_file(file_path, gemmi, error)
                except:
                    # Fall back to Scalepack
                    self._load_scalepack_file(file_path, error)

            elif suffix == '.shelx':
                self.format = 'shelx'
                self.merged = 'unk'
                self.knowncell = False
                self.knownwavelength = False
                object.__setattr__(self, '_file_path', file_path)

            else:
                # Try MTZ as default for unknown extensions
                try:
                    self._load_mtz_file(file_path, gemmi, error)
                except:
                    self.format = 'unk'
                    error.append(
                        klass=self.__class__.__name__,
                        code=103,
                        details=f"Unknown file format for: '{file_path}'",
                        name=self.object_name() if hasattr(self, 'object_name') else ''
                    )

        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=104,
                details=f"Error reading file '{file_path}': {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )

        return error

    def _load_mtz_file(self, file_path: str, gemmi, error):
        """Load MTZ file and extract metadata."""
        mtz = gemmi.read_mtz_file(str(file_path))

        # Set format
        self.format = 'mtz'

        # Determine if merged (unmerged MTZ have batch columns)
        has_batch = any(col.label == 'BATCH' for col in mtz.columns)
        self.merged = 'unmerged' if has_batch else 'merged'

        # Extract cell using smart assignment
        if hasattr(self, 'cell') and self.cell is not None:
            self.cell.a = mtz.cell.a
            self.cell.b = mtz.cell.b
            self.cell.c = mtz.cell.c
            self.cell.alpha = mtz.cell.alpha
            self.cell.beta = mtz.cell.beta
            self.cell.gamma = mtz.cell.gamma

        # Extract space group
        if hasattr(self, 'spaceGroup') and self.spaceGroup is not None:
            self.spaceGroup = mtz.spacegroup.hm

        # Extract resolution range
        mtz.update_reso()
        if hasattr(self, 'lowRes') and self.lowRes is not None:
            self.lowRes = mtz.resolution_low()
        if hasattr(self, 'highRes') and self.highRes is not None:
            self.highRes = mtz.resolution_high()

        # Extract dataset information
        if len(mtz.datasets) > 0:
            # Use first non-HKL_base dataset
            for ds in mtz.datasets:
                if ds.dataset_name != 'HKL_base':
                    if hasattr(self, 'datasetName') and self.datasetName is not None:
                        self.datasetName = ds.dataset_name
                    if hasattr(self, 'crystalName') and self.crystalName is not None:
                        self.crystalName = ds.crystal_name
                    if hasattr(self, 'wavelength') and self.wavelength is not None:
                        self.wavelength = ds.wavelength
                    break

        # Number of datasets
        if hasattr(self, 'numberofdatasets') and self.numberofdatasets is not None:
            self.numberofdatasets = len(mtz.datasets)

        # Number of lattices (count of batches) and batch numbers
        if has_batch:
            if hasattr(self, 'numberLattices') and self.numberLattices is not None:
                self.numberLattices = len(mtz.batches)
            # Populate batchs with compact range notation (e.g., "1-100, 302-305")
            if hasattr(self, 'batchs') and mtz.batches:
                batch_numbers = sorted([batch.number for batch in mtz.batches])
                self.batchs = _compact_batch_ranges(batch_numbers)

        # Extract column information as CMtzColumn objects
        # Skip H-type columns (Miller indices H, K, L)
        if hasattr(self, 'listOfColumns') and self.listOfColumns is not None:
            # Clear existing columns to prevent duplication on re-load
            self.listOfColumns.clear()
            # Create CMtzColumn objects with columnLabel and columnType
            mtz_columns = []
            for col in mtz.columns:
                # Skip Miller indices (type 'H')
                if col.type == 'H':
                    continue

                mtz_col = CMtzColumn(name=col.label)
                mtz_col.columnLabel = col.label
                mtz_col.columnType = col.type
                # Get dataset name and use dataset ID as groupIndex
                if col.dataset:
                    mtz_col.dataset = col.dataset.dataset_name
                    mtz_col.groupIndex = col.dataset.id
                mtz_columns.append(mtz_col)
            # Populate the CList using append (direct assignment to CList is silently ignored)
            for mtz_col in mtz_columns:
                self.listOfColumns.append(mtz_col)

        # Extract datasets, crystal names, and wavelengths lists
        if hasattr(self, 'datasets') and self.datasets is not None:
            dataset_names = [ds.dataset_name for ds in mtz.datasets]
            self.datasets = dataset_names
        if hasattr(self, 'crystalNames') and self.crystalNames is not None:
            crystal_names = [ds.crystal_name for ds in mtz.datasets]
            self.crystalNames = crystal_names
        if hasattr(self, 'wavelengths') and self.wavelengths is not None:
            wavelength_list = [ds.wavelength for ds in mtz.datasets]
            self.wavelengths = wavelength_list

        # Cell and wavelength are known for MTZ
        self.knowncell = True
        self.knownwavelength = True

        # Store gemmi Mtz object
        object.__setattr__(self, '_gemmi_mtz', mtz)

    def _load_mmcif_file(self, file_path: str, gemmi, error):
        """Load mmCIF file and extract metadata."""
        cif_doc = gemmi.cif.read(str(file_path))

        if len(cif_doc) == 0:
            error.append(
                klass=self.__class__.__name__,
                code=105,
                details=f"mmCIF file contains no data blocks: '{file_path}'",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )
            return

        block = cif_doc[0]

        # Set format
        self.format = 'mmcif'

        # mmCIF reflection files are typically merged
        self.merged = 'merged'

        # Extract cell
        if hasattr(self, 'cell') and self.cell is not None:
            try:
                self.cell.a = float(block.find_value('_cell.length_a'))
                self.cell.b = float(block.find_value('_cell.length_b'))
                self.cell.c = float(block.find_value('_cell.length_c'))
                self.cell.alpha = float(block.find_value('_cell.angle_alpha'))
                self.cell.beta = float(block.find_value('_cell.angle_beta'))
                self.cell.gamma = float(block.find_value('_cell.angle_gamma'))
                self.knowncell = True
            except (ValueError, RuntimeError):
                self.knowncell = False

        # Extract space group
        if hasattr(self, 'spaceGroup') and self.spaceGroup is not None:
            try:
                sg_name = block.find_value('_symmetry.space_group_name_H-M')
                if sg_name:
                    self.spaceGroup = sg_name.strip('"').strip("'").strip()
            except RuntimeError:
                pass

        # Extract wavelength
        if hasattr(self, 'wavelength') and self.wavelength is not None:
            try:
                wavelength_val = block.find_value('_diffrn_radiation_wavelength.wavelength')
                if wavelength_val:
                    self.wavelength = float(wavelength_val)
                    self.knownwavelength = True
                else:
                    self.knownwavelength = False
            except (ValueError, RuntimeError):
                self.knownwavelength = False

        # Store gemmi CIF document
        object.__setattr__(self, '_gemmi_cif_doc', cif_doc)

    def _load_xds_file(self, file_path: str, gemmi, error):
        """Load XDS file (INTEGRATE.HKL or XDS_ASCII.HKL) using gemmi."""
        from pathlib import Path

        xds = gemmi.read_xds_ascii(str(file_path))

        # Set format
        if 'INTEGRATE' in Path(file_path).name:
            self.format = 'xds'
        else:
            self.format = 'xds'

        # XDS files are always unmerged
        self.merged = 'unmerged'

        # Extract cell using smart assignment
        if hasattr(self, 'cell') and self.cell is not None:
            cell_params = xds.cell_constants
            self.cell.a = cell_params[0]
            self.cell.b = cell_params[1]
            self.cell.c = cell_params[2]
            self.cell.alpha = cell_params[3]
            self.cell.beta = cell_params[4]
            self.cell.gamma = cell_params[5]
            self.knowncell = True
        else:
            self.knowncell = False

        # Extract space group from number
        if hasattr(self, 'spaceGroup') and self.spaceGroup is not None:
            sg_num = xds.spacegroup_number
            # Convert to Hermann-Mauguin symbol if possible
            try:
                from gemmi import SpaceGroup
                sg = SpaceGroup(sg_num)
                self.spaceGroup = sg.hm
            except:
                self.spaceGroup = f"Space group {sg_num}"

        # Extract wavelength
        if hasattr(self, 'wavelength') and self.wavelength is not None:
            self.wavelength = xds.wavelength
            self.knownwavelength = True
        else:
            self.knownwavelength = False

        # Number of reflections
        data_size = xds.data_size
        if data_size > 0:
            # Store data size info
            object.__setattr__(self, '_xds_data_size', data_size)

        # Store gemmi XDS object
        object.__setattr__(self, '_gemmi_xds', xds)

    def _load_scalepack_file(self, file_path: str, error):
        """Load Scalepack (.sca) file header to extract metadata.

        Scalepack has two formats:

        Merged format:
          Line 1: Version (1 = merged)
          Line 2: -987 (magic number)
          Line 3: a b c alpha beta gamma space_group

        Unmerged format:
          Line 1: nsyms space_group (e.g., "8 C2221")
          Line 2+: Symmetry matrices (nsyms * 3 lines)
          No cell parameters in file
        """
        try:
            with open(file_path, 'r') as f:
                line1 = f.readline().strip()
                line2 = f.readline().strip()
                line3 = f.readline().strip()

            # Set format
            self.format = 'sca'

            # Detect file format by checking line 2
            if line2 == '-987':
                # MERGED FORMAT: line1 = version, line2 = -987, line3 = cell
                try:
                    version = int(line1)
                    if version == 1:
                        self.merged = 'merged'
                    elif version in [2, 3]:
                        self.merged = 'unmerged'
                    else:
                        self.merged = 'unk'
                except ValueError:
                    self.merged = 'unk'

                # Parse cell parameters from line 3
                parts = line3.split()
                if len(parts) >= 6:
                    try:
                        self.cell.a = float(parts[0])
                        self.cell.b = float(parts[1])
                        self.cell.c = float(parts[2])
                        self.cell.alpha = float(parts[3])
                        self.cell.beta = float(parts[4])
                        self.cell.gamma = float(parts[5])
                        self.knowncell = True

                        # Space group is after cell parameters
                        if len(parts) > 6:
                            sg_text = ' '.join(parts[6:])
                            self.spaceGroup = sg_text.strip()
                    except (ValueError, IndexError):
                        self.knowncell = False
                else:
                    self.knowncell = False

                # Merged format doesn't have wavelength
                self.knownwavelength = False

            else:
                # UNMERGED FORMAT: line1 = "nsyms spacegroup"
                self.merged = 'unmerged'

                # Parse space group from line 1
                parts = line1.split(None, 1)  # Split on first whitespace
                if len(parts) >= 2:
                    try:
                        # parts[0] is nsyms (number of symmetry operations)
                        spacegroup = parts[1]
                        self.spaceGroup = spacegroup
                    except (ValueError, IndexError):
                        pass

                # Unmerged Scalepack files do not contain cell or wavelength
                self.knowncell = False
                self.knownwavelength = False

            # Store file path for potential future full parsing
            object.__setattr__(self, '_file_path', file_path)

        except Exception as e:
            error.append(
                klass=self.__class__.__name__,
                code=107,
                details=f"Error parsing Scalepack file '{file_path}': {e}",
                name=self.object_name() if hasattr(self, 'object_name') else ''
            )

    def getColumnGroups(self):
        """
        Build list of column groups from MTZ column info (for mini-mtz's).

        This method is identical to CMtzData.getColumnGroups() since CUnmergedDataContent
        also has listOfColumns populated with CMtzColumn objects.

        Returns:
            list: List of CColumnGroup objects with columnList, columnGroupType, contentFlag

        Example:
            >>> groups = unmerged.getColumnGroups()
            >>> for group in groups:
            ...     print(f"{group.columnGroupType}: {group.columnList}")
        """
        groupIndex = 0  # MapCoeffs broken KJS
        groupList = []

        # Get listOfColumns - handle both CList and plain list
        if not hasattr(self, 'listOfColumns') or self.listOfColumns is None:
            return groupList

        columns = self.listOfColumns.value if hasattr(self.listOfColumns, 'value') else self.listOfColumns
        if not columns:
            return groupList

        # Sort listOfColumns into a list of grouped columns
        for ii in range(len(columns)):
            fileColumn = columns[ii]
            col_group_idx = fileColumn.groupIndex.value if hasattr(fileColumn.groupIndex, 'value') else fileColumn.groupIndex

            if col_group_idx != groupIndex:
                col_dataset = fileColumn.dataset.value if hasattr(fileColumn.dataset, 'value') else fileColumn.dataset
                groupList.append(CColumnGroup(dataset=str(col_dataset)))
                groupIndex = int(col_group_idx)

            # Ensure we have at least one group (defensive check for first column)
            if not groupList:
                col_dataset = fileColumn.dataset.value if hasattr(fileColumn.dataset, 'value') else fileColumn.dataset
                groupList.append(CColumnGroup(dataset=str(col_dataset)))
                groupIndex = int(col_group_idx)

            col_label = fileColumn.columnLabel.value if hasattr(fileColumn.columnLabel, 'value') else fileColumn.columnLabel
            col_type = fileColumn.columnType.value if hasattr(fileColumn.columnType, 'value') else fileColumn.columnType
            col_dataset = fileColumn.dataset.value if hasattr(fileColumn.dataset, 'value') else fileColumn.dataset
            col_group_idx = fileColumn.groupIndex.value if hasattr(fileColumn.groupIndex, 'value') else fileColumn.groupIndex

            groupList[-1].columnList.append(CMtzColumn(
                columnLabel=col_label,
                columnType=col_type,
                groupIndex=col_group_idx,
                dataset=col_dataset
            ))

        # Loop over and catch map-coefs (F-Phi pairs)
        for i, col1 in enumerate(columns[:-1]):
            col2 = columns[i+1]
            col1_type = col1.columnType.value if hasattr(col1.columnType, 'value') else col1.columnType
            col2_type = col2.columnType.value if hasattr(col2.columnType, 'value') else col2.columnType
            MapT = (col1_type == 'F' and col2_type == 'P')  # Map (FP): F-Phi / Phases (PW): Phi-FOM

            if MapT:
                col1_dataset = col1.dataset.value if hasattr(col1.dataset, 'value') else col1.dataset
                groupList.append(CColumnGroup(dataset=str(col1_dataset)))

                col1_label = col1.columnLabel.value if hasattr(col1.columnLabel, 'value') else col1.columnLabel
                col1_group = col1.groupIndex.value if hasattr(col1.groupIndex, 'value') else col1.groupIndex
                col2_label = col2.columnLabel.value if hasattr(col2.columnLabel, 'value') else col2.columnLabel
                col2_group = col2.groupIndex.value if hasattr(col2.groupIndex, 'value') else col2.groupIndex

                groupList[-1].columnList.append(CMtzColumn(
                    columnLabel=col1_label,
                    columnType=col1_type,
                    groupIndex=col1_group,
                    dataset=col1_dataset
                ))
                groupList[-1].columnList.append(CMtzColumn(
                    columnLabel=col2_label,
                    columnType=col2_type,
                    groupIndex=col2_group,
                    dataset=col1_dataset
                ))

        # Assign the columnGroupType and contentFlag based on column signatures
        for group in groupList:
            signature = ''
            labels = []
            col_list = group.columnList.value if hasattr(group.columnList, 'value') else group.columnList
            for col in col_list:
                col_type = col.columnType.value if hasattr(col.columnType, 'value') else col.columnType
                col_label = col.columnLabel.value if hasattr(col.columnLabel, 'value') else col.columnLabel

                # Skip H-type columns (Miller indices H, K, L) when building signature
                if col_type != 'H':
                    signature = signature + str(col_type)
                    labels.append(str(col_label))

            # Check against known column signatures for each data type
            from ccp4i2.core.base_object.class_metadata import get_class_metadata_by_type
            for cls, label in [[CObsDataFile, 'Obs'], [CPhsDataFile, 'Phs'],
                              [CMapCoeffsDataFile, 'MapCoeffs'], [CFreeRDataFile, 'FreeR']]:
                meta = get_class_metadata_by_type(cls)
                correct_columns = meta.qualifiers.get('correctColumns') if meta and meta.qualifiers else None
                if correct_columns and signature in correct_columns:
                    if label == 'FreeR':
                        # Is this really FreeR? check column label
                        if (len(labels) == 1) and ('free' in str(labels[0]).lower()):
                            print(f"Column recognised as FreeR, label: {str(labels[0])}")
                            group.columnGroupType = label
                            group.contentFlag = correct_columns.index(signature) + 1
                    else:
                        group.columnGroupType = label
                        group.contentFlag = correct_columns.index(signature) + 1

            # Handle special cases if columnGroupType not set
            if not group.columnGroupType.isSet():
                if signature == 'FPW':
                    group.columnGroupType = 'Phs'
                    group.contentFlag = 2
                    col_list.pop(0)
                elif signature == 'PWF':
                    group.columnGroupType = 'Phs'
                    group.contentFlag = 2
                    col_list.pop(2)
                elif signature == 'FQDQ':
                    group.columnGroupType = 'Obs'
                    col_list.pop(3)
                    col_list.pop(2)
                    group.contentFlag = 4

        return groupList

    def clipperSameCell(self, other_content, tolerance=None):
        """
        Compare unit cells using Clipper's reciprocal space algorithm.

        Two cells disagree if the difference in their orthogonalisation matrices
        is sufficient to map a reflection from one cell onto a different reflection
        in the other cell at the given tolerance (resolution in Angstroms).

        This implements the Clipper Cell::equals() algorithm which considers
        reciprocal space vectors and finds the resolution at which reflections
        would be mis-indexed by more than 0.5 reciprocal lattice units.

        Args:
            other_content: Another CUnmergedDataContent or CMtzData instance
            tolerance: Resolution tolerance in Angstroms (default 1.0)
                      Cells are compatible if mis-indexing doesn't occur
                      until this resolution or higher.

        Returns:
            dict with keys:
                'validity': bool - True if cells are compatible at tolerance
                'tolerance': float - the resolution tolerance value
                'difference': float - resolution where mis-indexing starts (Å)
                'maximumResolution1': float - max resolution for cell1
                'maximumResolution2': float - max resolution for cell2
        """
        import numpy as np

        # Use default tolerance if None (Clipper default is 1.0 Angstrom)
        if tolerance is None:
            tolerance = 1.0

        # Get cell parameters
        cell1 = self.cell
        cell2 = other_content.cell

        if not cell1 or not cell2:
            return {
                'validity': False,
                'tolerance': tolerance,
                'difference': float('inf'),
                'maximumResolution1': 0.0,
                'maximumResolution2': 0.0
            }

        # Build orthogonalization matrices for both cells
        # Orthogonalization matrix converts fractional to Cartesian coordinates
        def build_orth_matrix(cell):
            """Build orthogonalization matrix from cell parameters."""
            import math
            a, b, c = float(cell.a), float(cell.b), float(cell.c)
            alpha = math.radians(float(cell.alpha))
            beta = math.radians(float(cell.beta))
            gamma = math.radians(float(cell.gamma))

            # Volume calculation
            cos_alpha = math.cos(alpha)
            cos_beta = math.cos(beta)
            cos_gamma = math.cos(gamma)
            sin_gamma = math.sin(gamma)

            volume = a * b * c * math.sqrt(
                1 - cos_alpha**2 - cos_beta**2 - cos_gamma**2
                + 2 * cos_alpha * cos_beta * cos_gamma
            )

            # Orthogonalization matrix (fractional to Cartesian)
            orth = np.array([
                [a, b * cos_gamma, c * cos_beta],
                [0, b * sin_gamma, c * (cos_alpha - cos_beta * cos_gamma) / sin_gamma],
                [0, 0, volume / (a * b * sin_gamma)]
            ])
            return orth

        orth1 = build_orth_matrix(cell1)
        orth2 = build_orth_matrix(cell2)

        # Calculate the difference matrix
        diff_matrix = orth1 - orth2

        # Find the resolution at which a reflection would be mis-indexed
        # by 0.5 reciprocal lattice units or more
        #
        # We test reflections along the crystallographic axes (h00, 0k0, 00l)
        # at increasing resolution until we find where the difference exceeds
        # 0.5 reciprocal lattice units

        max_resolution = 1000.0  # Start with very high resolution (small d)
        axes = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  # Test along a*, b*, c*

        for axis in axes:
            h, k, l = axis
            # Compute reciprocal space vector for this reflection
            hkl = np.array([h, k, l])

            # Transform to Cartesian reciprocal space
            # (reciprocal of orth matrix)
            orth1_inv = np.linalg.inv(orth1).T
            orth2_inv = np.linalg.inv(orth2).T

            # Reciprocal space vectors
            s1 = orth1_inv @ hkl
            s2 = orth2_inv @ hkl

            # Difference in reciprocal space (in Å⁻¹)
            diff = np.linalg.norm(s1 - s2)

            # Resolution of this reflection
            d_spacing = 1.0 / np.linalg.norm(s1)

            # At what resolution would this axis cause 0.5 r.l.u. mis-indexing?
            # diff * d = 0.5 / |hkl|
            # So critical resolution = 0.5 / (diff * |hkl|)
            if diff > 1e-10:  # Avoid division by zero
                critical_res = 0.5 / (diff * np.linalg.norm(hkl))
                max_resolution = min(max_resolution, critical_res)

        # Calculate maximum resolutions based on cell dimensions
        a1, b1, c1 = float(cell1.a), float(cell1.b), float(cell1.c)
        a2, b2, c2 = float(cell2.a), float(cell2.b), float(cell2.c)
        max_res1 = min(a1, b1, c1) / 2.0
        max_res2 = min(a2, b2, c2) / 2.0

        # Cells are compatible if mis-indexing doesn't occur until
        # resolution better than (smaller than) tolerance
        validity = max_resolution <= tolerance

        return {
            'validity': validity,
            'tolerance': tolerance,
            'difference': max_resolution,  # Resolution where mis-indexing starts
            'maximumResolution1': max_res1,
            'maximumResolution2': max_res2
        }


class CUnmergedDataFile(CUnmergedDataFileStub):
    """
    Handle MTZ, XDS and scalepack files. Allow wildcard filename

    Extends CUnmergedDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, file_path: str = None, parent=None, name=None, **kwargs):
        super().__init__(file_path=file_path, parent=parent, name=name, **kwargs)
        # Note: fileContentClassName='CUnmergedDataContent' is already set in CUnmergedDataFileStub decorator


class CUnmergedDataFileList(CUnmergedDataFileListStub):
    """
    A list with all items of one CData sub-class

    Extends CUnmergedDataFileListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    def __init__(self, parent=None, name=None, **kwargs):
        super().__init__(parent=parent, name=name, **kwargs)
        # Set the subItem qualifier to specify CUnmergedDataFile as item type
        # This allows aimless_pipe to use list items as CUnmergedDataFile objects
        self.set_qualifier('subItem', {'class': CUnmergedDataFile, 'qualifiers': {}})


class CUnmergedMtzDataFile(CUnmergedMtzDataFileStub, CMtzDataFile):
    """
    
    Inherits from:
    - CUnmergedMtzDataFileStub: Metadata and structure
    - CMtzDataFile: Shared full-fat methods
    An MTZ experimental data file
    
    Extends CUnmergedMtzDataFileStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CWavelength(CWavelengthStub):
    """
    Wavelength in Angstrom
    
    Extends CWavelengthStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CXia2ImageSelection(CXia2ImageSelectionStub):
    """
    QObject(self, parent: typing.Optional[PySide2.QtCore.QObject] = None) -> None
    
    Extends CXia2ImageSelectionStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CXia2ImageSelectionList(CXia2ImageSelectionListStub):
    """
    A list with all items of one CData sub-class
    
    Extends CXia2ImageSelectionListStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass

