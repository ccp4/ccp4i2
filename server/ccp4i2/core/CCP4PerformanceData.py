"""
CCP4PerformanceData.py --- CData classes.

These were once two classes each: a generated stub carrying the data
model, and an implementation carrying the methods. The generator is no
longer run and the split cost more than it saved --- the two halves
interleaved in the MRO, so an implementation could drop out of its own
subclass's ancestry and `isinstance` would say no to a file that plainly
was one. They are one class now.
"""

from typing import Optional, Any
from typing import TYPE_CHECKING, Optional, Any
from ccp4i2.core.base_object.class_metadata import cdata_class, attribute, AttributeType, content
from ccp4i2.core.base_object.base_classes import CData
from ccp4i2.core.base_object.fundamental_types import CFloat, CInt, CString
from ccp4i2.core.CCP4XtalData import CSpaceGroup



# Define CPerformanceIndicator FIRST since other classes inherit from it
class CPerformanceIndicator(CData):


    class Meta:
        error_codes = {
            "301": {
                "description": "Data value not set"
            },
            "302": {
                "description": "Performance indicator value difference greater than tolereance"
            },
            "303": {
                "description": "Performance indicator value different"
            },
            "304": {
                "description": "Performance indicator value difference greater than tolereance - but improved",
                "severity": 2
            },
            "305": {
                "description": "Performance indicator not used",
                "severity": 0
            }
        }
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
    value = content("CFloat", min=0.0)
    annotation = content("CString")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPerformanceIndicator.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Extends CPerformanceIndicator with implementation-specific methods.

    Overrides isSet() so that a performance indicator is considered "set"
    when any of its KPI children (beyond the inherited 'value' and
    'annotation' fields) have been explicitly set.  This ensures
    exclude_unset serialization preserves PERFORMANCE elements that
    contain real metrics like highResLimit, rMeas, RFactor, etc.
    """

    # Fields inherited from CPerformanceIndicator that are structural,
    # not KPI values — skip these when deciding if the indicator is "set".
    _STRUCTURAL_FIELDS = frozenset(('value', 'annotation'))

    def isSet(self, field_name=None, allowUndefined=False,
              allowDefault=False, allSet=False):
        # When checking the 'value' field or the object as a whole,
        # a performance indicator is "set" if any KPI child is set.
        # (getEtree calls isSet('value', ...) which would otherwise
        # only check the inherited CFloat 'value' attribute.)
        if field_name is None or field_name == 'value':
            from ccp4i2.core.base_object.cdata import CData
            for child in self.children():
                if not isinstance(child, CData):
                    continue
                name = child.objectName() if hasattr(child, 'objectName') else None
                if name in self._STRUCTURAL_FIELDS:
                    continue
                if child.isSet(allowDefault=allowDefault,
                               allowUndefined=allowUndefined):
                    return True
            return False

        # For any other specific field, delegate to base
        return super().isSet(field_name=field_name,
                             allowUndefined=allowUndefined,
                             allowDefault=allowDefault, allSet=allSet)

    def _value_explicitly_set(self):
        """Check if the inherited 'value' CFloat was explicitly set (not just default 0.0)."""
        return super().isSet(field_name='value', allowDefault=False)

    def getEtree(self, name=None, excludeUnset=False, allSet=False):
        """Override to avoid writing default 0.0 as element text for multi-metric indicators."""
        elem = super().getEtree(name=name, excludeUnset=excludeUnset, allSet=allSet)
        # If 'value' was never explicitly set, remove the spurious 0.0 text
        if not self._value_explicitly_set() and elem.text is not None:
            elem.text = None
        return elem


class CDataReductionPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['spaceGroup', 'highResLimit', 'rMeas']
    spaceGroup = content("CSpaceGroup")
    highResLimit = content("CFloat", min=0.0)
    rMeas = content("CFloat")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDataReductionPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CDataReductionPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CDataReductionPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CPairefPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['cutoff']
    cutoff = content("CFloat", min=0.0)

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPairefPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CPairefPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CPairefPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDataReductionCCPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['spaceGroup', 'highResLimit', 'ccHalf']
    spaceGroup = content("CSpaceGroup")
    highResLimit = content("CFloat", min=0.0)
    ccHalf = content("CFloat")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDataReductionCCPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CDataReductionCCPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CDataReductionCCPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CServalcatPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['RFactor', 'RFree', 'R', 'R1Factor', 'R1Free', 'R1', 'FSCaverage', 'annotation']
    RFactor = content("CFloat", min=0.0)
    RFree = content("CFloat", min=0.0)
    R = content("CFloat", min=0.0)
    R1Factor = content("CFloat", min=0.0)
    R1Free = content("CFloat", min=0.0)
    R1 = content("CFloat", min=0.0)
    FSCaverage = content("CFloat", min=-1.0)
    CCFwork_avg = content("CFloat", min=-1.0)
    CCFfree_avg = content("CFloat", min=-1.0)
    CCF_avg = content("CFloat", min=-1.0)
    CCIwork_avg = content("CFloat", min=-1.0)
    CCIfree_avg = content("CFloat", min=-1.0)
    CCI_avg = content("CFloat", min=-1.0)

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CServalcatPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CServalcatPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CServalcatPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CExpPhasPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['FOM', 'CFOM', 'Hand1Score', 'Hand2Score', 'CC', 'RFactor', 'RFree', 'annotation']
    FOM = content("CFloat", min=0.0)
    CFOM = content("CFloat", min=0.0)
    Hand1Score = content("CFloat")
    Hand2Score = content("CFloat")
    CC = content("CFloat", min=0.0)
    RFactor = content("CFloat", min=0.0)
    RFree = content("CFloat", min=0.0)

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExpPhasPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CExpPhasPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CExpPhasPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CAtomCountPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['nAtoms', 'nResidues']
    nAtoms = content("CInt")
    nResidues = content("CInt")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAtomCountPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """How much model a task produced: ``nAtoms`` and ``nResidues``.

    Inherits from:
    - CAtomCountPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """

    def setFromPdbDataFile(self, pdbDataFile):
        """Count the atoms and residues in *pdbDataFile* into this indicator.

        Both ``chainsaw`` and ``sculptor`` have called this from
        ``processOutputFiles()`` since the Qt days; it was never ported, so
        both raised ``AttributeError`` on every run and both jobs were
        recorded as Finished with the exception printed to the server console.
        See C1 in ``docs/error-handling-remediation.md``.

        Counts every atom and every residue of the first model, waters and
        ligands included -- these are model-preparation tasks, and what the
        indicator reports is the size of the file they produced. The format is
        determined from the file's content, not its name.

        Args:
            pdbDataFile: a ``CPdbDataFile`` (or anything whose ``str()`` is a
                path to a coordinate file).

        Returns:
            self, so a caller may chain.
        """
        import gemmi

        path = None
        if hasattr(pdbDataFile, 'getFullPath'):
            full_path = pdbDataFile.getFullPath()
            if full_path:
                path = str(full_path)
        if not path:
            path = str(pdbDataFile)

        structure = gemmi.read_structure(path, format=gemmi.CoorFormat.Detect)

        nAtoms = 0
        nResidues = 0
        if len(structure):
            for chain in structure[0]:
                nResidues += len(chain)
                for residue in chain:
                    nAtoms += len(residue)

        self.nAtoms.set(nAtoms)
        self.nResidues.set(nResidues)
        return self


class CModelBuildPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['RFactor', 'completeness', 'annotation']
    RFactor = content("CFloat", min=0.0)
    completeness = content("CFloat", min=0.0)

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CModelBuildPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CModelBuildPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CModelBuildPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CRefinementPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['RFactor', 'RFree', 'RMSBond', 'RMSAngle', 'weightUsed', 'annotation']
    RFactor = content("CFloat", min=0.0)
    RFree = content("CFloat", min=0.0)
    RMSBond = content("CFloat", min=0.0)
    RMSAngle = content("CFloat", min=0.0)
    weightUsed = content("CFloat")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefinementPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CRefinementPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CRefinementPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSuperposePerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['RMSxyz', 'nResidues']
    RMSxyz = content("CFloat")
    nResidues = content("CInt")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSuperposePerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CSuperposePerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CSuperposePerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDataReductionRefinementPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['spaceGroup', 'highResLimit', 'rMeas', 'RFactor', 'RFree']
    spaceGroup = content("CSpaceGroup")
    highResLimit = content("CFloat", min=0.0)
    rMeas = content("CFloat")
    RFactor = content("CFloat", min=0.0)
    RFree = content("CFloat", min=0.0)

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDataReductionRefinementPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """
    Composite performance indicator for pipelines spanning data reduction
    and refinement (e.g. dr_mr_modelbuild, substitute_ligand).

    Inherits from:
    - CDataReductionRefinementPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """

    pass


class CPhaseErrorPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['phaseError', 'weightedPhaseError', 'reflectionCorrelation']
    phaseError = content("CFloat", min=0.0)
    weightedPhaseError = content("CFloat", min=0.0)
    reflectionCorrelation = content("CFloat", min=0.0)

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPhaseErrorPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CPhaseErrorPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CPhaseErrorPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CTestObsConversionsPerformance(CPerformanceIndicator):


    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['columnLabelsString']
    columnLabelsString = content("CString")

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CTestObsConversionsPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)
    """

    Inherits from:
    - CTestObsConversionsPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CTestObsConversionsPerformance with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass
