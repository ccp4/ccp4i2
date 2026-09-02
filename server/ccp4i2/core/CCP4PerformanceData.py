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


    """
    A key performance indicator harvested from a job into the project
    database. The base case is a single float with a label; subclasses
    declare named metrics instead.
    """

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
    value = content(
        "CFloat",
        min=0.0,
        guiLabel='Value',
        toolTip='The indicator itself, where a single number suffices; typed subclasses declare named metrics instead')
    annotation = content(
        "CString",
        guiLabel='Annotation',
        toolTip='Human-readable description')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPerformanceIndicator.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

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


    """
    Inherits from:
    - CDataReductionPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['spaceGroup', 'highResLimit', 'rMeas']
    spaceGroup = content(
        "CSpaceGroup",
        guiLabel='Space group',
        toolTip='Space group of the crystal')
    highResLimit = content(
        "CFloat",
        min=0.0,
        guiLabel='High resolution',
        toolTip='High resolution limit of the reduced data, in Angstroms')
    rMeas = content(
        "CFloat",
        guiLabel='Rmeas',
        toolTip='Redundancy-independent merging R factor')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDataReductionPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CPairefPerformance(CPerformanceIndicator):


    """
    Inherits from:
    - CPairefPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['cutoff']
    cutoff = content(
        "CFloat",
        min=0.0,
        guiLabel='Resolution cutoff',
        toolTip='Resolution at which paired refinement stopped improving, in Angstroms')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPairefPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CDataReductionCCPerformance(CPerformanceIndicator):


    """
    Inherits from:
    - CDataReductionCCPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['spaceGroup', 'highResLimit', 'ccHalf']
    spaceGroup = content(
        "CSpaceGroup",
        guiLabel='Space group',
        toolTip='Space group of the crystal')
    highResLimit = content(
        "CFloat",
        min=0.0,
        guiLabel='High resolution',
        toolTip='High resolution limit of the reduced data, in Angstroms')
    ccHalf = content(
        "CFloat",
        guiLabel='CC half',
        toolTip='Correlation between two half-datasets, the usual resolution criterion')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDataReductionCCPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CServalcatPerformance(CPerformanceIndicator):


    """
    Inherits from:
    - CServalcatPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['RFactor', 'RFree', 'R', 'R1Factor', 'R1Free', 'R1', 'FSCaverage', 'annotation']
    RFactor = content(
        "CFloat",
        min=0.0,
        guiLabel='R factor',
        toolTip='R factor against the working set')
    RFree = content(
        "CFloat",
        min=0.0,
        guiLabel='R free',
        toolTip='R factor against the free set, the unbiased measure of fit')
    R = content(
        "CFloat",
        min=0.0,
        guiLabel='R',
        toolTip='R factor over all reflections')
    R1Factor = content(
        "CFloat",
        min=0.0,
        guiLabel='R1 factor',
        toolTip='R1 factor against the working set')
    R1Free = content(
        "CFloat",
        min=0.0,
        guiLabel='R1 free',
        toolTip='R1 factor against the free set')
    R1 = content(
        "CFloat",
        min=0.0,
        guiLabel='R1',
        toolTip='R1 factor over all reflections')
    FSCaverage = content(
        "CFloat",
        min=-1.0,
        guiLabel='FSC average',
        toolTip='Mean Fourier shell correlation between model and map')
    CCFwork_avg = content(
        "CFloat",
        min=-1.0,
        guiLabel='CCF work',
        toolTip='Mean correlation of amplitudes for the working set')
    CCFfree_avg = content(
        "CFloat",
        min=-1.0,
        guiLabel='CCF free',
        toolTip='Mean correlation of amplitudes for the free set')
    CCF_avg = content(
        "CFloat",
        min=-1.0,
        guiLabel='CCF',
        toolTip='Mean correlation of amplitudes over all reflections')
    CCIwork_avg = content(
        "CFloat",
        min=-1.0,
        guiLabel='CCI work',
        toolTip='Mean correlation of intensities for the working set')
    CCIfree_avg = content(
        "CFloat",
        min=-1.0,
        guiLabel='CCI free',
        toolTip='Mean correlation of intensities for the free set')
    CCI_avg = content(
        "CFloat",
        min=-1.0,
        guiLabel='CCI',
        toolTip='Mean correlation of intensities over all reflections')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CServalcatPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CExpPhasPerformance(CPerformanceIndicator):


    """
    Inherits from:
    - CExpPhasPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['FOM', 'CFOM', 'Hand1Score', 'Hand2Score', 'CC', 'RFactor', 'RFree', 'annotation']
    FOM = content(
        "CFloat",
        min=0.0,
        guiLabel='FOM',
        toolTip='Mean figure of merit of the phases')
    CFOM = content(
        "CFloat",
        min=0.0,
        guiLabel='Combined FOM',
        toolTip='Combined figure of merit over both hands')
    Hand1Score = content(
        "CFloat",
        guiLabel='Hand 1 score',
        toolTip='Score for the first of the two possible hands')
    Hand2Score = content(
        "CFloat",
        guiLabel='Hand 2 score',
        toolTip='Score for the second of the two possible hands')
    CC = content(
        "CFloat",
        min=0.0,
        guiLabel='CC',
        toolTip='Correlation coefficient of the experimental phasing')
    RFactor = content(
        "CFloat",
        min=0.0,
        guiLabel='R factor',
        toolTip='R factor of the phasing model against the data')
    RFree = content(
        "CFloat",
        min=0.0,
        guiLabel='R free',
        toolTip='R factor against the free set')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CExpPhasPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CAtomCountPerformance(CPerformanceIndicator):


    """How much model a task produced: ``nAtoms`` and ``nResidues``.

    Inherits from:
    - CAtomCountPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['nAtoms', 'nResidues']
    nAtoms = content(
        "CInt",
        guiLabel='Atoms',
        toolTip='Number of atoms in the model produced')
    nResidues = content(
        "CInt",
        guiLabel='Residues',
        toolTip='Number of residues in the model produced')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CAtomCountPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

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


    """
    Inherits from:
    - CModelBuildPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['RFactor', 'completeness', 'annotation']
    RFactor = content(
        "CFloat",
        min=0.0,
        guiLabel='R factor',
        toolTip='R factor of the built model against the data')
    completeness = content(
        "CFloat",
        min=0.0,
        guiLabel='Completeness',
        toolTip='Fraction of the expected model that was built')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CModelBuildPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CRefinementPerformance(CPerformanceIndicator):


    """
    Inherits from:
    - CRefinementPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['RFactor', 'RFree', 'RMSBond', 'RMSAngle', 'weightUsed', 'annotation']
    RFactor = content(
        "CFloat",
        min=0.0,
        guiLabel='R factor',
        toolTip='R factor against the working set')
    RFree = content(
        "CFloat",
        min=0.0,
        guiLabel='R free',
        toolTip='R factor against the free set, the unbiased measure of fit')
    RMSBond = content(
        "CFloat",
        min=0.0,
        guiLabel='RMS bond',
        toolTip='Root-mean-square deviation of bond lengths from ideal, in Angstroms')
    RMSAngle = content(
        "CFloat",
        min=0.0,
        guiLabel='RMS angle',
        toolTip='Root-mean-square deviation of bond angles from ideal, in degrees')
    weightUsed = content(
        "CFloat",
        guiLabel='Weight',
        toolTip='Weight applied between the X-ray and geometry terms')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CRefinementPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CSuperposePerformance(CPerformanceIndicator):


    """
    Inherits from:
    - CSuperposePerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['RMSxyz', 'nResidues']
    RMSxyz = content(
        "CFloat",
        guiLabel='RMS deviation',
        toolTip='Root-mean-square distance between superposed atoms, in Angstroms')
    nResidues = content(
        "CInt",
        guiLabel='Residues',
        toolTip='Number of residues matched in the superposition')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CSuperposePerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CDataReductionRefinementPerformance(CPerformanceIndicator):


    """
    Composite performance indicator for pipelines spanning data reduction
    and refinement (e.g. dr_mr_modelbuild, substitute_ligand).

    Inherits from:
    - CDataReductionRefinementPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['spaceGroup', 'highResLimit', 'rMeas', 'RFactor', 'RFree']
    spaceGroup = content(
        "CSpaceGroup",
        guiLabel='Space group',
        toolTip='Space group of the crystal')
    highResLimit = content(
        "CFloat",
        min=0.0,
        guiLabel='High resolution',
        toolTip='High resolution limit of the data, in Angstroms')
    rMeas = content(
        "CFloat",
        guiLabel='Rmeas',
        toolTip='Redundancy-independent merging R factor')
    RFactor = content(
        "CFloat",
        min=0.0,
        guiLabel='R factor',
        toolTip='R factor against the working set after refinement')
    RFree = content(
        "CFloat",
        min=0.0,
        guiLabel='R free',
        toolTip='R factor against the free set after refinement')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CDataReductionRefinementPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    pass


class CPhaseErrorPerformance(CPerformanceIndicator):


    """
    Inherits from:
    - CPhaseErrorPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['phaseError', 'weightedPhaseError', 'reflectionCorrelation']
    phaseError = content(
        "CFloat",
        min=0.0,
        guiLabel='Phase error',
        toolTip='Mean difference between calculated and reference phases, in degrees')
    weightedPhaseError = content(
        "CFloat",
        min=0.0,
        guiLabel='Weighted phase error',
        toolTip='Phase error weighted by structure factor amplitude, in degrees')
    reflectionCorrelation = content(
        "CFloat",
        min=0.0,
        guiLabel='Reflection correlation',
        toolTip='Correlation between calculated and reference structure factors')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CPhaseErrorPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass


class CTestObsConversionsPerformance(CPerformanceIndicator):


    """
    Inherits from:
    - CTestObsConversionsPerformance: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """
    class Meta:
        qualifiers = {
            "allowUndefined": True,
            "guiDefinition": {},
            "saveToDb": False,
        }
        contents_order = ['columnLabelsString']
    columnLabelsString = content(
        "CString",
        guiLabel='Columns',
        toolTip='Column labels produced by the conversion')

    def __init__(self, parent=None, name=None, **kwargs):
        """
        Initialize CTestObsConversionsPerformance.

        Args:
            parent: Parent object in hierarchy
            name: Object name
            **kwargs: Additional keyword arguments
        """
        super().__init__(parent=parent, name=name, **kwargs)

    # Add your methods here
    pass
