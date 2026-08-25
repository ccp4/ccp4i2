"""
Implementation classes for CCP4PerformanceData.py

Extends stub classes from ccp4i2.core.cdata_stubs with methods and business logic.
This file is safe to edit - add your implementation code here.
"""

from typing import Optional, Any

from ccp4i2.core.cdata_stubs.CCP4PerformanceData import CAtomCountPerformanceStub, CDataReductionCCPerformanceStub, CDataReductionPerformanceStub, CDataReductionRefinementPerformanceStub, CExpPhasPerformanceStub, CModelBuildPerformanceStub, CPairefPerformanceStub, CPerformanceIndicatorStub, CPhaseErrorPerformanceStub, CRefinementPerformanceStub, CServalcatPerformanceStub, CSuperposePerformanceStub, CTestObsConversionsPerformanceStub


# Define CPerformanceIndicator FIRST since other classes inherit from it
class CPerformanceIndicator(CPerformanceIndicatorStub):
    """
    Extends CPerformanceIndicatorStub with implementation-specific methods.

    Overrides isSet() so that a performance indicator is considered "set"
    when any of its KPI children (beyond the inherited 'value' and
    'annotation' fields) have been explicitly set.  This ensures
    exclude_unset serialization preserves PERFORMANCE elements that
    contain real metrics like highResLimit, rMeas, RFactor, etc.
    """

    # Fields inherited from CPerformanceIndicatorStub that are structural,
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


class CAtomCountPerformance(CAtomCountPerformanceStub, CPerformanceIndicator):
    """How much model a task produced: ``nAtoms`` and ``nResidues``.

    Inherits from:
    - CAtomCountPerformanceStub: Metadata and structure
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


class CDataReductionCCPerformance(CDataReductionCCPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CDataReductionCCPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CDataReductionCCPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDataReductionPerformance(CDataReductionPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CDataReductionPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CDataReductionPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CDataReductionRefinementPerformance(CDataReductionRefinementPerformanceStub, CPerformanceIndicator):
    """
    Composite performance indicator for pipelines spanning data reduction
    and refinement (e.g. dr_mr_modelbuild, substitute_ligand).

    Inherits from:
    - CDataReductionRefinementPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    """

    pass


class CExpPhasPerformance(CExpPhasPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CExpPhasPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CExpPhasPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CModelBuildPerformance(CModelBuildPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CModelBuildPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CModelBuildPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CPairefPerformance(CPairefPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CPairefPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CPairefPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CPhaseErrorPerformance(CPhaseErrorPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CPhaseErrorPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CPhaseErrorPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CRefinementPerformance(CRefinementPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CRefinementPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CRefinementPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CServalcatPerformance(CServalcatPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CServalcatPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CServalcatPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CSuperposePerformance(CSuperposePerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CSuperposePerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CSuperposePerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass


class CTestObsConversionsPerformance(CTestObsConversionsPerformanceStub, CPerformanceIndicator):
    """

    Inherits from:
    - CTestObsConversionsPerformanceStub: Metadata and structure
    - CPerformanceIndicator: Shared full-fat methods
    Extends CTestObsConversionsPerformanceStub with implementation-specific methods.
    Add file I/O, validation, and business logic here.
    """

    # Add your methods here
    pass
