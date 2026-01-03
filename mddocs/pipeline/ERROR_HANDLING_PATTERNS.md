# Error Handling Patterns for CCP4i2 Pipelines

This guide describes the error handling and reporting patterns used in CCP4i2 pipelines and wrappers. Following these patterns ensures consistent error reporting, better debugging, and a better user experience.

## Overview

CCP4i2 provides two complementary systems for error handling:

1. **CErrorReport** - Runtime error tracking during execution
2. **validity()** - Content-aware validation before execution

These systems work together to provide comprehensive error handling.

---

## CErrorReport: Runtime Error Tracking

### What is CErrorReport?

`CErrorReport` is a container for tracking errors and warnings that occur during pipeline execution. It provides:

- Multiple severity levels (OK, UNDEFINED, WARNING, UNDEFINED_ERROR, ERROR)
- Error codes for programmatic handling
- Human-readable descriptions
- XML serialization for reports and diagnostics

### Location

```python
from ccp4i2.core.base_object.error_reporting import CErrorReport, SEVERITY_ERROR, SEVERITY_WARNING
# Or via the legacy import:
from ccp4i2.core import CCP4ErrorHandling
```

### Severity Levels

```python
from ccp4i2.core.base_object.error_reporting import (
    SEVERITY_OK,            # 0 - No error
    SEVERITY_UNDEFINED,     # 1 - Value not set
    SEVERITY_WARNING,       # 2 - Warning (non-blocking)
    SEVERITY_UNDEFINED_ERROR,  # 3 - Required value missing
    SEVERITY_ERROR          # 4 - Fatal error
)
```

### Basic Usage

```python
# Create an error report
error = CErrorReport()

# Add an error
error.append(
    klass="servalcat_pipe",      # Class name for identification
    code=102,                     # Error code from ERROR_CODES
    details="ProSMART failed: connection timeout",
    name="PROSMART_PROTEIN",      # Optional: specific parameter name
    severity=SEVERITY_ERROR       # Default: SEVERITY_ERROR
)

# Check maximum severity
if error.maxSeverity() >= SEVERITY_ERROR:
    print("Fatal errors occurred")

# Get formatted report
print(error.report())
# Output: ERROR in servalcat_pipe 'PROSMART_PROTEIN': ProSMART failed: connection timeout (code 102)

# Merge reports
other_errors = CErrorReport()
other_errors.append(...)
error.extend(other_errors)
```

---

## ERROR_CODES Dictionary Pattern

### Defining Error Codes

Every pipeline should define an `ERROR_CODES` dictionary mapping error codes to descriptions:

```python
class servalcat_pipe(CPluginScript):
    TASKNAME = 'servalcat_pipe'

    ERROR_CODES = {
        101: {'description': 'Failed to generate restraints'},
        102: {'description': 'ProSMART protein restraints failed'},
        103: {'description': 'ProSMART nucleic acid restraints failed'},
        104: {'description': 'MetalCoord restraints failed'},
        105: {'description': 'Servalcat refinement failed'},
        106: {'description': 'Invalid input structure'},
    }
```

### Best Practices for Error Codes

1. **Use ranges** - Assign error code ranges per class (e.g., 100-199 for a pipeline)
2. **Be specific** - Each error code should identify a specific failure mode
3. **Include context** - The description should identify what failed, not why
4. **Document in def.xml** - Error codes should match those in the .def.xml file

---

## try/except Pattern for Pipeline Steps

### Standard Pattern

Wrap each major pipeline step in try/except to capture and report errors:

```python
def startProcess(self, processId):
    """Execute pipeline phases."""

    # Phase 1: Generate protein restraints
    try:
        print("[servalcat_pipe] Phase 1: ProSMART protein restraints")
        self.executeProsmart()
    except Exception as e:
        tb = traceback.format_exc()
        print(f"[servalcat_pipe] Phase 1 FAILED: {e}")
        self.appendErrorReport(102, f'ProSMART protein restraints: {e}\n{tb}')
        self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.FAILED

    # Phase 2: Generate nucleic acid restraints
    try:
        print("[servalcat_pipe] Phase 2: ProSMART nucleic acid restraints")
        self.executeProsmart_na()
    except Exception as e:
        tb = traceback.format_exc()
        print(f"[servalcat_pipe] Phase 2 FAILED: {e}")
        self.appendErrorReport(103, f'ProSMART nucleic acid restraints: {e}\n{tb}')
        self.reportStatus(CPluginScript.FAILED)
        return CPluginScript.FAILED

    # Continue to next phases...
    return CPluginScript.SUCCEEDED
```

### Key Elements

1. **Phase logging** - Print phase start for debugging
2. **Capture traceback** - Use `traceback.format_exc()` for full stack trace
3. **appendErrorReport()** - Record the error with code and details
4. **reportStatus(FAILED)** - Signal failure to pipeline infrastructure
5. **Return early** - Don't continue if a critical step fails

### appendErrorReport() Method

```python
def appendErrorReport(self, code=0, details='', name=None, label=None, cls=None,
                     recordTime=False, stack=True, exc_info=None):
    """
    Append an error to the plugin's error report.

    Args:
        code: Error code number (from ERROR_CODES dictionary)
        details: Error message details (can include exception message and traceback)
        name: Error name (defaults to wrapper name)
        label: Error label
        cls: Class where error occurred
        recordTime: Whether to record timestamp
        stack: Whether to include stack trace
        exc_info: Exception info tuple
    """
```

---

## Real-World Examples

### Example 1: phaser_simple - Setup with Error Handling

```python
def createEnsembleElements(self):
    try:
        from ccp4i2.core.CCP4ModelData import CPdbDataFile, CAtomSelection, CPdbEnsembleItem
        elements = self.container.inputData.ENSEMBLES
        # ... setup ensemble elements ...

    except Exception as e:
        self.appendErrorReport(302, 'Exception setting up search model ensemble: ' + str(e))
        self.reportStatus(CPluginScript.FAILED)
        raise  # Re-raise to stop pipeline

    try:
        if self.container.inputData.INPUT_FIXED.isSet() and ...:
            # ... setup fixed structure ...

    except Exception as e:
        self.appendErrorReport(303, 'Exception setting up fixed structure ensemble: ' + str(e))
        self.reportStatus(CPluginScript.FAILED)
        raise
```

### Example 2: MakeLink - Multiple Validation Checks

```python
def createLinkInstruction(self):
    instruct = "LINK:"

    if not self.container.inputData.ATOM_NAME_1.isSet():
        print("Error - required parameter is not set: ATOM_NAME_1")
        return CPluginScript.FAILED

    if self.container.inputData.MON_1_TYPE.__str__() == 'CIF':
        if not self.container.inputData.RES_NAME_1_CIF.isSet():
            print("Error - required parameter is not set: RES_NAME_1_CIF")
            return CPluginScript.FAILED
        if not self.container.inputData.DICT_1.isSet():
            print("Error - required parameter is not set: DICT_1")
            return CPluginScript.FAILED
        instruct += " RES-NAME-1 " + self.container.inputData.RES_NAME_1_CIF.__str__()
        # ...

    return instruct
```

---

## Recommendations

### Do's

- **Define ERROR_CODES** for all pipelines
- **Wrap major phases** in try/except blocks
- **Include traceback** in error details for debugging
- **Log phase transitions** for debugging
- **Return immediately** after fatal errors
- **Use specific error codes** for different failure modes

### Don'ts

- **Don't swallow exceptions** - Always record them
- **Don't continue** after fatal errors
- **Don't use generic error codes** - Be specific
- **Don't lose context** - Include relevant variable values in error messages

---

## Integration with Diagnostics

Errors are serialized to diagnostic.xml using the `getEtree()` method:

```python
# CErrorReport generates XML like:
# <errorReportList>
#   <errorReport>
#     <class>servalcat_pipe</class>
#     <code>102</code>
#     <details>ProSMART protein restraints: Connection timeout</details>
#     <name>PROSMART_PROTEIN</name>
#     <severity>4</severity>
#     <severityName>ERROR</severityName>
#   </errorReport>
# </errorReportList>
```

This XML appears in the job's diagnostic.xml file and is displayed in the CCP4i2 GUI Diagnostics tab.

---

## See Also

- [VALIDITY_PATTERNS.md](VALIDITY_PATTERNS.md) - Content-aware validation before execution
- [QUICK_REFERENCE.md](../QUICK_REFERENCE.md) - Plugin development patterns
- `core/base_object/error_reporting.py` - CErrorReport implementation
