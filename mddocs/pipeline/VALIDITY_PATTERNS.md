# Validity Method Patterns for CCP4i2 Pipelines

This guide describes how to implement content-aware validation in CCP4i2 pipelines using the `validity()` method override pattern.

## Overview

The `validity()` method provides pre-execution validation that can:

1. Check parameter constraints defined in .def.xml
2. Verify input file existence (mustExist files)
3. Apply **content-aware logic** based on actual parameter values
4. Adjust validation rules dynamically (allowUndefined, etc.)
5. Filter or suppress expected validation errors

---

## Base validity() Method

The base implementation in `CPluginScript.validity()`:

```python
def validity(self) -> CErrorReport:
    """
    Validate the plugin's container and return an error report.

    The validation includes:
    - Container validity (all child objects' validity() methods)
    - Input file existence checks (mustExist files)

    Returns:
        CErrorReport containing all validation errors/warnings
    """
    error = CErrorReport()

    # First, get container validation
    if hasattr(self, 'container') and self.container is not None:
        container_errors = self.container.validity()
        if container_errors:
            error.extend(container_errors)

    # Then add input data checks (file existence etc.)
    input_errors = self.checkInputData()
    if input_errors:
        error.extend(input_errors)

    return error
```

---

## When to Override validity()

Override `validity()` when:

1. **Conditional field requirements** - Some fields are required only based on other parameter values
2. **Programmatically populated fields** - Fields that are empty at validation time but filled during process()
3. **Embedded wrapper adjustments** - Child wrappers need qualifier adjustments in pipeline context
4. **Custom validation logic** - Domain-specific validation not expressible in .def.xml

---

## Pattern 1: Adjusting allowUndefined Based on Mode Selection

When a pipeline has multiple modes (e.g., different input types), adjust which fields are required based on the selected mode.

### Example: MakeLink - TLC vs CIF Mode

```python
class MakeLink(CPluginScript):
    TASKNAME = 'MakeLink'

    def validity(self):
        """Override to adjust allowUndefined based on MON_TYPE selection.

        MakeLink has conditional field requirements:
        - When MON_1_TYPE='TLC', the TLC fields are required and CIF fields are optional
        - When MON_1_TYPE='CIF', the CIF fields are required and TLC fields are optional
        """
        inp = self.container.inputData
        ctrl = self.container.controlParameters

        # Determine which mode we're in (defaults to TLC from def.xml)
        mon_1_type = str(inp.MON_1_TYPE) if inp.MON_1_TYPE.isSet() else 'TLC'
        mon_2_type = str(inp.MON_2_TYPE) if inp.MON_2_TYPE.isSet() else 'TLC'

        # Set allowUndefined for fields based on mode
        # For monomer 1
        if mon_1_type == 'TLC':
            # TLC mode: CIF fields are optional
            inp.RES_NAME_1_CIF.setQualifier('allowUndefined', True)
            inp.ATOM_NAME_1_CIF.setQualifier('allowUndefined', True)
        else:
            # CIF mode: TLC fields are optional
            inp.RES_NAME_1_TLC.setQualifier('allowUndefined', True)
            inp.ATOM_NAME_1_TLC.setQualifier('allowUndefined', True)

        # For monomer 2 (similar pattern)
        if mon_2_type == 'TLC':
            inp.RES_NAME_2_CIF.setQualifier('allowUndefined', True)
            inp.ATOM_NAME_2_CIF.setQualifier('allowUndefined', True)
        else:
            inp.RES_NAME_2_TLC.setQualifier('allowUndefined', True)
            inp.ATOM_NAME_2_TLC.setQualifier('allowUndefined', True)

        # LIST fields are populated by GUI dropdowns - make optional for CLI use
        inp.DELETE_1_LIST.setQualifier('allowUndefined', True)
        inp.DELETE_2_LIST.setQualifier('allowUndefined', True)
        inp.CHARGE_1_LIST.setQualifier('allowUndefined', True)
        inp.CHARGE_2_LIST.setQualifier('allowUndefined', True)
        inp.CHANGE_BOND_1_LIST.setQualifier('allowUndefined', True)
        inp.CHANGE_BOND_2_LIST.setQualifier('allowUndefined', True)
        ctrl.MODEL_RES_LIST.setQualifier('allowUndefined', True)

        # Now call parent validity() which will use our updated allowUndefined settings
        return super(MakeLink, self).validity()
```

### Key Points

1. **Adjust qualifiers BEFORE calling super()** - The parent's validity() uses these qualifiers
2. **Check isSet()** - Use safe defaults when parameters aren't set
3. **Handle both directions** - Make one set optional when the other is selected

---

## Pattern 2: Filtering Expected Validation Errors

When a field is intentionally empty at validation time because it will be populated programmatically during process(), filter out that specific error.

### Example: phaser_simple - Programmatically Populated List

```python
class phaser_simple(phaser_pipeline.phaser_pipeline):
    TASKNAME = 'phaser_simple'

    ERROR_CODES = {
        301: {'description': 'Exception in createEnsembleElements'},
        302: {'description': 'Exception setting up search model ensemble'},
        303: {'description': 'Exception setting up fixed structure ensemble'},
    }

    def validity(self):
        """Override to filter out ENSEMBLES list length error.

        ENSEMBLES is intentionally empty at validation time because it's
        populated programmatically by createEnsembleElements() during process().
        """
        from ccp4i2.core import CCP4ErrorHandling

        # Get parent validation
        error = super(phaser_simple, self).validity()

        # Filter out the ENSEMBLES minimum length error (code 101)
        # This error is expected since ENSEMBLES is populated in createEnsembleElements()
        filtered = CCP4ErrorHandling.CErrorReport()
        for err in error.getErrors():
            # Skip error code 101 (min list length) for ENSEMBLES
            if err.get('code') == 101 and 'ENSEMBLES' in err.get('name', ''):
                continue
            filtered.append(
                klass=err.get('class', ''),
                code=err.get('code', 0),
                details=err.get('details', ''),
                name=err.get('name', ''),
                severity=err.get('severity', 0)
            )

        return filtered

    def process(self):
        self.createEnsembleElements()  # Populates ENSEMBLES
        super(phaser_simple, self).process()
```

### Key Points

1. **Call super() FIRST** - Get the full validation report
2. **Create new filtered report** - Copy valid errors, skip expected ones
3. **Filter specifically** - Match both error code AND parameter name
4. **Document why** - Explain why the error is expected

---

## Pattern 3: Adjusting Embedded Wrapper Qualifiers

When a pipeline embeds another wrapper, adjust the embedded wrapper's validation to account for the pipeline context.

### Example: servalcat_pipe - Embedded metalCoordWrapper

```python
class servalcat_pipe(CPluginScript):
    TASKNAME = 'servalcat_pipe'

    def validity(self):
        """
        Validate plugin, adjusting qualifiers for embedded wrappers first.

        The metalCoordWrapper embeds metalCoord, which requires XYZIN when run standalone.
        However, in the context of servalcat_pipe, XYZIN is programmatically filled from
        the pipeline's own XYZIN in executeMetalCoord(), so we set allowUndefined=True.

        This override is called by the i2run workflow after parameters are loaded from
        input_params.xml, allowing qualifier adjustments to be based on actual parameter values.
        """
        # Adjust qualifiers for embedded metalCoordWrapper inputs
        if hasattr(self.container, 'metalCoordWrapper'):
            self.container.metalCoordWrapper.inputData.XYZIN.set_qualifier('allowUndefined', True)

        # Call parent validity
        return super(servalcat_pipe, self).validity()
```

### Key Points

1. **Check hasattr()** - The embedded wrapper may not always exist
2. **Use set_qualifier()** - Adjust the specific qualifier
3. **Context matters** - Fields required standalone may be optional in pipeline context

---

## Pattern 4: checkInputData Override

For file-specific validation beyond mustExist, override `checkInputData()`:

### Example: phaser_simple - Conditional File Validation

```python
def checkInputData(self):
    invalidFiles = super(phaser_simple, self).checkInputData()

    # If INPUT_FIXED is not set, don't require XYZIN_FIXED
    if (not self.container.inputData.INPUT_FIXED) and ('XYZIN_FIXED' in invalidFiles):
        invalidFiles.remove('XYZIN_FIXED')

    return invalidFiles
```

---

## Validation Flow

```
Plugin.validity() called (e.g., by validate_job API)
        ↓
Override adjusts qualifiers (allowUndefined, etc.)
        ↓
super().validity() called
        ↓
Container.validity() validates all children
        ↓
checkInputData() verifies file existence
        ↓
Override filters/modifies error report
        ↓
Final CErrorReport returned
```

---

## Common Qualifiers

| Qualifier | Type | Description |
|-----------|------|-------------|
| `allowUndefined` | bool | If False, unset value is an error |
| `mustExist` | bool | File must exist on filesystem |
| `min` / `max` | number | Range constraints for numeric values |
| `listMinLength` / `listMaxLength` | int | List length constraints |
| `onlyEnumerators` | bool | Value must be from enumerated list |

### Setting Qualifiers

```python
# Single qualifier
obj.setQualifier('allowUndefined', True)
# or
obj.set_qualifier('allowUndefined', True)

# Multiple qualifiers
obj.setQualifiers({
    'allowUndefined': True,
    'mustExist': False
})
```

---

## Testing validity() Overrides

```python
def test_validity_mode_switching():
    """Test that validity adjusts requirements based on mode."""
    plugin = MakeLink(workDirectory='/tmp/test', name='test')

    # Set to TLC mode
    plugin.container.inputData.MON_1_TYPE.set('TLC')
    plugin.container.inputData.RES_NAME_1_TLC.set('ALA')
    # Don't set RES_NAME_1_CIF - should be optional

    errors = plugin.validity()

    # Should not have error for missing RES_NAME_1_CIF
    for err in errors.getErrors():
        assert 'RES_NAME_1_CIF' not in err.get('name', ''), \
            "RES_NAME_1_CIF should be optional in TLC mode"
```

---

## Best Practices

### Do's

- **Document the override** - Explain why custom validation is needed
- **Call super()** - Preserve parent validation logic
- **Be specific** - Filter only the exact errors that are expected
- **Test both modes** - Verify validation works for all configurations

### Don'ts

- **Don't skip super()** - Unless you're completely replacing validation
- **Don't blindly filter** - Only filter errors you explicitly expect
- **Don't modify global state** - Qualifier changes should be local to validation
- **Don't ignore warnings** - They may indicate real issues

---

## See Also

- [ERROR_HANDLING_PATTERNS.md](ERROR_HANDLING_PATTERNS.md) - Runtime error handling with CErrorReport
- [QUICK_REFERENCE.md](../QUICK_REFERENCE.md) - Plugin development patterns
- `core/CCP4PluginScript.py` - Base validity() implementation
- `core/base_object/error_reporting.py` - CErrorReport class
