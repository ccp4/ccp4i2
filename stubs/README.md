# Legacy CCP4-Python Module Stubs

This directory contains minimal stub implementations for legacy ccp4-python modules that are no longer available in modern Python environments.

## Purpose

These stubs allow legacy plugins (particularly `acedrgNew`) to **import successfully** without requiring the full legacy ccp4-python distribution. The stubs provide:

1. **Import compatibility** - Modules can be imported without `ModuleNotFoundError`
2. **Basic constants and classes** - Minimal definitions for code that accesses attributes
3. **Clear error messages** - Methods raise `NotImplementedError` with descriptive messages if called

## Included Stubs

### `ccp4mg.py`
CCP4 Molecular Graphics library stub.
- **Status**: Imported but never actually used by acedrgNew
- **Implementation**: Empty stub (only needs to exist for import)

### `mmdb2.py`
Macromolecular Database library stub.
- **Provides**: Constants (`nAminoacidNames`, `STYPE_RESIDUE`, etc.)
- **Provides**: Stub classes (`Manager`, `intp`)
- **Provides**: Stub functions (`InitMatType`, `getAAProperty`, etc.)
- **Methods**: Raise `NotImplementedError` if called

### `ccp4srs.py`
CCP4 Structure Refinement Suite library stub.
- **Provides**: Constants (`EXTTYPE_Ignore`)
- **Provides**: Stub classes (`Manager`, `Graph`, `GraphMatch`)
- **Methods**: Raise `NotImplementedError` if called

## Installation

The stubs are automatically added to Python's path via `.pth` file in the virtual environment:

```bash
# This is done automatically during setup
echo "/Users/nmemn/Developer/cdata-codegen/stubs" > .venv/lib/python3.11/site-packages/ccp4_stubs.pth
```

## Usage

The stubs are transparent - code imports these modules normally:

```python
import ccp4mg      # Works - stub provides empty module
import mmdb2       # Works - stub provides constants
import ccp4srs     # Works - stub provides classes

# This works (constants and attribute access):
n = mmdb2.nAminoacidNames  # Returns 20

# This fails (actual functionality not implemented):
mgr = mmdb2.Manager()
mgr.ReadCoorFile("file.pdb")  # Raises NotImplementedError
```

## Limitations

### What Works
- **Importing modules** - All three modules can be imported
- **Accessing constants** - Constants like `mmdb2.STYPE_RESIDUE` are available
- **Creating stub objects** - Classes can be instantiated
- **Attribute access** - Basic attributes are available

### What Doesn't Work
- **Atom matching** - `ccp4srs.GraphMatch` matching functionality not implemented
- **PDB reading** - `mmdb2.Manager.ReadCoorFile()` not implemented
- **Structure manipulation** - `ccp4srs.Graph` operations not implemented
- **SRS database** - `ccp4srs.Manager.loadIndex()` does nothing

## Impact on Tests

### Tests That Work
Most acedrg tests work fine because they don't use atom matching:
- `test_from_cif_monomer_library` - Basic dictionary generation
- `test_from_cif_rcsb` - Download and process CIF
- `test_from_smiles` - SMILES to dictionary
- `test_from_mol` - MOL file processing

### Tests That May Fail
Tests requiring atom matching functionality will fail:
- `test_from_smiles_atom_name_matching` - Uses `ATOMMATCHOPTION`
- Any test with `MATCHTLC` parameter

The failure occurs when `atomMatching.matchAtoms()` is called and tries to use real ccp4srs graph matching functionality.

## Maintenance

### Adding New Stubs
If other legacy modules are needed:

1. Create stub file in this directory (e.g., `new_module.py`)
2. Add minimal class/function definitions based on import errors
3. Raise `NotImplementedError` for unimplemented functionality
4. Document limitations in this README

### Updating Existing Stubs
If plugins need additional constants or classes:

1. Check error messages to see what's missing
2. Add minimal definition to appropriate stub
3. Keep implementation minimal - just enough to allow import
4. Update documentation

## Design Philosophy

These stubs follow the principle of **minimum viable compatibility**:

- ✅ **DO** provide enough to allow imports
- ✅ **DO** provide constants and simple attributes
- ✅ **DO** raise clear errors if functionality is invoked
- ❌ **DON'T** try to replicate actual functionality
- ❌ **DON'T** add dependencies to implement features
- ❌ **DON'T** hide errors - fail fast with clear messages

This ensures the codebase remains lightweight while supporting legacy plugin imports.
