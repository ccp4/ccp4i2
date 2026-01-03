# CData Code Generation System

This directory contains the **modern production code generator** for CData classes.

## Active Files

### Core Generator Files
- **`production_generator.py`** - Main production code generator
  - Generates complete, production-ready Python classes from `cdata.json`
  - Handles imports, topological sorting, and complete class bodies
  - Run with: `python production_generator.py [--output DIR] [--report]`

- **`type_resolver.py`** - Type resolution system
  - Resolves type references (CUUID, CFilePath, etc.) to proper imports
  - Determines which module each type lives in
  - Generates correct import statements

- **`class_graph.py`** - Dependency graph and topological sorting
  - Builds complete dependency graph across all classes
  - Performs topological sorting (global and per-file)
  - Ensures parent classes always come before children

### Data Files
- **`cdata.json`** - Source metadata scanned from CCP4i2 codebase
  - Contains definitions for 204+ classes
  - Updated from the moving CCP4i2 codebase

## Usage

```bash
# Generate production code (default: core/generated/)
python production_generator.py

# Generate with custom output directory
python production_generator.py --output core/custom_output

# Show dependency analysis report
python production_generator.py --report

# Verify imports after generation
python production_generator.py --verify
```

## Generated Output

The generator creates:
- 13 Python files in the output directory
- Each file contains properly sorted classes with:
  - Correct imports for all type references
  - Complete `__init__` methods with docstrings
  - Type-annotated attributes
  - Complete `@cdata_class` decorators
  - Auto-formatting with autopep8

Example: `core/generated/CCP4Data.py` contains 15 classes with all dependencies resolved.

## Obsolete Files (Archived)

The following obsolete files have been moved to `_obsolete_backup/`:
- `build_lookup.py` - Old lookup building system
- `cdata_class_generator.py` - Early generator attempt
- `generate_new_files.py` - Old stub generator (superseded)
- `metadata_loader.py` - Old metadata loading
- `test_cdata_class_generator.py` - Test for obsolete generator
- `test_cdata_class_generator_multifile.py` - Test for obsolete generator
- `cdata_requirements.txt` - Empty requirements file

These files are preserved for reference but are no longer used in the build process.

## Migration History

**Original Problem:** The old generator (`generate_new_files.py`) produced incomplete stubs with missing imports, causing `NameError: name 'CUUID' is not defined`.

**Solution:** Complete rewrite as `production_generator.py` with:
- ✅ Automatic type resolution and import generation
- ✅ Global + per-file topological sorting
- ✅ Complete class generation (not stubs)
- ✅ Zero manual patching required
- ✅ Works with moving CCP4i2 codebase

**Results:**
- Before: 1 import error, manual patching required
- After: 0 import errors, fully automated generation

## Dependencies

The production generator requires:
- Python 3.11+
- `autopep8` for code formatting (optional but recommended)

Install with: `pip install autopep8`
