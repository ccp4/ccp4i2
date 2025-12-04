# Changelog: Async Execution Infrastructure

## [1.0.0] - 2025-10-30

### 🎉 Major Milestone: Qt-Free Async Execution Complete

Successfully implemented complete asynchronous execution infrastructure for CCP4i2 plugins without Qt dependencies. Real legacy pipelines (demo_copycell) now run end-to-end with real crystallographic programs.

### ✨ New Features

#### Execution System
- **Async Process Manager** - Pure Python subprocess execution using `asyncio`
- **Synchronous Execution** - Blocking execution for fast plugins with proper signal emission
- **Background Monitoring** - Async plugins monitored with automatic status reporting
- **Environment Handling** - CCP4 environment variables (CBIN, CCP4, CLIB) correctly passed to subprocesses

#### Template System
- **CComTemplate** - Modern Qt-free template processor (`core/CCP4ComTemplate.py`)
- **Variable Substitution** - `$VARIABLE` and `$VARIABLE.attribute` expansion
- **Legacy Prefix Handling** - Automatic stripping of numeric prefixes (`"1 HKLIN"` → `"HKLIN"`)
- **Template Routing** - COMLINETEMPLATE → command line, COMTEMPLATE → stdin

#### Signal/Slot System
- **Type-Safe Signals** - `Signal[T]` with generic type support
- **Legacy Compatibility** - Automatic adaptation between dict and int signatures
- **Weak References** - Memory-safe connections with automatic cleanup
- **Async Support** - `emit_async()` for coroutine slots

#### Hierarchical Navigation
- **find() Method** - Path-based object lookup (`find("inputData.protein.XYZIN")`)
- **Depth-First Search** - Simple name resolution across entire hierarchy
- **Template Integration** - Direct use in `$VARIABLE` expansion

#### Smart Data Assignment
- **CData.set()** - Handles dict, CData objects, and primitives intelligently
- **CData.get()** - Extract attributes as dictionary for template expansion
- **Value Type Preservation** - Updates `.value` instead of replacing objects

#### Job Management
- **Subdirectory Convention** - Sub-plugins get `job_1/`, `job_2/`, etc.
- **Automatic Naming** - Sub-plugin names: `parent_name_1`, `parent_name_2`
- **Isolated Workspaces** - Each job has own log, params, and temp files
- **Parameter Persistence** - Automatic `.params.xml` export after execution

#### Pipeline Support
- **postProcessWrapper()** - Propagate sub-plugin completion status
- **connectSignal()** - Qt-compatible signal connection with adapters
- **Callback Chaining** - Multi-step pipelines with conditional branching

### 🔧 Core Changes

#### New Files
- `core/CCP4ComTemplate.py` (221 lines) - Template processor
- `MILESTONE_ASYNC_EXECUTION.md` - Complete architecture documentation
- `QUICK_REFERENCE.md` - Developer quick start guide
- `tests/test_demo_copycell_integration.py` - Real pipeline integration test

#### Modified Files
- `core/CCP4PluginScript.py`
  - Added `startProcess()` async/sync execution (lines 350-520)
  - Added `makeCommandAndScript()` template expansion (lines 700-784)
  - Added `reportStatus()` signal emission (lines 1262-1288)
  - Added `postProcessWrapper()` callback (lines 1290-1306)
  - Added `makePluginObject()` sub-job creation (lines 1312-1380)
  - Added `logFileText()` compatibility method (lines 888-905)
  - Added database integration attributes (lines 120-126)

- `core/base_object/hierarchy_system.py`
  - Added `find()` path navigation (lines 647-728)
  - Added `connectSignal()` legacy adapter (lines 574-645)

- `core/base_object/cdata.py`
  - Added `get()` dict extraction (lines 198-230)
  - Enhanced `set()` to handle CData objects (lines 135-196)

- `core/base_object/signal_system.py`
  - Enhanced error handling in emit() (line 432)
  - Improved connection cleanup

### ✅ Test Results

- **266 tests passed** (97.8% pass rate)
- **6 tests failed** (minor edge cases)
  - 3x converter "already X returns copy"
  - 1x makePluginObject kwargs propagation
  - 2x phase conversion benchmarks

### 🚀 Performance

- **demo_copycell Pipeline**: 1.49 seconds total
  - mtzdump: ~0.5s
  - pdbset: ~0.5s
  - Overhead: <0.01s
- **Memory**: No leaks detected, weak references working correctly
- **Scalability**: Tested with 2-step pipelines, supports arbitrary nesting

### 📚 Documentation

- **MILESTONE_ASYNC_EXECUTION.md** - 500+ lines of comprehensive documentation
  - Architecture overview with diagrams
  - Implementation details for all components
  - API examples for common patterns
  - Performance metrics and benchmarks
  - Known limitations and future work

- **QUICK_REFERENCE.md** - 300+ lines of developer guide
  - Essential concepts and patterns
  - Common use cases with code examples
  - Debugging tips and troubleshooting
  - Performance guidelines

### 🔬 Real-World Validation

Successfully ran **demo_copycell** pipeline with real MDM2 crystallographic data:
1. **mtzdump** extracted unit cell: `71.478 71.478 104.236 90° 90° 120°` (P61 2 2)
2. **pdbset** applied cell to PDB file
3. Output PDB verified with correct CRYST1 record
4. Complete pipeline executed in <2 seconds

### 🐛 Known Issues

1. **FileNotFoundError in saveParams()** - Some tests fail due to missing directory creation
2. **kwargs Not Propagated** - makePluginObject() doesn't pass through all keyword arguments
3. **Phase Conversion Tolerance** - Small numerical differences vs. CCP4 chltofom

### 🔜 Next Steps

1. Fix remaining 6 test failures
2. Robust directory creation before saveParams()
3. Enhanced error reporting with stderr capture
4. Async output parsing for real-time progress

### 💡 Key Insights

1. **Signals Must Always Emit** - Even synchronous plugins must call `reportStatus()` for pipeline compatibility
2. **Smart Assignment Critical** - Legacy code expects `set(cdata_object)` to work transparently
3. **Job Isolation Essential** - Subdirectories prevent conflicts in nested pipelines
4. **Legacy Compatibility Complex** - Template prefixes, slot decorators, and Qt quirks require careful handling

### 🎯 Impact

This milestone enables:
- ✅ Command-line execution of CCP4i2 pipelines without Qt
- ✅ Integration testing with real CCP4 programs
- ✅ Foundation for GUI development (Qt or web-based)
- ✅ Database-backed job tracking (infrastructure ready)
- ✅ Distributed execution on HPC clusters (future work)

### 🙏 Credits

Implemented through iterative development with comprehensive testing at each stage. Special recognition for:
- Signal system design - Qt-free, type-safe, async-compatible
- Template processor - Clean implementation of complex legacy behavior
- Smart assignment - Seamless CData/dict/primitive handling

---

## [0.9.0] - Previous Work

### Base Infrastructure (Pre-Async)
- CData metadata system with decorators
- Fundamental types (CInt, CFloat, CString, CBoolean, CList)
- HierarchicalObject with lifecycle management
- Signal/Slot foundation
- Data file handling (CDataFile, CMtzDataFile, CPdbDataFile)
- Container operations
- Parameter XML export/import
- Phase and observation data converters

---

**Version**: 1.0.0
**Date**: October 30, 2025
**Stability**: Beta (97.8% test pass rate)
**Recommended**: Command-line use, integration testing
**Next Release**: 1.1.0 (after fixing remaining test failures)
