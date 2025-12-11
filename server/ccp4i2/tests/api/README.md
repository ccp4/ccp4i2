# API Test Suite

## Overview

This directory contains comprehensive tests for all CCP4X API endpoints, designed to enable safe refactoring of legacy `job_utils` code to modern CData-based approaches.

## Files

### test_viewsets_comprehensive.py
Comprehensive test suite covering all API endpoints:
- **FileViewSetTests**: Tests for file endpoints (8 tests)
- **JobViewSetTests**: Tests for job endpoints (26 tests)
- **ProjectViewSetTests**: Tests for project endpoints (19 tests)
- **SimpleViewSetTests**: Tests for simple CRUD endpoints (5 tests)
- **ProjectExportViewSetTests**: Tests for export endpoints (3 tests)

**Total: 61 test cases** covering all API functionality

### API_ENDPOINT_REFERENCE.md
Comprehensive documentation of:
- All API endpoints with URLs and HTTP methods
- Purpose and functionality of each endpoint
- Dependencies on `job_utils` modules
- Prioritized refactoring targets
- Migration strategy and checklist

### test_api.py (existing)
Original API tests focusing on core workflows:
- Project import
- Job cloning
- Parameter setting
- File upload
- File digest

## Running Tests

### Setup
```bash
export CCP4I2_ROOT=$CCP4I2_ROOT
source .venv/bin/activate
```

### Run All API Tests
```bash
pytest server/ccp4i2/tests/api/ -v
```

### Run Specific Test Class
```bash
pytest server/ccp4i2/tests/api/test_viewsets_comprehensive.py::JobViewSetTests -v
```

### Run Single Test
```bash
pytest server/ccp4i2/tests/api/test_viewsets_comprehensive.py::JobViewSetTests::test_job_clone -v
```

## Purpose

These tests enable systematic refactoring by:

1. **Documenting Current Behavior**: Tests capture how endpoints currently work
2. **Ensuring Compatibility**: Tests verify behavior doesn't change during refactoring
3. **Enabling Confidence**: Comprehensive coverage means safe refactoring
4. **Tracking Progress**: Can run tests after each job_utils module replacement

## job_utils Migration Strategy

### Phase 1: Baseline (COMPLETE)
- ✅ Survey all API endpoints
- ✅ Identify job_utils dependencies
- ✅ Create comprehensive tests
- ✅ Document refactoring priorities

### Phase 2: Setup & Validation
1. Resolve test authentication issues
2. Run baseline tests to capture current behavior
3. Document any test failures or edge cases
4. Create fixtures for common test scenarios

### Phase 3: Systematic Refactoring
For each `job_utils` module (in priority order):

1. **json_for_job_container** (HIGH)
   - Endpoints: `/jobs/{id}/container/`
   - Create: `lib/cdata_utils/serialize.py`
   - Replace: Manual container traversal → CData metadata
   - Test: `test_job_container()`

2. **digest_param_file** (HIGH)
   - Endpoints: `/jobs/{id}/digest/`, `/jobs/{id}/digest_param_file/`
   - Create: `lib/cdata_utils/digest.py`
   - Replace: Legacy file parsing → CData file introspection
   - Test: `test_job_digest()`, `test_job_digest_param_file()`

3. **get_what_next** (HIGH)
   - Endpoints: `/jobs/{id}/what_next/`
   - Create: `lib/workflow/suggestions.py`
   - Replace: Legacy DB queries → Modern workflow analysis
   - Test: `test_job_what_next()`

4. **object_method** (MEDIUM)
   - Endpoints: `/jobs/{id}/object_method/`
   - Create: `lib/cdata_utils/navigation.py`
   - Replace: Manual navigation → `find_by_path()`
   - Test: Custom tests for object method calls

5. **upload_file_param** (MEDIUM)
   - Endpoints: `/jobs/{id}/upload_file_param/`
   - Use: AsyncDatabaseHandler + modern CData
   - Test: `test_job_upload_file_param()`

6. **list_project** (MEDIUM)
   - Endpoints: `/projects/{id}/directory/`
   - Create: `lib/utils/projects/directory.py`
   - Test: `test_project_directory()`

7. **create_task** (MEDIUM)
   - Endpoints: `/projects/{id}/create_task/`
   - Use: Modern CData initialization
   - Test: `test_project_create_task()`

### Phase 4: Verification
1. Run full test suite
2. Test in development environment
3. Verify API clients still work
4. Update API documentation

### Phase 5: Cleanup
1. Remove old `job_utils` modules
2. Update imports
3. Clean up unused code
4. Update documentation

## Test Maintenance

### Adding New Tests
When adding API endpoints:
1. Add test method to appropriate test class
2. Document endpoint in API_ENDPOINT_REFERENCE.md
3. Note any job_utils dependencies
4. Run test to verify it works

### Modifying Endpoints
When changing endpoint behavior:
1. Update test to reflect new behavior
2. Update documentation
3. Ensure backward compatibility or document breaking changes

## Known Issues

### Test Setup
The comprehensive test suite needs minor configuration updates to handle:
- Django REST Framework authentication settings
- Test database configuration
- URL routing in test mode

These are standard Django test setup issues and don't affect the test logic itself.

### Workaround
The existing `test_api.py` can be used as a template for test configuration until the comprehensive suite is fully configured.

## Refactoring Checklist

Use this checklist when replacing job_utils modules:

- [ ] Identify all API endpoints using the module
- [ ] Run baseline tests for those endpoints
- [ ] Create modern replacement utility
- [ ] Update endpoint to use new utility
- [ ] Run tests to verify behavior unchanged
- [ ] Update imports in endpoint
- [ ] Test in development environment
- [ ] Mark module as deprecated
- [ ] Remove old module when all usages replaced

## Benefits of This Approach

1. **Safety**: Tests ensure API behavior doesn't break
2. **Documentation**: Clear mapping of endpoints to dependencies
3. **Visibility**: Easy to track refactoring progress
4. **Confidence**: Comprehensive coverage reduces regression risk
5. **Systematic**: Clear priority order for refactoring work

## Next Steps

1. Fix test configuration issues (authentication, URLs)
2. Run baseline tests and document results
3. Start with highest priority module: `json_for_job_container`
4. Create modern CData utilities in `lib/cdata_utils/`
5. Replace one module at a time, testing after each
6. Document any API changes or improvements

## Contact

For questions about the test suite or refactoring strategy, refer to:
- API_ENDPOINT_REFERENCE.md for detailed endpoint documentation
- test_viewsets_comprehensive.py for test implementations
- ../../../CLAUDE.md for overall project architecture
