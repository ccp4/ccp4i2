"""
API-driven tests for CCP4i2 pipelines and wrappers.

This package provides comprehensive API tests that exercise the full
Django REST API workflow:

1. Create project
2. Create job for a specific task
3. Set parameters and upload input files
4. Run job and wait for completion
5. Validate outputs via digest endpoints

Test Modules:
- test_data_reduction_api: aimless_pipe, import_merged, freerflag
- test_mr_pipelines_api: phaser_simple, molrep_pipe, i2Dimple, mrbump_basic
- test_refinement_api: prosmart_refmac, servalcat_pipe, sheetbend
- test_ep_pipelines_api: shelx, crank2, phaser_EP
- test_model_building_api: modelcraft, parrot, shelxeMR, acorn
- test_utilities_api: LidiaAcedrgNew, validate_protein, AUSPEX, etc.

Running Tests:
    # Run all API tests
    pytest tests/api/ -v

    # Run specific test file
    pytest tests/api/test_data_reduction_api.py -v

    # Run specific test
    pytest tests/api/test_data_reduction_api.py::TestAimlessPipeAPI::test_gamma -v

Infrastructure:
- base.py: APITestBase class with helper methods
- conftest.py: Fixtures for demo data and file-based database
"""
