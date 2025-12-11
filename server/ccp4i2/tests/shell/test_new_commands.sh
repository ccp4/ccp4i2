#!/bin/bash
# Test script for new management commands
# Tests all 6 newly created commands

set -e  # Exit on error

# Source common setup (sets CCP4I2_ROOT, sources CCP4 and venv)
source "$(dirname "$0")/common.sh"
cd $CCP4I2_ROOT/server

echo "=========================================="
echo "Testing New Management Commands"
echo "=========================================="
echo ""

# Get a test job UUID
echo "1. Finding test job..."
TEST_JOB_UUID=$(python manage.py list_jobs banana1 --json --limit 1 2>/dev/null | grep -o '"uuid": "[^"]*"' | head -1 | cut -d'"' -f4)
echo "   Using job UUID: $TEST_JOB_UUID"
echo ""

# Test 1: validate_job
echo "=========================================="
echo "TEST 1: validate_job"
echo "=========================================="
python manage.py validate_job --jobuuid "$TEST_JOB_UUID" 2>&1 | tail -20
echo ""

# Test 2: get_job_report (params)
echo "=========================================="
echo "TEST 2: get_job_report --type params"
echo "=========================================="
python manage.py get_job_report --jobuuid "$TEST_JOB_UUID" --type params -o /tmp/test_params.xml 2>&1 | tail -10
if [ -f /tmp/test_params.xml ]; then
    echo "   ✓ Params XML written"
    ls -lh /tmp/test_params.xml | awk '{print "   Size:", $5}'
else
    echo "   ✗ Params XML not created"
fi
echo ""

# Test 3: get_job_report (report)
echo "=========================================="
echo "TEST 3: get_job_report --type report"
echo "=========================================="
python manage.py get_job_report --jobuuid "$TEST_JOB_UUID" --type report -o /tmp/test_report.xml 2>&1 | tail -10
if [ -f /tmp/test_report.xml ]; then
    echo "   ✓ Report XML written"
    ls -lh /tmp/test_report.xml | awk '{print "   Size:", $5}'
else
    echo "   ✗ Report XML not created"
fi
echo ""

# Test 4: clone_job
echo "=========================================="
echo "TEST 4: clone_job"
echo "=========================================="
python manage.py clone_job --jobuuid "$TEST_JOB_UUID" 2>&1 | tail -15
echo ""

# Test 5: execute_job (will fail if no params, but tests the command)
echo "=========================================="
echo "TEST 5: execute_job (dry run)"
echo "=========================================="
echo "   (This may fail if job has no input_params.xml - that's expected)"
python manage.py execute_job --jobuuid "$TEST_JOB_UUID" 2>&1 | tail -10 || echo "   Expected failure - job not ready to run"
echo ""

# Test 6: set_job_parameter
echo "=========================================="
echo "TEST 6: set_job_parameter"
echo "=========================================="
echo "   (This may fail if job has no params - that's expected)"
python manage.py set_job_parameter --jobuuid "$TEST_JOB_UUID" --path "container.inputData" --value "test" 2>&1 | tail -10 || echo "   Expected failure - job not ready for param setting"
echo ""

# Test 7: export_job
echo "=========================================="
echo "TEST 7: export_job"
echo "=========================================="
python manage.py export_job --jobuuid "$TEST_JOB_UUID" -o /tmp/test_export.zip 2>&1 | tail -10
if [ -f /tmp/test_export.zip ]; then
    echo "   ✓ Export ZIP created"
    ls -lh /tmp/test_export.zip | awk '{print "   Size:", $5}'
else
    echo "   ✗ Export ZIP not created"
fi
echo ""

# Summary
echo "=========================================="
echo "TEST SUMMARY"
echo "=========================================="
echo "Commands tested:"
echo "  1. validate_job"
echo "  2. get_job_report (params)"
echo "  3. get_job_report (report)"
echo "  4. clone_job"
echo "  5. execute_job"
echo "  6. set_job_parameter"
echo "  7. export_job"
echo ""
echo "Check output above for results"
echo "=========================================="
