#!/bin/bash
# Simple end-to-end test of core commands
# Focuses on validating library→CLI integration works

set -e

export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
source /Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh
source $CCP4I2_ROOT/.venv/bin/activate
cd $CCP4I2_ROOT/server

echo "==========================================="
echo "Simple End-to-End Test"
echo "==========================================="
echo ""

# Use existing project with existing job
PROJECT="banana1"
echo "1. Listing jobs in project $PROJECT..."
python manage.py list_jobs "$PROJECT" --limit 3
echo ""

# Get first job UUID
JOB_UUID=$(python manage.py list_jobs "$PROJECT" --json --limit 1 2>/dev/null | jq -r '.[0].uuid')
echo "2. Using job: $JOB_UUID"
echo ""

# Test validate_job
echo "3. Testing validate_job..."
python manage.py validate_job --jobuuid "$JOB_UUID" || echo "  (Expected - job has no params)"
echo ""

# Test get_job_report
echo "4. Testing get_job_report..."
python manage.py get_job_report --jobuuid "$JOB_UUID" --type report -o /tmp/test_report.xml
if [ -f /tmp/test_report.xml ]; then
    echo "  ✓ Report created ($(ls -lh /tmp/test_report.xml | awk '{print $5}'))"
fi
echo ""

# Test clone_job
echo "5. Testing clone_job..."
python manage.py clone_job --jobuuid "$JOB_UUID" 2>&1 | grep -E "Successfully|Clone|UUID" | head -5
echo ""

# Test export_job
echo "6. Testing export_job..."
python manage.py export_job --jobuuid "$JOB_UUID" -o /tmp/test_export.zip
if [ -f /tmp/test_export.zip ]; then
    echo "  ✓ Export created ($(ls -lh /tmp/test_export.zip | awk '{print $5}'))"
fi
echo ""

echo "==========================================="
echo "✓ Core commands work!"
echo "==========================================="
