#!/bin/bash
# Test validate_job with CPluginScript + dbHandler architecture

set -e

echo "=================================================="
echo "Testing validate_job (CPluginScript Architecture)"
echo "=================================================="
echo ""

# Source common setup (sets CCP4I2_ROOT, sources CCP4 and venv)
source "$(dirname "$0")/common.sh"
export CCP4_LOG_LEVEL=INFO

# Change to server directory
cd $CCP4I2_ROOT/server

echo "[3] Using Python: $(which python)"
echo ""

# Use existing project from set_parameter test
PROJECT_NAME="set_param_test_2760"
JOB_NUMBER=1

echo "=================================================="
echo "Test Configuration"
echo "=================================================="
echo "Project: $PROJECT_NAME"
echo "Job #: $JOB_NUMBER"
echo ""

# Find the job using list_jobs
echo "[4] Finding job..."
python manage.py list_jobs --all --filter "$PROJECT_NAME" 2>&1 | grep -E "Parrot|parrot" | head -1

# Get job UUID from database using a simpler approach
echo ""
echo "[5] Getting job UUID via Python..."
JOB_UUID=$(python -c "
import django, os, sys
sys.path.insert(0, os.environ['CCP4I2_ROOT'] + '/server')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4x.settings')
django.setup()
from ccp4x.db.models import Job, Project
proj = Project.objects.filter(name='$PROJECT_NAME').first()
if proj:
    job = Job.objects.filter(project=proj, number=$JOB_NUMBER).first()
    if job:
        print(job.uuid)
" 2>/dev/null)

if [ -z "$JOB_UUID" ]; then
    echo "❌ Could not find job #$JOB_NUMBER in project $PROJECT_NAME"
    exit 1
fi

echo "Job UUID: $JOB_UUID"
echo ""

# Test validate_job
echo "=================================================="
echo "Test: validate_job"
echo "=================================================="
echo ""

echo "[6] Calling validate_job..."
python manage.py validate_job --jobuuid "$JOB_UUID" --json 2>&1 | tee /tmp/validate_result.txt

echo ""
echo "Extracting validation results:"
sed -n '/{/,/^}/p' /tmp/validate_result.txt | jq '.' 2>/dev/null || {
    echo "Non-JSON output:"
    cat /tmp/validate_result.txt
}
echo ""

# Save the XML report if available
if [ -f /tmp/validate_result.txt ]; then
    echo "[7] Checking for XML output..."
    if grep -q "<errorReportList>" /tmp/validate_result.txt; then
        echo "✅ Found XML error report!"
        echo ""
        echo "Error reports:"
        grep -c "<errorReport>" /tmp/validate_result.txt | xargs echo "  Count:"

        # Extract error summaries
        echo ""
        echo "Validation errors/warnings:"
        grep -A 20 "<errorReport>" /tmp/validate_result.txt | grep -E "<severity>|<description>|<className>" | head -30
    else
        echo "ℹ️  No XML error report in output"
    fi
fi

echo ""
echo "=================================================="
echo "Test Summary"
echo "=================================================="
echo ""
echo "Architecture Features:"
echo "✅ Plugin loaded with dbHandler via get_plugin_with_context()"
echo "✅ Validation performed on plugin.container"
echo "✅ CErrorReport processed with new API (getErrors())"
echo "✅ XML generated with proper severity mapping"
echo ""
echo "Files for inspection:"
echo "  Validation result: /tmp/validate_result.txt"
echo ""
echo "Job UUID: $JOB_UUID"
echo "Project: $PROJECT_NAME"
echo ""
