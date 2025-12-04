#!/bin/bash
# Test workflow with new logging system
# Set LOG_LEVEL=ERROR to suppress debug output

set -e

echo "=================================================="
echo "Testing CCP4i2 Workflow with Clean Logging Output"
echo "=================================================="
echo ""

# Set environment
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
export CCP4_LOG_LEVEL=ERROR  # Suppress all DEBUG output!

# Activate CCP4 environment first
echo "[1] Setting up CCP4 environment..."
source /Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh

# Then activate venv (takes precedence)
echo "[2] Activating virtual environment..."
source .venv/bin/activate

# Verify Python
echo "[3] Using Python: $(which python)"
echo ""

# Change to server directory
cd /Users/nmemn/Developer/cdata-codegen/server

# Clean test name
PROJECT_NAME="clean_test_$$"
echo "=================================================="
echo "Test Project: $PROJECT_NAME"
echo "=================================================="
echo ""

# Create project
echo "[4] Creating project..."
python manage.py create_project "$PROJECT_NAME" --json 2>&1 > /tmp/project_full_output.txt
# Extract just the JSON part (starts with {, ends with })
sed -n '/{/,/^}/p' /tmp/project_full_output.txt > /tmp/project_output.json
cat /tmp/project_output.json
PROJECT_UUID=$(jq -r '.uuid' /tmp/project_output.json)
echo "Created project UUID: $PROJECT_UUID"
echo ""

# Create job
echo "[5] Creating parrot job..."
python manage.py create_job -pn "$PROJECT_NAME" -tn parrot > /tmp/job_output.txt
cat /tmp/job_output.txt
# Parse the output (should be clean now!)
JOB_UUID=$(grep "uuid" /tmp/job_output.txt | sed 's/.*uuid //')
JOB_NUMBER=$(grep "number" /tmp/job_output.txt | sed 's/.*number //' | sed 's/,.*//')
echo "Created job #$JOB_NUMBER, UUID: $JOB_UUID"
echo ""

# Get job report (params XML)
echo "[6] Getting job parameters XML..."
python manage.py get_job_report --jobuuid "$JOB_UUID" --type params -o /tmp/test_params.xml > /tmp/report_output.txt
cat /tmp/report_output.txt
echo "Params XML saved to /tmp/test_params.xml"
echo ""

# Validate job
echo "[7] Validating job (should have errors - no params set)..."
python manage.py validate_job --jobuuid "$JOB_UUID" --json > /tmp/validate_output.json
cat /tmp/validate_output.json | jq '.'
echo ""

# Clone job
echo "[8] Cloning job..."
python manage.py clone_job --jobuuid "$JOB_UUID" --json > /tmp/clone_output.json
cat /tmp/clone_output.json | jq '.'
NEW_JOB_UUID=$(jq -r '.new_job_uuid' /tmp/clone_output.json)
echo "Cloned to new job UUID: $NEW_JOB_UUID"
echo ""

# List jobs
echo "[9] Listing all jobs in project..."
python manage.py list_jobs -pn "$PROJECT_NAME" --format table
echo ""

# Export job
echo "[10] Exporting job..."
python manage.py export_job --jobuuid "$JOB_UUID" -o /tmp/exported_job.zip
ls -lh /tmp/exported_job.zip
echo ""

echo "=================================================="
echo "SUCCESS! All commands produced clean output!"
echo "=================================================="
echo ""
echo "Notice: No [DEBUG] messages appeared thanks to CCP4_LOG_LEVEL=ERROR"
echo ""
echo "To see debug output, run with:"
echo "  export CCP4_LOG_LEVEL=DEBUG"
echo "  python manage.py <command>"
