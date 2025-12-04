#!/bin/bash
# Test CPluginScript + dbHandler architecture for parameter setting

set -e

echo "=================================================="
echo "Testing CPluginScript + dbHandler Architecture"
echo "=================================================="
echo ""

# Set environment
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
export CCP4_LOG_LEVEL=INFO  # Show info messages to see what's happening

# Activate CCP4 environment
echo "[1] Setting up CCP4 environment..."
source /Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh

# Activate venv
echo "[2] Activating virtual environment..."
source .venv/bin/activate

# Change to server directory
cd /Users/nmemn/Developer/cdata-codegen/server

# Verify Python
echo "[3] Using Python: $(which python)"
echo ""

# Create test project and job
PROJECT_NAME="plugin_test_$$"
echo "=================================================="
echo "Test Project: $PROJECT_NAME"
echo "=================================================="
echo ""

# Create project
echo "[4] Creating project..."
python manage.py create_project "$PROJECT_NAME" --json 2>&1 | sed -n '/{/,/^}/p' > /tmp/project.json
PROJECT_UUID=$(jq -r '.uuid' /tmp/project.json)
echo "Created project UUID: $PROJECT_UUID"
echo ""

# Create parrot job
echo "[5] Creating parrot job..."
python manage.py create_job -pn "$PROJECT_NAME" -tn parrot 2>&1 | grep -E "uuid|number" | tee /tmp/job_output.txt
JOB_UUID=$(grep "uuid" /tmp/job_output.txt | sed 's/.*uuid //')
echo "Created job UUID: $JOB_UUID"
echo ""

# Now test parameter setting with the NEW architecture
echo "=================================================="
echo "Testing NEW CPluginScript Architecture"
echo "=================================================="
echo ""

# Test 1: Set a file parameter
echo "[6] Setting HKLIN file parameter..."
TEST_FILE="$CCP4I2_ROOT/demo_data/gamma/gamma_Xe_mosflm.mtz"
if [ -f "$TEST_FILE" ]; then
    echo "Using test file: $TEST_FILE"
    python manage.py set_job_parameter \
        --jobuuid "$JOB_UUID" \
        --path "inputData.F_SIGF" \
        --value "$TEST_FILE" \
        --json-output 2>&1 | tee /tmp/set_param_result.json | sed -n '/{/,/^}/p'
    echo ""
    echo "Result:"
    cat /tmp/set_param_result.json | jq '.'
    echo ""
else
    echo "Test file not found: $TEST_FILE"
    echo "Skipping file parameter test"
fi

# Test 2: Get params XML to verify it persisted
echo "[7] Verifying parameter persisted in params.xml..."
python manage.py get_job_report \
    --jobuuid "$JOB_UUID" \
    --type params \
    -o /tmp/test_params.xml 2>&1 | grep -v "^Using\|^Debug\|^status:"

if [ -f /tmp/test_params.xml ]; then
    echo ""
    echo "Checking params.xml for F_SIGF..."
    if grep -q "F_SIGF" /tmp/test_params.xml; then
        echo "✅ F_SIGF found in params.xml"
        # Show the relevant section
        grep -A 5 "F_SIGF" /tmp/test_params.xml | head -10
    else
        echo "❌ F_SIGF NOT found in params.xml"
    fi
fi

echo ""
echo "=================================================="
echo "Architecture Test Complete!"
echo "=================================================="
echo ""
echo "Key points to verify:"
echo "1. ✅ Plugin loaded with dbHandler attached"
echo "2. ✅ Parameter set through plugin.container"
echo "3. ✅ File path should show in result"
echo "4. ✅ dbFileId should be present (if DB integration works)"
echo "5. ✅ Parameter should persist in params.xml"
echo ""
echo "Check the output above to confirm all features are working!"
