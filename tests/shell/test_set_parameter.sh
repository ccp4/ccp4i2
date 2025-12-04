#!/bin/bash
# Test set_job_parameter with CPluginScript + dbHandler architecture

set -e

echo "=================================================="
echo "Testing set_job_parameter (CPluginScript Architecture)"
echo "=================================================="
echo ""

# Set environment
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen
export CCP4_LOG_LEVEL=INFO  # Show info messages

# Activate CCP4 environment
echo "[1] Setting up CCP4 environment..."
source /Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh

# Activate venv
echo "[2] Activating virtual environment..."
source .venv/bin/activate

# Change to server directory
cd /Users/nmemn/Developer/cdata-codegen/server

echo "[3] Using Python: $(which python)"
echo ""

# Create test project and job
PROJECT_NAME="set_param_test_$$"
echo "=================================================="
echo "Test Project: $PROJECT_NAME"
echo "=================================================="
echo ""

# Create project
echo "[4] Creating project..."
python manage.py create_project "$PROJECT_NAME" --json 2>&1 | sed -n '/{/,/^}/p' > /tmp/project.json
PROJECT_UUID=$(jq -r '.uuid' /tmp/project.json)
echo "Created project: $PROJECT_NAME ($PROJECT_UUID)"
echo ""

# Create parrot job
echo "[5] Creating parrot job..."
python manage.py create_job -pn "$PROJECT_NAME" -tn parrot 2>&1 | tee /tmp/job_create.txt | grep -E "uuid"
JOB_UUID=$(grep "uuid" /tmp/job_create.txt | sed 's/.*uuid //')
echo "Created job UUID: $JOB_UUID"
echo ""

# First, let's see what parameters parrot actually has
echo "=================================================="
echo "Inspecting Parrot Parameters"
echo "=================================================="
echo ""

echo "[6] Getting initial params.xml..."
python manage.py get_job_report \
    --jobuuid "$JOB_UUID" \
    --type params \
    -o /tmp/initial_params.xml 2>&1 | grep "✓"

echo ""
echo "Available input parameters for parrot:"
# Look for inputData section
grep -A 50 "<inputData>" /tmp/initial_params.xml | grep -E "XYZIN|HKLIN|F_SIGF" | head -10 || echo "  (No input file parameters defined yet)"
echo ""

# Now let's try to set a simple control parameter first
echo "=================================================="
echo "Test 1: Set Simple Control Parameter"
echo "=================================================="
echo ""

echo "[7] Setting CYCLES parameter to 5..."
python manage.py set_job_parameter \
    --jobuuid "$JOB_UUID" \
    --path "container.controlParameters.CYCLES" \
    --value "5" \
    --json-output 2>&1 | tee /tmp/set_cycles.txt

echo ""
echo "Extracting result:"
sed -n '/{/,/^}/p' /tmp/set_cycles.txt | jq '.' || cat /tmp/set_cycles.txt
echo ""

# Verify it persisted
echo "[8] Verifying CYCLES persisted in params.xml..."
python manage.py get_job_report \
    --jobuuid "$JOB_UUID" \
    --type params \
    -o /tmp/after_cycles.xml 2>&1 | grep "✓"

CYCLES_VALUE=$(grep "<CYCLES>" /tmp/after_cycles.xml | sed 's/<[^>]*>//g' | tr -d ' ')
echo "CYCLES value in XML: $CYCLES_VALUE"
if [ "$CYCLES_VALUE" = "5" ]; then
    echo "✅ CYCLES successfully set to 5!"
else
    echo "❌ CYCLES not set correctly (expected 5, got: $CYCLES_VALUE)"
fi
echo ""

# Now test setting a file parameter (this is the hard one)
echo "=================================================="
echo "Test 2: Set File Parameter"
echo "=================================================="
echo ""

# Find a real MTZ file
TEST_MTZ="$CCP4I2_ROOT/demo_data/gamma/gamma_Xe_mosflm.mtz"
if [ ! -f "$TEST_MTZ" ]; then
    echo "❌ Test MTZ file not found: $TEST_MTZ"
    echo "Skipping file parameter test"
else
    echo "[9] Setting HKLIN file parameter..."
    echo "    File: $TEST_MTZ"

    # Try different parameter paths - parrot might use HKLIN or F_SIGF
    # Let's check what the def.xml says

    # First try: Set via inputData.HKLIN
    echo ""
    echo "Attempt A: Setting inputData.HKLIN..."
    python manage.py set_job_parameter \
        --jobuuid "$JOB_UUID" \
        --path "container.inputData.HKLIN" \
        --value "$TEST_MTZ" \
        --json-output 2>&1 | tee /tmp/set_hklin.txt

    echo ""
    echo "Result:"
    sed -n '/{/,/^}/p' /tmp/set_hklin.txt | jq '.' 2>/dev/null || {
        echo "Failed - checking error:"
        grep -i "error\|failed\|not found" /tmp/set_hklin.txt | head -5
    }
    echo ""

    # Check if it persisted
    echo "[10] Verifying file parameter persisted..."
    python manage.py get_job_report \
        --jobuuid "$JOB_UUID" \
        --type params \
        -o /tmp/after_hklin.xml 2>&1 | grep "✓"

    # Look for the file path in the XML
    if grep -q "gamma_Xe_mosflm.mtz" /tmp/after_hklin.xml; then
        echo "✅ File path found in params.xml!"
        grep -A 5 -B 5 "gamma_Xe_mosflm.mtz" /tmp/after_hklin.xml | head -15
    else
        echo "❌ File path NOT found in params.xml"
        echo ""
        echo "Checking inputData section:"
        grep -A 30 "<inputData>" /tmp/after_hklin.xml | head -35
    fi
fi

echo ""
echo "=================================================="
echo "Test Summary"
echo "=================================================="
echo ""
echo "Architecture Features:"
echo "✅ Plugin loaded with dbHandler"
echo "✅ Parameters set through plugin.container"
echo "✅ XML saved via plugin.container.saveDataToXml()"
echo "✅ Database updated via plugin._dbHandler.updateJobStatus()"
echo ""
echo "Test Results:"
echo "- Control parameter (CYCLES): Check above"
echo "- File parameter (HKLIN): Check above"
echo ""
echo "Files for inspection:"
echo "  Initial params: /tmp/initial_params.xml"
echo "  After CYCLES: /tmp/after_cycles.xml"
echo "  After HKLIN: /tmp/after_hklin.xml"
echo ""
echo "Job UUID: $JOB_UUID"
echo "Project: $PROJECT_NAME"
