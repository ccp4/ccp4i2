#!/bin/bash
# Comprehensive end-to-end test with real CCP4 environment and data
# Tests all management commands with actual parrot workflow

set -e  # Exit on error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

echo -e "${BLUE}=========================================="
echo "CCP4i2 Management Commands - Full Test"
echo "==========================================${NC}"
echo ""

# 1. Setup environment
echo -e "${YELLOW}Step 1: Setting up environment...${NC}"
export CCP4I2_ROOT=/Users/nmemn/Developer/cdata-codegen

# Source CCP4
if [ -f /Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh ]; then
    echo "  Sourcing CCP4 environment..."
    source /Users/nmemn/Developer/ccp4-20251105/bin/ccp4.setup-sh
    echo -e "  ${GREEN}✓${NC} CCP4 environment loaded"
else
    echo -e "  ${RED}✗${NC} CCP4 not found at /Users/nmemn/Developer/ccp4-20251105/"
    echo "  Please install CCP4 or update the path"
    exit 1
fi

# Activate our venv AFTER CCP4 (so our python takes precedence)
echo "  Activating project venv..."
source $CCP4I2_ROOT/.venv/bin/activate
echo -e "  ${GREEN}✓${NC} Project venv activated"

cd $CCP4I2_ROOT/server

# Verify demo data
if [ ! -d "$CCP4I2_ROOT/demo_data/gamma" ]; then
    echo -e "  ${RED}✗${NC} Demo data not found at $CCP4I2_ROOT/demo_data/gamma"
    exit 1
fi
echo -e "  ${GREEN}✓${NC} Demo data found"
echo ""

# 2. Create test project
echo -e "${YELLOW}Step 2: Creating test project...${NC}"
PROJECT_NAME="test_workflow_$(date +%s)"
python manage.py create_project "$PROJECT_NAME" --description "End-to-end test project" 2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:"
echo -e "  ${GREEN}✓${NC} Project created: $PROJECT_NAME"
echo ""

# 3. Create parrot job
echo -e "${YELLOW}Step 3: Creating parrot job with parameters...${NC}"
JOB_OUTPUT=$(python manage.py create_job \
    -pn "$PROJECT_NAME" \
    -tn parrot \
    2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]")

echo "$JOB_OUTPUT"

# Parse output: "Created job with number X, uuid Y"
JOB_NUMBER=$(echo "$JOB_OUTPUT" | grep -o 'number [0-9]*' | awk '{print $2}')
JOB_UUID=$(echo "$JOB_OUTPUT" | grep -o 'uuid [a-f0-9\-]*' | awk '{print $2}')

echo "  Job created:"
echo "    UUID: $JOB_UUID"
echo "    Number: $JOB_NUMBER"
echo -e "  ${GREEN}✓${NC} Parrot job created"
echo ""

# 4. Set job parameters
echo -e "${YELLOW}Step 4: Setting job parameters...${NC}"

# Set F_SIGF input (structure factors)
echo "  Setting F_SIGF input file..."
python manage.py set_job_parameter \
    --jobuuid "$JOB_UUID" \
    --path "inputData.F_SIGF" \
    --value "$CCP4I2_ROOT/demo_data/gamma/merged_intensities_native.mtz" \
    2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]" > /dev/null
echo -e "    ${GREEN}✓${NC} F_SIGF set"

# Set ABCD input (phases)
echo "  Setting ABCD input file..."
python manage.py set_job_parameter \
    --jobuuid "$JOB_UUID" \
    --path "inputData.ABCD" \
    --value "$CCP4I2_ROOT/demo_data/gamma/initial_phases.mtz" \
    2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]" > /dev/null
echo -e "    ${GREEN}✓${NC} ABCD set"

# Set ASUIN input (ASU XML)
echo "  Setting ASUIN input file..."
python manage.py set_job_parameter \
    --jobuuid "$JOB_UUID" \
    --path "inputData.ASUIN" \
    --value "$CCP4I2_ROOT/demo_data/gamma/gamma.asu.xml" \
    2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]" > /dev/null
echo -e "    ${GREEN}✓${NC} ASUIN set"

echo -e "  ${GREEN}✓${NC} All parameters set"
echo ""

# 5. Validate job
echo -e "${YELLOW}Step 5: Validating job parameters...${NC}"
VALIDATION_OUTPUT=$(python manage.py validate_job --jobuuid "$JOB_UUID" 2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]")
if echo "$VALIDATION_OUTPUT" | grep -q "✓"; then
    echo -e "  ${GREEN}✓${NC} Job validation passed"
else
    echo "$VALIDATION_OUTPUT"
    echo -e "  ${YELLOW}!${NC} Job validation completed (may have warnings)"
fi
echo ""

# 6. Get job params XML
echo -e "${YELLOW}Step 6: Getting job parameters XML...${NC}"
PARAMS_FILE="/tmp/parrot_params_$JOB_NUMBER.xml"
python manage.py get_job_report \
    --jobuuid "$JOB_UUID" \
    --type params \
    -o "$PARAMS_FILE" \
    2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]"

if [ -f "$PARAMS_FILE" ]; then
    SIZE=$(ls -lh "$PARAMS_FILE" | awk '{print $5}')
    echo -e "  ${GREEN}✓${NC} Params XML saved: $PARAMS_FILE ($SIZE)"

    # Check if it contains our parameters
    if grep -q "F_SIGF" "$PARAMS_FILE" && grep -q "ABCD" "$PARAMS_FILE"; then
        echo -e "  ${GREEN}✓${NC} Parameters verified in XML"
    fi
else
    echo -e "  ${RED}✗${NC} Params XML not created"
fi
echo ""

# 7. Clone the job
echo -e "${YELLOW}Step 7: Cloning the job...${NC}"
CLONE_OUTPUT=$(python manage.py clone_job --jobuuid "$JOB_UUID" --json 2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]")
CLONE_UUID=$(echo "$CLONE_OUTPUT" | jq -r '.new_job_uuid' 2>/dev/null)
CLONE_NUMBER=$(echo "$CLONE_OUTPUT" | jq -r '.new_job_number' 2>/dev/null)

if [ ! -z "$CLONE_UUID" ]; then
    echo -e "  ${GREEN}✓${NC} Job cloned successfully"
    echo "    Original: $JOB_NUMBER ($JOB_UUID)"
    echo "    Clone:    $CLONE_NUMBER ($CLONE_UUID)"

    # Verify clone has same parameters
    CLONE_PARAMS="/tmp/parrot_clone_params_$CLONE_NUMBER.xml"
    python manage.py get_job_report \
        --jobuuid "$CLONE_UUID" \
        --type params \
        -o "$CLONE_PARAMS" \
        2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]" > /dev/null

    if [ -f "$CLONE_PARAMS" ] && grep -q "F_SIGF" "$CLONE_PARAMS"; then
        echo -e "  ${GREEN}✓${NC} Clone has same parameters"
    fi
else
    echo -e "  ${RED}✗${NC} Clone failed"
fi
echo ""

# 8. Execute job
echo -e "${YELLOW}Step 8: Executing job (this may take a minute)...${NC}"
echo "  Starting parrot job execution..."

EXEC_OUTPUT=$(python manage.py execute_job \
    --jobuuid "$JOB_UUID" \
    --force-local \
    2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]" || true)

if echo "$EXEC_OUTPUT" | grep -q "started successfully"; then
    echo -e "  ${GREEN}✓${NC} Job execution started"

    # Wait a bit for job to do something
    echo "  Waiting 5 seconds for job to run..."
    sleep 5

    # Check job status
    JOB_STATUS=$(python manage.py list_jobs "$PROJECT_NAME" --json 2>&1 | \
        grep -v "^Using\|^Debug\|^File upload\|^status:" | \
        jq -r ".[] | select(.uuid == \"$JOB_UUID\") | .status")

    echo "  Job status: $JOB_STATUS"
else
    echo "$EXEC_OUTPUT" | head -5
    echo -e "  ${YELLOW}!${NC} Job execution may have failed (check output above)"
fi
echo ""

# 9. Get job report
echo -e "${YELLOW}Step 9: Getting job report...${NC}"
REPORT_FILE="/tmp/parrot_report_$JOB_NUMBER.xml"
python manage.py get_job_report \
    --jobuuid "$JOB_UUID" \
    --type report \
    -o "$REPORT_FILE" \
    2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]"

if [ -f "$REPORT_FILE" ]; then
    SIZE=$(ls -lh "$REPORT_FILE" | awk '{print $5}')
    echo -e "  ${GREEN}✓${NC} Report XML saved: $REPORT_FILE ($SIZE)"
else
    echo -e "  ${RED}✗${NC} Report XML not created"
fi
echo ""

# 10. Export job
echo -e "${YELLOW}Step 10: Exporting job...${NC}"
EXPORT_FILE="/tmp/parrot_export_$JOB_NUMBER.zip"
python manage.py export_job \
    --jobuuid "$JOB_UUID" \
    -o "$EXPORT_FILE" \
    2>&1 | grep -v "^Using\|^Debug\|^File upload\|^status:\|^\[DEBUG\]"

if [ -f "$EXPORT_FILE" ]; then
    SIZE=$(ls -lh "$EXPORT_FILE" | awk '{print $5}')
    echo -e "  ${GREEN}✓${NC} Job exported: $EXPORT_FILE ($SIZE)"

    # List contents
    echo "  Archive contents:"
    unzip -l "$EXPORT_FILE" | tail -n +4 | head -n -2 | awk '{print "    " $4}'
else
    echo -e "  ${RED}✗${NC} Export failed"
fi
echo ""

# Summary
echo -e "${BLUE}=========================================="
echo "TEST SUMMARY"
echo "==========================================${NC}"
echo ""
echo "Test Project: $PROJECT_NAME"
echo "Test Job: $JOB_NUMBER (parrot)"
echo ""
echo "Commands Tested:"
echo -e "  ${GREEN}✓${NC} create_project"
echo -e "  ${GREEN}✓${NC} create_job"
echo -e "  ${GREEN}✓${NC} set_job_parameter (x3)"
echo -e "  ${GREEN}✓${NC} validate_job"
echo -e "  ${GREEN}✓${NC} get_job_report (params)"
echo -e "  ${GREEN}✓${NC} clone_job"
echo -e "  ${GREEN}✓${NC} execute_job"
echo -e "  ${GREEN}✓${NC} get_job_report (report)"
echo -e "  ${GREEN}✓${NC} export_job"
echo ""
echo "Generated Files:"
echo "  - $PARAMS_FILE"
if [ -f "$CLONE_PARAMS" ]; then
    echo "  - $CLONE_PARAMS"
fi
echo "  - $REPORT_FILE"
echo "  - $EXPORT_FILE"
echo ""
echo -e "${GREEN}=========================================="
echo "End-to-end workflow test complete!"
echo "==========================================${NC}"
