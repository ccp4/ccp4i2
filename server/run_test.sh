#!/bin/bash
# Test runner for ccp4i2 tests
#
# Prerequisites:
#   1. CCP4 installed with ccp4-python
#   2. ccp4i2 installed in editable mode: ccp4-python -m pip install -e ".[full]"
#
# Usage: ./run_test.sh <test_file_or_dir> [pytest_args...]
#
# Examples:
#   ./run_test.sh ccp4i2/tests/i2run/test_parrot.py -v
#   ./run_test.sh ccp4i2/tests/i2run/test_servalcat.py::test_servalcat_basic
#   ./run_test.sh ccp4i2/tests/i2run/ -n 4
#   ./run_test.sh ccp4i2/tests/i2run/ --ignore=test_mrbump.py -n auto

set -e  # Exit on error

# Get the directory where this script is located (server/)
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Project root is the parent of server/
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Determine CCP4 root from environment or .env file
if [ -z "$CCP4" ]; then
    # Check for .env file in server/ directory or project root
    if [ -f "$SCRIPT_DIR/.env" ]; then
        export $(grep -v '^#' "$SCRIPT_DIR/.env" | grep CCP4_ROOT | xargs)
        CCP4_ROOT="${CCP4_ROOT:-}"
    elif [ -f "$PROJECT_ROOT/.env" ]; then
        export $(grep -v '^#' "$PROJECT_ROOT/.env" | grep CCP4_ROOT | xargs)
        CCP4_ROOT="${CCP4_ROOT:-}"
    fi

    # Try to find CCP4 in sibling directories of project root
    # e.g., if project is at ~/Developer/ccp4i2, look for ~/Developer/ccp4-*
    if [ -z "$CCP4_ROOT" ]; then
        for dir in "$PROJECT_ROOT"/../ccp4-*/bin/ccp4.setup-sh; do
            if [ -f "$dir" ]; then
                CCP4_ROOT="$(dirname $(dirname $dir))"
                break
            fi
        done
    fi

    if [ -z "$CCP4_ROOT" ]; then
        echo "ERROR: CCP4 not found. Either:"
        echo "  1. Source CCP4 setup: source /path/to/ccp4/bin/ccp4.setup-sh"
        echo "  2. Set CCP4_ROOT in .env file (in server/ or project root)"
        echo "  3. Place CCP4 in a sibling directory of the project (../ccp4-*)"
        exit 1
    fi

    echo "Sourcing CCP4 from: $CCP4_ROOT"
    source "$CCP4_ROOT/bin/ccp4.setup-sh"
fi

# Verify ccp4-python is available
if ! command -v ccp4-python &> /dev/null; then
    echo "ERROR: ccp4-python not found in PATH"
    echo "Ensure CCP4 environment is sourced"
    exit 1
fi

# Set environment for tests
# CCP4I2_ROOT should point to the server/ directory where ccp4i2 package lives
export CCP4I2_ROOT="$SCRIPT_DIR"
export DJANGO_SETTINGS_MODULE="${DJANGO_SETTINGS_MODULE:-ccp4i2.config.test_settings}"

# Change to server directory so relative paths in tests work correctly
cd "$SCRIPT_DIR"

# Verify ccp4i2 is installed
if ! ccp4-python -c "import ccp4i2" 2>/dev/null; then
    echo "ERROR: ccp4i2 not installed in ccp4-python"
    echo "Run: ccp4-python -m pip install -e \".[full]\""
    exit 1
fi

# Get test path from first argument
TEST_PATH="${1:-ccp4i2/tests/i2run/}"
shift 2>/dev/null || true  # Remove first argument if present

echo ""
echo "CCP4I2_ROOT: $CCP4I2_ROOT"
echo "DJANGO_SETTINGS_MODULE: $DJANGO_SETTINGS_MODULE"
echo "Python: $(which ccp4-python)"
echo ""
echo "Running: pytest $TEST_PATH $@"
echo ""

ccp4-python -m pytest "$TEST_PATH" "$@"

echo ""
echo "Test projects are stored in: ~/.cache/ccp4i2-tests/"
echo "To clean up: rm -rf ~/.cache/ccp4i2-tests/*"
