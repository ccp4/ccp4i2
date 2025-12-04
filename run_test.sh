#!/bin/bash
# Test dispatcher for i2run tests
# Usage: ./run_test.sh <test_file_or_dir> [pytest_args...]
# Example: ./run_test.sh server/ccp4x/tests/i2run/test_i2run.py
# Example: ./run_test.sh server/ccp4x/tests/i2run/test_i2run.py::test_prosmart_refmac
# Example: ./run_test.sh i2run/ --ignore=test_mrbump.py -n auto
# Example: ./run_test.sh i2run/ -n 4 -v

# Load environment configuration from .env file
if [ -f .env ]; then
    echo "Loading environment from .env..."
    export $(grep -v '^#' .env | xargs)
    echo "Using CCP4: $CCP4_VERSION (Python $PYTHON_VERSION)"
    echo "CCP4 Root: $CCP4_ROOT"
    echo "Virtual env: $VENV_DIR"
else
    echo "WARNING: .env file not found, using defaults"
    export CCP4_ROOT=/Users/nmemn/Developer/ccp4-20251105
    export VENV_DIR=.venv
    export CCP4_VERSION=ccp4-20251105
    export PYTHON_VERSION=3.11
    export DJANGO_SETTINGS_MODULE=ccp4x.config.settings
    export CCP4I2_ROOT=$(pwd)
fi

# Source CCP4 setup script
if [ -f "$CCP4_ROOT/bin/ccp4.setup-sh" ]; then
    source "$CCP4_ROOT/bin/ccp4.setup-sh"
else
    echo "ERROR: CCP4 setup script not found at $CCP4_ROOT/bin/ccp4.setup-sh"
    exit 1
fi

set -e  # Exit on error

# Set up environment
export CCP4I2_ROOT=$(pwd)
export PYTHONPATH=$(pwd):$(pwd)/server:$PYTHONPATH
# Use test settings for pytest, production settings otherwise
export DJANGO_SETTINGS_MODULE="${DJANGO_SETTINGS_MODULE:-ccp4x.config.test_settings}"

# Determine Python interpreter
# Priority 1: ccp4-python (after sourcing ccp4.setup-sh)
# Priority 2: .venv virtual environment (development fallback)
if command -v ccp4-python &> /dev/null; then
    echo "Using ccp4-python: $(which ccp4-python)"
    PYTHON_CMD="ccp4-python"
elif [ -L .venv ] || [ -d .venv ]; then
    echo "Activating virtual environment: .venv -> $(readlink .venv 2>/dev/null || echo 'direct')"
    source .venv/bin/activate
    PYTHON_CMD="python"
else
    echo "ERROR: No Python interpreter found"
    echo "Either source CCP4 setup (for ccp4-python) or create .venv"
    exit 1
fi

# Get test path from first argument
TEST_PATH="${1:-server/ccp4x/tests/i2run/test_i2run.py}"
shift  # Remove first argument, leaving remaining args for pytest

# Run pytest with all remaining arguments passed through
echo ""
echo "Running pytest $TEST_PATH $@"
echo "Using Python: $(which $PYTHON_CMD)"
$PYTHON_CMD -m pytest "$TEST_PATH" "$@"
