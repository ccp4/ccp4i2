#!/bin/bash
# Wrapper script for filepicker that sets up the environment automatically.
#
# Usage:
#   ./run_filepicker.sh fedid@nx.diamond.ac.uk:/dls/.../model_building info
#   ./run_filepicker.sh fedid@nx.diamond.ac.uk:/dls/.../model_building pull -c x0001-x0050 -d ./local_data
#   ./run_filepicker.sh fedid@nx.diamond.ac.uk:/dls/.../model_building resolve -c x0053 -v
#
# CCP4 detection (in order):
#   1. Already sourced (CCP4 env var set)
#   2. CCP4_ROOT in .env file (this dir or project root)
#   3. Sibling directory of project root (../ccp4-*)

set -e  # Exit on error

# Get the directory where this script is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Project root is Docker/contrib/file_picker_project -> ../../.. -> ccp4i2
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Determine CCP4 root from environment or .env file
if [ -z "$CCP4" ]; then
    # Check for .env file in this directory or project root
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
        echo "  2. Set CCP4_ROOT in .env file (in this dir or project root)"
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

# Set up PYTHONPATH for ccp4i2 imports and filepicker module
export PYTHONPATH="$SCRIPT_DIR:$PROJECT_ROOT:$PROJECT_ROOT/server:$PYTHONPATH"

# Run filepicker
exec ccp4-python -m filepicker.cli "$@"
