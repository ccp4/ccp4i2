#!/bin/bash
# Common setup for shell-based tests
# Source this at the start of each test script:
#   source "$(dirname "$0")/common.sh"

# Determine project root from script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../.." && pwd)"

# Set CCP4I2_ROOT if not already set
export CCP4I2_ROOT="${CCP4I2_ROOT:-$PROJECT_ROOT}"

# Find and source CCP4 environment if not already sourced
if [ -z "$CCP4" ]; then
    # Try to find CCP4 in sibling directories
    for dir in "$PROJECT_ROOT"/../ccp4-*/bin/ccp4.setup-sh; do
        if [ -f "$dir" ]; then
            CCP4_ROOT="$(dirname "$(dirname "$dir")")"
            break
        fi
    done

    if [ -n "$CCP4_ROOT" ] && [ -f "$CCP4_ROOT/bin/ccp4.setup-sh" ]; then
        echo "Sourcing CCP4 from: $CCP4_ROOT"
        source "$CCP4_ROOT/bin/ccp4.setup-sh"
    else
        echo "ERROR: CCP4 not found. Either:"
        echo "  1. Source CCP4 setup: source /path/to/ccp4/bin/ccp4.setup-sh"
        echo "  2. Place CCP4 in a sibling directory (../ccp4-*)"
        exit 1
    fi
fi

# Activate virtual environment if present
if [ -d "$CCP4I2_ROOT/.venv" ]; then
    source "$CCP4I2_ROOT/.venv/bin/activate"
fi

# Verify environment
echo "Environment:"
echo "  CCP4I2_ROOT: $CCP4I2_ROOT"
echo "  CCP4: $CCP4"
echo "  Python: $(which python 2>/dev/null || which ccp4-python 2>/dev/null)"
echo ""
