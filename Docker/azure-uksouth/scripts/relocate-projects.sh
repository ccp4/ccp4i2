#!/bin/bash

# Relocate CCP4i2 Projects from /mnt/projects/* to /mnt/projects/CCP4I2_PROJECTS/*
#
# This script is designed to run INSIDE a server container instance to fix
# projects that were incorrectly migrated to the root of /mnt/projects.
#
# Usage: ./relocate-projects.sh [check|move]
#   check - Show what would be moved (dry run)
#   move  - Perform the actual relocation

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Source and destination paths
PROJECTS_ROOT="/mnt/projects"
DEST_SUBDIR="CCP4I2_PROJECTS"
DEST_PATH="$PROJECTS_ROOT/$DEST_SUBDIR"

# Directories to exclude from relocation (these are not project directories)
EXCLUDE_DIRS=("fixtures" "$DEST_SUBDIR")

# Check if a directory should be excluded
should_exclude() {
    local dir_name="$1"
    for exclude in "${EXCLUDE_DIRS[@]}"; do
        if [ "$dir_name" == "$exclude" ]; then
            return 0  # true, should exclude
        fi
    done
    return 1  # false, should not exclude
}

# Check what would be moved
check_relocations() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Checking Project Relocations (Dry Run)${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    if [ ! -d "$PROJECTS_ROOT" ]; then
        echo -e "${RED}Error: $PROJECTS_ROOT does not exist${NC}"
        exit 1
    fi

    echo -e "${BLUE}Source: $PROJECTS_ROOT${NC}"
    echo -e "${BLUE}Destination: $DEST_PATH${NC}"
    echo ""

    # Count directories to move
    local count=0
    local total_size=0

    echo -e "${YELLOW}Directories to be relocated:${NC}"
    echo ""

    for dir in "$PROJECTS_ROOT"/*/; do
        # Skip if not a directory
        [ -d "$dir" ] || continue

        # Get just the directory name
        dir_name=$(basename "$dir")

        # Skip excluded directories
        if should_exclude "$dir_name"; then
            echo -e "${BLUE}  [SKIP] $dir_name (excluded)${NC}"
            continue
        fi

        # Check if already exists in destination
        if [ -d "$DEST_PATH/$dir_name" ]; then
            echo -e "${RED}  [CLASH] $dir_name (already exists in destination)${NC}"
            continue
        fi

        # Get size
        size=$(du -sh "$dir" 2>/dev/null | cut -f1)
        echo -e "${GREEN}  [MOVE] $dir_name ($size)${NC}"
        ((count++))
    done

    echo ""
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Summary${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""
    echo "Directories to relocate: $count"
    echo "Excluded directories: ${EXCLUDE_DIRS[*]}"
    echo ""

    if [ $count -gt 0 ]; then
        echo -e "${YELLOW}Run '$0 move' to perform the relocation${NC}"
    else
        echo -e "${GREEN}No directories need to be relocated${NC}"
    fi
}

# Perform the relocation
do_move() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Relocating Project Directories${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    if [ ! -d "$PROJECTS_ROOT" ]; then
        echo -e "${RED}Error: $PROJECTS_ROOT does not exist${NC}"
        exit 1
    fi

    # Create destination directory if it doesn't exist
    if [ ! -d "$DEST_PATH" ]; then
        echo -e "${YELLOW}Creating destination directory: $DEST_PATH${NC}"
        mkdir -p "$DEST_PATH"
    fi

    local count=0
    local errors=0

    for dir in "$PROJECTS_ROOT"/*/; do
        # Skip if not a directory
        [ -d "$dir" ] || continue

        # Get just the directory name
        dir_name=$(basename "$dir")

        # Skip excluded directories
        if should_exclude "$dir_name"; then
            echo -e "${BLUE}Skipping: $dir_name (excluded)${NC}"
            continue
        fi

        # Check if already exists in destination
        if [ -d "$DEST_PATH/$dir_name" ]; then
            echo -e "${RED}Error: $dir_name already exists in destination - skipping${NC}"
            ((errors++))
            continue
        fi

        # Move the directory
        echo -e "${YELLOW}Moving: $dir_name${NC}"
        if mv "$dir" "$DEST_PATH/$dir_name"; then
            echo -e "${GREEN}  Success: $dir_name${NC}"
            ((count++))
        else
            echo -e "${RED}  Failed: $dir_name${NC}"
            ((errors++))
        fi
    done

    echo ""
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Relocation Complete${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""
    echo "Directories moved: $count"
    echo "Errors: $errors"
    echo ""

    if [ $errors -eq 0 ]; then
        echo -e "${GREEN}All relocations completed successfully!${NC}"
    else
        echo -e "${RED}Some relocations failed - please check the errors above${NC}"
    fi

    echo ""
    echo -e "${YELLOW}Next steps:${NC}"
    echo "1. Run the database patch script to update Project.directory paths"
    echo "2. Restart the server container to pick up changes"
}

# Show usage
show_usage() {
    echo "CCP4i2 Project Directory Relocation Script"
    echo ""
    echo "Relocates project directories from $PROJECTS_ROOT/* to $DEST_PATH/*"
    echo ""
    echo "This script should be run INSIDE a server container instance."
    echo ""
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  check   Show what would be moved (dry run)"
    echo "  move    Perform the actual relocation"
    echo ""
    echo "Excluded directories: ${EXCLUDE_DIRS[*]}"
    echo ""
}

# Main command handler
case "${1:-}" in
    check)
        check_relocations
        ;;
    move)
        do_move
        ;;
    help|--help|-h)
        show_usage
        ;;
    *)
        if [ -n "$1" ]; then
            echo -e "${RED}Unknown command: $1${NC}"
            echo ""
        fi
        show_usage
        exit 1
        ;;
esac
