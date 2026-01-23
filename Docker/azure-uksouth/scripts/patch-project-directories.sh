#!/bin/bash

# Patch CCP4i2 Project Directory Paths in Database
#
# This script updates the Project.directory field in the Django database
# to reflect the new storage mount paths.
#
# Old path pattern: /usr/src/app/mydata/CompoundDatabaseData/CCP4I2_PROJECTS/{project_name}
# New path pattern: /mnt/projects/CCP4I2_PROJECTS/{project_name}
#
# This script should be run INSIDE a server container instance.
#
# Usage: ./patch-project-directories.sh [check|patch]
#   check - Show what would be changed (dry run)
#   patch - Perform the actual database update

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

# Path patterns
OLD_PATH_PREFIX="/usr/src/app/mydata/CompoundDatabaseData/CCP4I2_PROJECTS"
NEW_PATH_PREFIX="/mnt/projects/CCP4I2_PROJECTS"

# Alternative old path (some projects may have this)
ALT_OLD_PATH_PREFIX="/mnt/projects"

# Check what would be changed
check_paths() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Checking Project Directory Paths (Dry Run)${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    echo -e "${BLUE}Old path prefix: $OLD_PATH_PREFIX${NC}"
    echo -e "${BLUE}New path prefix: $NEW_PATH_PREFIX${NC}"
    echo ""

    # Run Django management command to check
    python /usr/src/app/manage.py shell << 'PYTHON_EOF'
import sys
from ccp4i2.db.models import Project

OLD_PREFIX = "/usr/src/app/mydata/CompoundDatabaseData/CCP4I2_PROJECTS"
NEW_PREFIX = "/mnt/projects/CCP4I2_PROJECTS"

# Find projects with old path
projects_to_update = Project.objects.filter(directory__startswith=OLD_PREFIX)
count = projects_to_update.count()

if count == 0:
    print(f"\033[0;32mNo projects found with old path prefix.\033[0m")
    print(f"\nChecking current path patterns...")

    # Show sample of current paths
    all_projects = Project.objects.all()[:10]
    if all_projects:
        print(f"\nSample project directories (first 10):")
        for p in all_projects:
            print(f"  {p.name}: {p.directory}")
else:
    print(f"\033[1;33mFound {count} projects to update:\033[0m\n")

    for project in projects_to_update[:20]:  # Show first 20
        old_dir = project.directory
        new_dir = old_dir.replace(OLD_PREFIX, NEW_PREFIX)
        print(f"  {project.name}:")
        print(f"    \033[0;31mOLD: {old_dir}\033[0m")
        print(f"    \033[0;32mNEW: {new_dir}\033[0m")
        print()

    if count > 20:
        print(f"  ... and {count - 20} more projects")

    print(f"\n\033[1;33mRun with 'patch' to apply these changes.\033[0m")

PYTHON_EOF

    echo ""
}

# Perform the patch
do_patch() {
    echo -e "${CYAN}========================================${NC}"
    echo -e "${CYAN}  Patching Project Directory Paths${NC}"
    echo -e "${CYAN}========================================${NC}"
    echo ""

    echo -e "${BLUE}Old path prefix: $OLD_PATH_PREFIX${NC}"
    echo -e "${BLUE}New path prefix: $NEW_PATH_PREFIX${NC}"
    echo ""

    # Run Django management command to patch
    python /usr/src/app/manage.py shell << 'PYTHON_EOF'
import sys
from django.db import transaction
from ccp4i2.db.models import Project

OLD_PREFIX = "/usr/src/app/mydata/CompoundDatabaseData/CCP4I2_PROJECTS"
NEW_PREFIX = "/mnt/projects/CCP4I2_PROJECTS"

# Find projects with old path
projects_to_update = Project.objects.filter(directory__startswith=OLD_PREFIX)
count = projects_to_update.count()

if count == 0:
    print(f"\033[0;32mNo projects found with old path prefix. Nothing to update.\033[0m")
else:
    print(f"\033[1;33mUpdating {count} projects...\033[0m\n")

    updated = 0
    errors = 0

    with transaction.atomic():
        for project in projects_to_update:
            old_dir = project.directory
            new_dir = old_dir.replace(OLD_PREFIX, NEW_PREFIX)

            try:
                project.directory = new_dir
                project.save(update_fields=['directory'])
                print(f"\033[0;32m  ✓ {project.name}\033[0m")
                updated += 1
            except Exception as e:
                print(f"\033[0;31m  ✗ {project.name}: {e}\033[0m")
                errors += 1

    print(f"\n\033[0;32mUpdated: {updated}\033[0m")
    if errors > 0:
        print(f"\033[0;31mErrors: {errors}\033[0m")

    print(f"\n\033[0;32mDatabase patch complete!\033[0m")

PYTHON_EOF

    echo ""
    echo -e "${YELLOW}Next steps:${NC}"
    echo "1. Verify the changes: $0 check"
    echo "2. Restart the server container if needed"
}

# Show usage
show_usage() {
    echo "CCP4i2 Project Directory Path Patch Script"
    echo ""
    echo "Updates Project.directory fields in the Django database."
    echo ""
    echo "Path transformation:"
    echo "  OLD: $OLD_PATH_PREFIX/{project_name}"
    echo "  NEW: $NEW_PATH_PREFIX/{project_name}"
    echo ""
    echo "This script should be run INSIDE a server container instance."
    echo ""
    echo "Usage: $0 <command>"
    echo ""
    echo "Commands:"
    echo "  check   Show what would be changed (dry run)"
    echo "  patch   Perform the actual database update"
    echo ""
}

# Main command handler
case "${1:-}" in
    check)
        check_paths
        ;;
    patch)
        do_patch
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
