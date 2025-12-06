#!/bin/bash
# Show database and file contents from last test run

set -e

# Determine project root from script location
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Find and source CCP4 environment
for dir in "$PROJECT_ROOT"/../ccp4-*/bin/ccp4.setup-sh; do
    if [ -f "$dir" ]; then
        source "$dir"
        break
    fi
done

export CCP4I2_ROOT="$PROJECT_ROOT"
export PYTHONPATH="$PROJECT_ROOT/server:$PROJECT_ROOT:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE=ccp4x.config.test_settings

# Create temp database and projects
TMP_DB=$(mktemp -t test_db_XXXXXX).sqlite3
TMP_PROJECTS=$(mktemp -d -t test_projects_XXXXXX)
export CCP4I2_DB_FILE=$TMP_DB
export CCP4I2_PROJECTS_DIR=$TMP_PROJECTS

echo "🗄️  Test database: $TMP_DB"
echo "📁 Test projects: $TMP_PROJECTS"
echo ""

# Run migrations
cd server
python manage.py migrate --run-syncdb 2>&1 > /dev/null
cd ..

echo "⚙️  Running ctruncate test with database tracking..."
echo ""

# Run the test (suppress debug output)
python -m pytest tests/test_async_plugin_with_database.py::TestAsyncPluginWithDatabase::test_ctruncate_with_database_tracking \
    -v -s 2>&1 | grep -E "(PASSED|FAILED|===|Running:|✓|Project UUID:|Project directory:)" | head -20

echo ""
echo "================================================================================
="
echo "DATABASE INSPECTION"
echo "================================================================================"

# Now inspect the database using Python
python << 'PYEOF'
import os
import sys
import django
from pathlib import Path

# Django is already set up from env vars
django.setup()

from ccp4x.db import models

# Get the most recent project
project = models.Project.objects.order_by('-creation_time').first()

if not project:
    print("No projects found in database")
    sys.exit(1)

print("\n📁 PROJECT")
print("-" * 80)
print(f"  UUID            : {project.uuid}")
print(f"  Name            : {project.name}")
print(f"  Directory       : {project.directory}")
print(f"  Created         : {project.creation_time.strftime('%Y-%m-%d %H:%M:%S')}")

print("\n⚙️  JOBS")
print("-" * 95)
print(f"{'Job #':<8} {'Task':<15} {'Title':<25} {'Status':<12} {'Started':<10} {'Finished':<10}")
print("-" * 95)
for job in models.Job.objects.filter(project=project):
    print(f"{job.number:<8} {job.task_name:<15} {(job.title or 'N/A')[:24]:<25} {job.get_status_display():<12} "
          f"{job.start_time.strftime('%H:%M:%S') if job.start_time else 'N/A':<10} {job.finish_time.strftime('%H:%M:%S') if job.finish_time else 'N/A':<10}")

print("\n📄 FILES")
print("-" * 120)
print(f"{'Job':<6} {'Parameter':<15} {'Filename':<20} {'Type':<30} {'Directory':<15} {'Size':<15}")
print("-" * 120)
for file in models.File.objects.filter(job__project=project).select_related('type', 'job'):
    size = f"{Path(file.path).stat().st_size:,}" if file.path and Path(file.path).exists() else "N/A"
    print(f"{file.job.number:<6} {(file.job_param_name or 'N/A')[:14]:<15} {file.name[:19]:<20} "
          f"{(file.type.name if file.type else 'unknown')[:29]:<30} {file.get_directory_display():<15} {size:<15}")

print("\n🔗 FILE USES (Input/Output)")
print("-" * 80)
print(f"{'Job':<6} {'Role':<10} {'Parameter':<20} {'Filename':<30}")
print("-" * 80)
for use in models.FileUse.objects.filter(job__project=project).select_related('file', 'job'):
    print(f"{use.job.number:<6} {use.get_role_display():<10} {(use.job_param_name or 'N/A')[:19]:<20} {use.file.name[:29]:<30}")

print("\n🏷️  FILE TYPES")
print("-" * 80)
print(f"{'MIME Type':<40} {'Description':<35}")
print("-" * 80)
for ft in models.FileType.objects.filter(file__job__project=project).distinct():
    print(f"{ft.name[:39]:<40} {(ft.description or 'N/A')[:34]:<35}")

print("\n" + "="*80)
print("END OF DATABASE INSPECTION")
print("="*80 + "\n")

# Show file tree
project_dir = Path(project.directory)
if project_dir.exists():
    print("="*80)
    print("FILE SYSTEM HIERARCHY")
    print("="*80)
    print(f"\n📂 {project_dir}\n")

    def show_tree(path, prefix="", is_last=True):
        items = sorted(path.iterdir(), key=lambda p: (not p.is_dir(), p.name))
        for i, item in enumerate(items):
            is_last_item = (i == len(items) - 1)
            connector = "└── " if is_last_item else "├── "
            new_prefix = prefix + ("    " if is_last_item else "│   ")

            if item.is_file():
                size = item.stat().st_size
                if size > 1024 * 1024:
                    size_str = f"{size / (1024*1024):.2f} MB"
                elif size > 1024:
                    size_str = f"{size / 1024:.2f} KB"
                else:
                    size_str = f"{size} bytes"
                print(f"{prefix}{connector}{item.name} ({size_str})")
            else:
                print(f"{prefix}{connector}{item.name}/")
                if item.is_dir():
                    show_tree(item, new_prefix, is_last_item)

    show_tree(project_dir)
    print("\n" + "="*80 + "\n")

PYEOF

echo ""
echo "🧹 Cleanup:"
echo "   rm -rf $TMP_PROJECTS"
echo "   rm $TMP_DB"
