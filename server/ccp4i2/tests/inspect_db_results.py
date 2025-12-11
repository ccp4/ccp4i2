#!/usr/bin/env python
"""
Inspect database contents after running ctruncate with database tracking.
"""
import os
import sys
from pathlib import Path

# Setup Django
sys.path.insert(0, str(Path(__file__).parent.parent / "server"))
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4x.config.test_settings')

import django
django.setup()

from ccp4x.db import models
from asgiref.sync import sync_to_async
import asyncio
import tempfile
import uuid as uuid_module


async def inspect_database(project_uuid: uuid_module.UUID):
    """Inspect and display all database entries for a project."""

    print("\n" + "="*80)
    print("DATABASE INSPECTION REPORT")
    print("="*80)

    # 1. PROJECT
    print("\n📁 PROJECT")
    print("-" * 80)
    project = await sync_to_async(models.Project.objects.get)(uuid=project_uuid)
    project_data = await sync_to_async(lambda: {
        "UUID": str(project.uuid),
        "Name": project.name,
        "Title": project.title or "N/A",
        "Directory": project.directory,
        "Created": project.created_at.strftime("%Y-%m-%d %H:%M:%S"),
    })()
    for key, value in project_data.items():
        print(f"  {key:15} : {value}")

    # 2. JOBS
    print("\n⚙️  JOBS")
    print("-" * 80)
    jobs = await sync_to_async(list)(
        models.Job.objects.filter(project=project).order_by('created_at')
    )

    print(f"{'Job #':<8} {'Task':<15} {'Title':<25} {'Status':<12} {'Started':<10} {'Finished':<10}")
    print("-" * 95)
    for job in jobs:
        job_data = await sync_to_async(lambda j=job: (
            str(j.number),
            j.task_name,
            (j.title or "N/A")[:24],
            j.get_status_display(),
            j.created_at.strftime("%H:%M:%S"),
            j.finish_time.strftime("%H:%M:%S") if j.finish_time else "N/A",
        ))()
        print(f"{job_data[0]:<8} {job_data[1]:<15} {job_data[2]:<25} {job_data[3]:<12} {job_data[4]:<10} {job_data[5]:<10}")

    # 3. FILES
    print("\n📄 FILES")
    print("-" * 120)
    files = await sync_to_async(list)(
        models.File.objects.filter(job__project=project).select_related('type', 'job')
    )

    print(f"{'Job':<6} {'Parameter':<15} {'Filename':<20} {'Type':<30} {'Directory':<15} {'Size':<15}")
    print("-" * 120)
    for file in files:
        file_data = await sync_to_async(lambda f=file: (
            f.job.number,
            (f.job_param_name or "N/A")[:14],
            f.name[:19],
            (f.type.name if f.type else "unknown")[:29],
            f.get_directory_display(),
            f"{Path(f.path).stat().st_size:,}" if f.path and Path(f.path).exists() else "N/A",
        ))()
        print(f"{file_data[0]:<6} {file_data[1]:<15} {file_data[2]:<20} {file_data[3]:<30} {file_data[4]:<15} {file_data[5]:<15}")

    # 4. FILE USES (Input/Output relationships)
    print("\n🔗 FILE USES (Input/Output)")
    print("-" * 80)
    file_uses = await sync_to_async(list)(
        models.FileUse.objects.filter(job__project=project).select_related('file', 'job')
    )

    print(f"{'Job':<6} {'Role':<10} {'Parameter':<20} {'Filename':<30}")
    print("-" * 80)
    for use in file_uses:
        use_data = await sync_to_async(lambda u=use: (
            u.job.number,
            u.get_role_display(),
            (u.job_param_name or "N/A")[:19],
            u.file.name[:29],
        ))()
        print(f"{use_data[0]:<6} {use_data[1]:<10} {use_data[2]:<20} {use_data[3]:<30}")

    # 5. KPIs - Float Values
    print("\n📊 KEY PERFORMANCE INDICATORS (Float Values)")
    print("-" * 100)
    float_kpis = await sync_to_async(list)(
        models.JobFloatValue.objects.filter(job__project=project).select_related('key', 'job')
    )

    if float_kpis:
        print(f"{'Job':<6} {'Key':<30} {'Value':<15} {'Description':<40}")
        print("-" * 100)
        for kpi in float_kpis:
            kpi_data = await sync_to_async(lambda k=kpi: (
                k.job.number,
                k.key.name[:29],
                f"{k.value:.4f}",
                (k.key.description or "N/A")[:39],
            ))()
            print(f"{kpi_data[0]:<6} {kpi_data[1]:<30} {kpi_data[2]:<15} {kpi_data[3]:<40}")
    else:
        print("No float KPIs registered")

    # 6. KPIs - Char Values
    print("\n📊 KEY PERFORMANCE INDICATORS (Text Values)")
    print("-" * 100)
    char_kpis = await sync_to_async(list)(
        models.JobCharValue.objects.filter(job__project=project).select_related('key', 'job')
    )

    if char_kpis:
        print(f"{'Job':<6} {'Key':<30} {'Value':<40} {'Description':<20}")
        print("-" * 100)
        for kpi in char_kpis:
            kpi_data = await sync_to_async(lambda k=kpi: (
                k.job.number,
                k.key.name[:29],
                (k.value[:37] + "...") if len(k.value) > 40 else k.value,
                (k.key.description or "N/A")[:19],
            ))()
            print(f"{kpi_data[0]:<6} {kpi_data[1]:<30} {kpi_data[2]:<40} {kpi_data[3]:<20}")
    else:
        print("No text KPIs registered")

    # 7. FILE TYPES
    print("\n🏷️  FILE TYPES")
    print("-" * 80)
    file_types = await sync_to_async(list)(
        models.FileType.objects.filter(
            file__job__project=project
        ).distinct()
    )

    print(f"{'MIME Type':<40} {'Description':<35}")
    print("-" * 80)
    for ft in file_types:
        type_data = await sync_to_async(lambda t=ft: (
            t.name[:39],
            (t.description or "N/A")[:34],
        ))()
        print(f"{type_data[0]:<40} {type_data[1]:<35}")

    print("\n" + "="*80)
    print("END OF DATABASE INSPECTION")
    print("="*80 + "\n")


def show_file_hierarchy(root_path: Path):
    """Show file hierarchy as a tree."""

    print("\n" + "="*80)
    print("FILE SYSTEM HIERARCHY")
    print("="*80)
    print(f"\n📂 {root_path}")

    def build_tree(path: Path, prefix: str = "", is_last: bool = True):
        """Recursively build tree representation."""
        if not path.exists():
            return

        items = sorted(path.iterdir(), key=lambda p: (not p.is_dir(), p.name))

        for i, item in enumerate(items):
            is_last_item = (i == len(items) - 1)

            # Determine connector
            if is_last_item:
                connector = "└── "
                new_prefix = prefix + "    "
            else:
                connector = "├── "
                new_prefix = prefix + "│   "

            # Get size info
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

                # Recurse for directories
                if item.is_dir():
                    build_tree(item, new_prefix, is_last_item)

    build_tree(root_path)

    print("\n" + "="*80 + "\n")


async def main():
    """Run the test and inspect results."""
    from ccp4x.db.async_db_handler import AsyncDatabaseHandler
    from ccp4i2.core.CCP4TaskManager import TASKMANAGER

    # Setup temporary database and project directory
    temp_db = tempfile.NamedTemporaryFile(suffix='.sqlite3', delete=False)
    temp_db.close()
    os.environ['CCP4I2_DB_FILE'] = temp_db.name

    temp_projects = tempfile.mkdtemp(prefix='test_projects_')
    os.environ['CCP4I2_PROJECTS_DIR'] = temp_projects

    print(f"\n🗄️  Test database: {temp_db.name}")
    print(f"📁 Test projects: {temp_projects}\n")

    # Re-setup Django with new DB
    from django.conf import settings
    from django.core.management import call_command

    # Run migrations
    print("Running migrations...")
    call_command('migrate', '--run-syncdb', verbosity=0)

    # Create test project
    project_name = f"test_project_{uuid_module.uuid4().hex[:8]}"
    project = await sync_to_async(models.Project.objects.create)(
        name=project_name,
        title="Database Inspection Test",
        directory=Path(temp_projects) / project_name
    )

    # Load task library
    TASKMANAGER.loadTaskLibrary(
        str(Path(__file__).parent.parent / "wrappers")
    )

    # Create and run plugin
    from ccp4i2.wrappers.ctruncate.script.ctruncate import ctruncate

    plugin = ctruncate(name='test_ctruncate', parent=None)

    # Configure input
    demo_data = Path(__file__).parent.parent / "demo_data/gamma/merged_intensities_native.mtz"
    plugin.container.inputData.HKLIN.setFullPath(str(demo_data))
    plugin.container.inputData.ISIGIanom.columnNames = ["Iplus", "SIGIplus", "Iminus", "SIGIminus"]

    # Setup database tracking
    handler = AsyncDatabaseHandler(project_uuid=project.uuid)
    job_uuid = await handler.create_job(
        task_name="ctruncate",
        title="Ctruncate with Database Tracking",
        job_number="1"
    )

    # Execute with tracking
    print(f"\n⚙️  Executing ctruncate (Job {job_uuid})...\n")
    await handler.track_job(job_uuid, plugin)

    print("\n✅ Execution complete!")

    # Inspect database
    await inspect_database(project.uuid)

    # Show file hierarchy
    show_file_hierarchy(project.directory)

    # Cleanup info
    print(f"\n🧹 Cleanup commands:")
    print(f"   rm -rf {temp_projects}")
    print(f"   rm {temp_db.name}\n")


if __name__ == "__main__":
    asyncio.run(main())
