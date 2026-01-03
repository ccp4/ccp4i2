#!/usr/bin/env python
"""
Standalone test that creates database entries and shows them.
Does NOT use pytest, so database persists for inspection.
"""
import os
import asyncio
import tempfile
import uuid as uuid_module
from pathlib import Path

from ccp4i2.core import CCP4Utils

# Setup temp database and projects
temp_db = tempfile.NamedTemporaryFile(suffix='.sqlite3', delete=False)
temp_db.close()
temp_projects = tempfile.mkdtemp(prefix='standalone_test_')

os.environ['DJANGO_SETTINGS_MODULE'] = 'ccp4i2.config.test_settings'
os.environ['CCP4I2_DB_FILE'] = temp_db.name
os.environ['CCP4I2_PROJECTS_DIR'] = temp_projects

# Import Django and setup
import django
django.setup()

from django.core.management import call_command
from ccp4i2.db import models
from ccp4i2.db.async_db_handler import AsyncDatabaseHandler
from asgiref.sync import sync_to_async

print(f"🗄️  Database: {temp_db.name}")
print(f"📁 Projects: {temp_projects}")
print()

# Run migrations
print("Running migrations...")
call_command('migrate', '--run-syncdb', verbosity=0)

async def main():
    # Create project
    project = await sync_to_async(models.Project.objects.create)(
        name="standalone_test",
        directory=Path(temp_projects) / "standalone_test"
    )
    print(f"✓ Created project: {project.name} ({project.uuid})\n")

    # Load task library and create plugin
    from ccp4i2.core.CCP4TaskManager import TASKMANAGER
    task_mgr = TASKMANAGER()
    task_mgr.loadTaskLibrary(str(Path(__file__).parent.parent / "wrappers"))

    # Get the plugin class for ctruncate
    ctruncate_class = task_mgr.get_plugin_class('ctruncate')
    plugin = ctruncate_class(name='test_ctruncate', parent=None)

    # Configure input
    demo_data = Path(__file__).parent.parent / "demo_data/gamma/merged_intensities_native.mtz"
    plugin.container.inputData.HKLIN.setFullPath(str(demo_data))
    plugin.container.inputData.ISIGIanom.columnNames = ["Iplus", "SIGIplus", "Iminus", "SIGIminus"]

    # Setup database tracking and execute
    handler = AsyncDatabaseHandler(project_uuid=project.uuid)
    job_uuid = await handler.create_job(
        task_name="ctruncate",
        title="Ctruncate with Database Tracking",
        job_number="1"
    )

    print(f"⚙️  Executing ctruncate (Job {job_uuid})...")
    await handler.track_job(job_uuid, plugin)
    print("✓ Execution complete!\n")

    # Now inspect the database - it should have data
    os.system(f"python {Path(__file__).parent}/inspect_test_db.py {temp_db.name}")

    print(f"\n🧹 Cleanup commands:")
    print(f"   rm -rf {temp_projects}")
    print(f"   rm {temp_db.name}")

if __name__ == "__main__":
    asyncio.run(main())
