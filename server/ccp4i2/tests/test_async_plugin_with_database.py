"""
End-to-end integration test: Async plugin execution with database tracking.

This test combines:
- Real CCP4 plugin execution (ctruncate)
- Django database tracking
- Modern async database handler
- Real crystallographic test data

Run with:
    ./run_db_tests.sh
"""

import pytest
import os
import sys
from pathlib import Path

# Check for CCP4I2_ROOT
CCP4I2_ROOT = os.environ.get('CCP4I2_ROOT')
if not CCP4I2_ROOT:
    pytest.skip("CCP4I2_ROOT not set", allow_module_level=True)

# Check for demo data
demo_data_path = Path(CCP4I2_ROOT) / 'demo_data' / 'gamma'
test_mtz = demo_data_path / 'merged_intensities_native.mtz'
if not test_mtz.exists():
    pytest.skip(f"Test data not found: {test_mtz}", allow_module_level=True)

# Add server directory to path for Django imports
server_dir = Path(__file__).parent.parent / 'server'
if str(server_dir) not in sys.path:
    sys.path.insert(0, str(server_dir))

# Set Django settings
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'ccp4i2.config.test_settings')


@pytest.fixture
def projects_base():
    """Get the base directory for test projects"""
    projects_dir = Path(os.environ.get('CCP4I2_PROJECTS_DIR', Path.home() / '.ccp4i2_test' / 'test_projects'))
    projects_dir.mkdir(parents=True, exist_ok=True)
    return projects_dir


@pytest.fixture
def test_project(projects_base, db):
    """Create a test project for testing"""
    from ccp4i2.db import models
    import shutil
    import uuid

    # Create unique project directory
    project_name = f"test_project_{uuid.uuid4().hex[:8]}"
    project_dir = projects_base / project_name
    project_dir.mkdir(parents=True, exist_ok=True)

    # Create project in database
    project = models.Project.objects.create(
        name=project_name,
        directory=str(project_dir),
        description="Test project for async plugin execution",
    )

    print(f"\n{'='*70}")
    print(f"Created test project: {project.name}")
    print(f"Project UUID: {project.uuid}")
    print(f"Project directory: {project_dir}")
    print(f"{'='*70}\n")

    yield project

    # Cleanup
    if project_dir.exists():
        shutil.rmtree(project_dir, ignore_errors=True)


@pytest.mark.django_db(transaction=True)
class TestAsyncPluginWithDatabase:
    """Test real plugin execution with database tracking."""

    @pytest.mark.asyncio
    async def test_ctruncate_with_database_tracking(self, test_project):
        """
        Run ctruncate plugin with full database tracking.

        This test:
        1. Creates a job in the database
        2. Runs ctruncate plugin (converts intensities to amplitudes)
        3. Tracks status changes automatically via signals
        4. Gleans output files and registers them in database
        5. Extracts and stores KPIs (performance indicators)

        Expected behavior:
        - Job status: PENDING → RUNNING → FINISHED
        - Output files: HKLOUT (MTZ file with F/SIGF)
        - KPIs: Various statistics from ctruncate
        """
        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler
        from ccp4i2.db import models
        from ccp4i2.core.CCP4PluginScript import CPluginScript
        from asgiref.sync import sync_to_async

        print(f"\n{'='*70}")
        print("RUNNING: Ctruncate with Async Database Tracking")
        print(f"{'='*70}")
        print(f"Input data: {test_mtz}")
        print(f"Project: {test_project.name} ({test_project.uuid})")
        print(f"{'='*70}\n")

        # Create handler
        handler = AsyncDatabaseHandler(project_uuid=test_project.uuid)

        # Create plugin using TASKMANAGER (loads task definition properly)
        from ccp4i2.core.CCP4TaskManager import TASKMANAGER

        # Get the plugin class for ctruncate
        ctruncate_class = TASKMANAGER().get_plugin_class('ctruncate')
        assert ctruncate_class is not None, "ctruncate plugin not found in registry"

        # Create plugin instance (this will load .def.xml via TASKNAME)
        plugin = ctruncate_class(name="test_ctruncate")

        # Set input parameters
        plugin.container.inputData.HKLIN.setFullPath(str(test_mtz))

        # Configure for anomalous intensity data (Iplus/Iminus)
        # The test data has Iplus, SIGIplus, Iminus, SIGIminus columns
        # Use .set() to properly configure the column mapping
        plugin.container.inputData.ISIGIanom.set({
            'Ip': 'Iplus',
            'SIGIp': 'SIGIplus',
            'Im': 'Iminus',
            'SIGIm': 'SIGIminus'
        })

        # Configure output to generate mini-MTZ with FMEAN format
        plugin.container.controlParameters.OUTPUTMINIMTZ.set(True)
        plugin.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG.set(4)  # FMEAN

        print(f"Plugin created: {plugin.container.objectName()}")
        print(f"Input file: {plugin.container.inputData.HKLIN}")
        print(f"ISIGIanom isSet: {plugin.container.inputData.ISIGIanom.isSet()}")
        print(f"ISIGIanom columnNames: {plugin.container.inputData.ISIGIanom.columnNames}")
        # Also check if ISIGI is set (might be needed by wrapper)
        print(f"ISIGI isSet: {plugin.container.inputData.ISIGI.isSet() if hasattr(plugin.container.inputData, 'ISIGI') else 'N/A'}")
        print()

        # Track job execution
        print("Starting job execution with automatic tracking...")

        # Debug: Check work directory before tracking
        print(f"DEBUG: plugin.workDirectory BEFORE tracking = {plugin.workDirectory}")

        async with handler.track_job(plugin):
            # Debug: Check work directory after tracking
            print(f"DEBUG: plugin.workDirectory AFTER tracking = {plugin.workDirectory}")

            # Execute the plugin (process() is synchronous, wrap with sync_to_async)
            result = await sync_to_async(plugin.process)()
            print(f"\n✓ Plugin execution completed")
            print(f"  Result: {result}")
            print(f"  Status: {plugin.get_status()}")

            # Debug: List files in work directory
            import os
            if os.path.exists(plugin.workDirectory):
                files_in_workdir = os.listdir(plugin.workDirectory)
                print(f"DEBUG: Files in work directory: {files_in_workdir[:10] if len(files_in_workdir) > 10 else files_in_workdir}")

        print(f"\n{'='*70}")
        print("VERIFICATION: Checking Database Records")
        print(f"{'='*70}\n")

        # Verify job was created
        jobs = await sync_to_async(list)(
            models.Job.objects.filter(project=test_project)
        )
        assert len(jobs) == 1, f"Expected 1 job, found {len(jobs)}"
        job = jobs[0]

        print(f"✓ Job created in database:")
        print(f"  UUID: {job.uuid}")
        print(f"  Number: {job.number}")
        print(f"  Task: {job.task_name}")
        print(f"  Title: {job.title}")
        print(f"  Status: {job.get_status_display()}")
        # Access directory property in sync context
        job_directory = await sync_to_async(lambda: job.directory)()
        print(f"  Directory: {job_directory}")

        # Verify job status
        assert job.task_name == "ctruncate"
        assert job.number == "1"

        # Debug: Show actual status values
        print(f"\nDEBUG: plugin.get_status() = {plugin.get_status()}")
        print(f"DEBUG: CPluginScript.SUCCEEDED = {CPluginScript.SUCCEEDED}")
        print(f"DEBUG: job.status = {job.status} (expected: {models.Job.Status.FINISHED})")
        print(f"DEBUG: models.Job.Status.RUNNING = {models.Job.Status.RUNNING}")

        # Refresh job from database to see if status was updated
        job = await sync_to_async(models.Job.objects.get)(uuid=job.uuid)
        print(f"DEBUG: After refresh, job.status = {job.status}")

        assert job.status == models.Job.Status.FINISHED, f"Expected FINISHED, got {job.status}"
        assert job.finish_time is not None
        print(f"  Finish time: {job.finish_time}")
        print()

        # Verify output files were gleaned
        # Note: role is on FileUse, not File, so we need to filter through the relationship
        output_file_uses = await sync_to_async(list)(
            models.FileUse.objects.filter(job=job, role=models.FileUse.Role.OUT).select_related('file')
        )
        # Access file relationship in sync context
        output_files = await sync_to_async(lambda: [file_use.file for file_use in output_file_uses])()
        print(f"✓ Output files registered: {len(output_files)}")

        # Access file properties in sync context
        @sync_to_async
        def print_file_info(file):
            print(f"  - {file.job_param_name}: {file.type.name if file.type else 'unknown'}")
            print(f"    Name: {file.name}")
            # Note: full_path is a property that might need path computation
            if hasattr(file, 'path') and file.path and Path(file.path).exists():
                size = Path(file.path).stat().st_size
                print(f"    Size: {size:,} bytes")
            print()

        for file in output_files:
            await print_file_info(file)

        # Verify at least HKLOUT exists
        hklout_files = await sync_to_async(lambda: [f for f in output_files if f.job_param_name == 'HKLOUT'])()
        assert len(hklout_files) >= 1, "Expected HKLOUT output file"

        # Verify KPIs were registered
        float_kpis = await sync_to_async(list)(
            models.JobFloatValue.objects.filter(job=job)
        )
        char_kpis = await sync_to_async(list)(
            models.JobCharValue.objects.filter(job=job)
        )

        print(f"✓ KPIs registered:")
        print(f"  Float values: {len(float_kpis)}")
        for kpi in float_kpis[:5]:  # Show first 5
            print(f"    {kpi.key.key_name}: {kpi.value}")
        if len(float_kpis) > 5:
            print(f"    ... and {len(float_kpis) - 5} more")

        print(f"  String values: {len(char_kpis)}")
        for kpi in char_kpis[:5]:  # Show first 5
            print(f"    {kpi.key.key_name}: {kpi.value}")
        if len(char_kpis) > 5:
            print(f"    ... and {len(char_kpis) - 5} more")
        print()

        # ========== DETAILED DATABASE INSPECTION ==========
        print("\n" + "="*120)
        print("DATABASE CONTENTS (DETAILED)")
        print("="*120)

        # Show all tables with data
        @sync_to_async
        def show_db_tables():
            print(f"\n📊 Projects: {models.Project.objects.count()}")
            print(f"📊 Jobs: {models.Job.objects.count()}")
            print(f"📊 Files: {models.File.objects.count()}")
            print(f"📊 File Uses: {models.FileUse.objects.count()}")
            print(f"📊 File Types: {models.FileType.objects.count()}")
            print(f"📊 Float KPIs: {models.JobFloatValue.objects.count()}")
            print(f"📊 Char KPIs: {models.JobCharValue.objects.count()}")

            # Show project details in table format
            project = models.Project.objects.order_by('-creation_time').first()
            if project:
                print(f"\n📁 PROJECT DETAILS")
                print(f"{'='*120}")
                print(f"  UUID       : {project.uuid}")
                print(f"  Name       : {project.name}")
                print(f"  Directory  : {project.directory}")

                # Show job details
                print(f"\n⚙️  JOBS")
                print(f"{'='*120}")
                print(f"{'#':<6} {'Task':<15} {'Status':<12} {'Created':<25} {'Finished':<25}")
                print("-" * 120)
                for job in models.Job.objects.filter(project=project):
                    print(f"{job.number:<6} {job.task_name:<15} {job.get_status_display():<12} "
                          f"{str(job.creation_time) if job.creation_time else 'N/A':<25} "
                          f"{str(job.finish_time) if job.finish_time else 'N/A':<25}")

                # Show file details
                print(f"\n📄 FILES")
                print(f"{'='*120}")
                print(f"{'Job':<6} {'Param':<15} {'Filename':<25} {'Type':<35} {'Directory':<15}")
                print("-" * 120)
                for file in models.File.objects.filter(job__project=project).select_related('type', 'job'):
                    print(f"{file.job.number:<6} {(file.job_param_name or 'N/A')[:14]:<15} "
                          f"{file.name[:24]:<25} {(file.type.name if file.type else 'unknown')[:34]:<35} "
                          f"{file.get_directory_display():<15}")

                # Show file uses
                print(f"\n🔗 FILE USES")
                print(f"{'='*120}")
                print(f"{'Job':<6} {'Role':<12} {'Param':<20} {'Filename':<30}")
                print("-" * 120)
                # Convert to list to eagerly evaluate the queryset before accessing relationships
                file_uses = list(models.FileUse.objects.filter(job__project=project).select_related('file', 'job'))
                for use in file_uses:
                    try:
                        filename = use.file.name if use.file else 'N/A'
                    except Exception:
                        filename = 'N/A'
                    print(f"{use.job.number:<6} {use.get_role_display():<12} "
                          f"{(use.job_param_name or 'N/A')[:19]:<20} {filename[:29]:<30}")

                # Show file hierarchy
                proj_dir = Path(project.directory)
                if proj_dir.exists():
                    print(f"\n📂 FILE SYSTEM HIERARCHY")
                    print(f"{'='*120}")
                    print(f"\n{proj_dir}\n")

                    def show_tree(path, prefix="", is_last=True):
                        items = sorted(path.iterdir(), key=lambda p: (not p.is_dir(), p.name))
                        for i, item in enumerate(items):
                            is_last_item = (i == len(items) - 1)
                            connector = "└── " if is_last_item else "├── "
                            new_prefix = prefix + ("    " if is_last_item else "│   ")
                            if item.is_file():
                                size = item.stat().st_size
                                size_str = f"{size/(1024*1024):.2f} MB" if size > 1024*1024 else f"{size/1024:.2f} KB" if size > 1024 else f"{size} bytes"
                                print(f"{prefix}{connector}{item.name} ({size_str})")
                            else:
                                print(f"{prefix}{connector}{item.name}/")
                                if item.is_dir():
                                    show_tree(item, new_prefix, is_last_item)

                    show_tree(proj_dir)

        await show_db_tables()
        print("\n" + "="*120 + "\n")

        print(f"{'='*70}")
        print("✓ TEST PASSED: Full async plugin execution with database tracking")
        print(f"{'='*70}\n")

    @pytest.mark.asyncio
    async def test_ctruncate_manual_control(self, test_project):
        """
        Test manual control of job lifecycle (without context manager).

        This demonstrates lower-level control where you manually:
        - Create the job
        - Update status
        - Import input files
        - Execute plugin
        - Glean output files
        - Register KPIs
        """
        from ccp4i2.db.async_db_handler import AsyncDatabaseHandler
        from ccp4i2.db import models
        from ccp4i2.core.CCP4PluginScript import CPluginScript
        from asgiref.sync import sync_to_async

        print(f"\n{'='*70}")
        print("RUNNING: Manual Job Control Test")
        print(f"{'='*70}\n")

        handler = AsyncDatabaseHandler(project_uuid=test_project.uuid)

        # Step 1: Create job manually
        print("Step 1: Creating job...")
        job = await handler.create_job(
            task_name="ctruncate",
            title="Manual control test"
        )
        print(f"  ✓ Created job {job.number}: {job.title}")
        print(f"    UUID: {job.uuid}")
        print()

        # Step 2: Create and configure plugin
        print("Step 2: Configuring plugin...")
        from ccp4i2.core.CCP4TaskManager import TASKMANAGER

        # Get the plugin class
        ctruncate_class = TASKMANAGER().get_plugin_class('ctruncate')
        assert ctruncate_class is not None, "ctruncate plugin not found"

        # Create plugin instance
        plugin = ctruncate_class(name="manual_ctruncate")
        plugin.container.inputData.HKLIN.setFullPath(str(test_mtz))

        # Configure for anomalous intensity data - use .set() method
        plugin.container.inputData.ISIGIanom.set({
            'Ip': 'Iplus',
            'SIGIp': 'SIGIplus',
            'Im': 'Iminus',
            'SIGIm': 'SIGIminus'
        })
        plugin.container.controlParameters.OUTPUTMINIMTZ.set(True)
        plugin.container.controlParameters.OUTPUTMINIMTZCONTENTFLAG.set(4)

        print(f"  ✓ Plugin configured")
        print()

        # Step 3: Update status to RUNNING
        print("Step 3: Starting execution...")
        await handler.update_job_status(job.uuid, models.Job.Status.RUNNING)
        print(f"  ✓ Job status: RUNNING")
        print()

        # Step 4: Execute plugin
        print("Step 4: Running plugin...")
        result = await sync_to_async(plugin.process)()
        print(f"  ✓ Plugin completed: {result}")
        print(f"    Status: {plugin.get_status()}")
        print()

        # Step 5: Glean files and KPIs
        print("Step 5: Gleaning files and KPIs...")

        # Debug: Show what's in outputData container
        print(f"DEBUG: outputData type = {type(plugin.container.outputData)}")
        print(f"DEBUG: outputData children = {plugin.container.outputData.children() if hasattr(plugin.container.outputData, 'children') else 'N/A'}")
        if hasattr(plugin.container.outputData, 'children'):
            for child in plugin.container.outputData.children():
                print(f"  - Child: {child.name if hasattr(child, 'name') else child}, type={type(child).__name__}")
                if hasattr(child, 'isSet'):
                    print(f"    isSet: {child.isSet()}")
                if hasattr(child, 'fullPath'):
                    # fullPath is a property, not a method
                    try:
                        print(f"    fullPath: {child.fullPath}")
                    except Exception as e:
                        print(f"    fullPath error: {e}")
                if hasattr(child, 'exists'):
                    try:
                        print(f"    exists(): {child.exists()}")
                    except Exception as e:
                        print(f"    exists() error: {e}")
                else:
                    print(f"    No exists() method")

        files = await handler.glean_job_files(job.uuid, plugin.container.outputData)
        print(f"  Gleaned {len(files)} files")
        await handler.glean_performance_indicators(job.uuid, plugin.container.outputData)
        print(f"  ✓ Files and KPIs gleaned")
        print()

        # Step 6: Update status to FINISHED
        print("Step 6: Marking job as finished...")
        await handler.update_job_status(job.uuid, models.Job.Status.FINISHED)
        print(f"  ✓ Job status: FINISHED")
        print()

        # Verify
        updated_job = await sync_to_async(models.Job.objects.get)(uuid=job.uuid)
        assert updated_job.status == models.Job.Status.FINISHED

        output_file_uses = await sync_to_async(list)(
            models.FileUse.objects.filter(job=job, role=models.FileUse.Role.OUT).select_related('file')
        )
        output_files = [file_use.file for file_use in output_file_uses]
        assert len(output_files) > 0, "Expected output files"

        print(f"{'='*70}")
        print("✓ TEST PASSED: Manual job control workflow")
        print(f"{'='*70}\n")
