# Database Handler Usage Examples

This document compares the legacy `CCP4i2DjangoDbHandler` with the new `AsyncDatabaseHandler` to illustrate the cleaner, more modern API.

## Legacy Handler (Old API)

```python
from server.ccp4i2.db.ccp4i2_django_db_handler import CCP4i2DjangoDbHandler

# Create handler
handler = CCP4i2DjangoDbHandler()

# Create job (synchronous, returns job object but requires UUID conversion)
job = handler.createJob(
    pluginName="ctruncate",
    jobTitle="Convert intensities",
    parentJobId=str(parent_uuid),  # Must convert to string
    jobNumber=None
)

# Manual status tracking - plugin must call updateJobStatus explicitly
handler.updateJobStatus(
    jobId=str(job.uuid),  # Must convert to string
    status=None,
    finishStatus=CPluginScript.SUCCEEDED,
    container=plugin.outputData,
    dbOutputData=None
)

# File gleaning happens inside updateJobStatus (coupled)
# No way to register files independently
# No async support - blocks on database operations
```

## Modern Handler (New API)

### Basic Usage: Manual Control

```python
from server.ccp4i2.db.async_db_handler import AsyncDatabaseHandler
import uuid

# Create handler for a project
handler = AsyncDatabaseHandler(project_uuid=uuid.UUID("..."))

# Create job (async, type-safe)
job = await handler.create_job(
    task_name="ctruncate",
    title="Convert intensities",
    parent_job_uuid=parent_uuid,  # Direct UUID, no string conversion
)

# Update status (async, non-blocking)
await handler.update_job_status(
    job_uuid=job.uuid,
    status=models.Job.Status.RUNNING
)

# Register output files independently
file_obj = await handler.register_output_file(
    job_uuid=job.uuid,
    file_path=Path("output.mtz"),
    file_type="hklout",
    param_name="HKLOUT",
    content_flag=123,
)

# Register input files
await handler.register_input_file(
    job_uuid=job.uuid,
    file_uuid=input_file.uuid,
    param_name="HKLIN"
)
```

### Advanced Usage: Automatic Tracking

The most powerful feature is the `track_job()` context manager that handles everything automatically:

```python
from server.ccp4i2.db.async_db_handler import AsyncDatabaseHandler
from ccp4i2.core.CCP4PluginScript import CPluginScript

# Create plugin instance
plugin = CPluginScript(
    taskName="ctruncate",
    name="my_ctruncate_job",
    parent=parent_plugin
)

# Set up plugin parameters
plugin.inputData.HKLIN.set("/path/to/input.mtz")
plugin.inputData.COLANO.set("/*/*/[IMEAN,SIGIMEAN]")

# Create database handler
handler = AsyncDatabaseHandler(project_uuid=project.uuid)

# Track job lifecycle automatically!
async with handler.track_job(plugin):
    # Job record created in database
    # Status automatically updated as plugin runs
    # Signals connected for real-time tracking
    result = await plugin.execute()
    # Output files automatically registered on completion

# That's it! No manual status updates, no manual file registration
```

### Nested Job Tracking

Tracking nested jobs is seamless:

```python
# Parent job
parent_plugin = CPluginScript(taskName="copycell", name="parent_job")
parent_handler = AsyncDatabaseHandler(project_uuid=project.uuid)

async with parent_handler.track_job(parent_plugin):
    # Parent tracked automatically

    # Create child plugin
    child_plugin = parent_plugin.makePluginObject(
        taskName="ctruncate",
        name="child_job_1"
    )

    child_handler = AsyncDatabaseHandler(project_uuid=project.uuid)
    async with child_handler.track_job(child_plugin):
        # Child job automatically linked to parent
        # Job number auto-generated (e.g., "1.1")
        await child_plugin.execute()

    # Create another child
    child_plugin_2 = parent_plugin.makePluginObject(
        taskName="refmac",
        name="child_job_2"
    )

    async with child_handler.track_job(child_plugin_2):
        # Job number: "1.2"
        await child_plugin_2.execute()

    # Parent continues
    await parent_plugin.execute()
```

## Key Improvements

### 1. **Async/Await Throughout**
- Old: Synchronous, blocks on database operations
- New: Fully async, non-blocking I/O

### 2. **Type Safety**
- Old: Everything is strings, error-prone UUID conversions
- New: Type hints everywhere, UUID objects remain UUIDs

### 3. **Automatic Tracking via Context Manager**
- Old: Manual `updateJobStatus()` calls at every status change
- New: Connect to signals once, automatic updates

### 4. **Decoupled File Management**
- Old: File gleaning coupled to status updates
- New: Independent methods for files, or automatic via `glean_job_files()`

### 5. **Signal Integration**
- Old: No signal support, must manually propagate status
- New: Connects to plugin's `statusChanged` signal for automatic updates

### 6. **Cleaner Job Creation**
- Old: String-based parent IDs, manual job number tracking
- New: UUID-based relationships, automatic job numbering

### 7. **Modern Python Patterns**
- Old: Legacy method names (camelCase), no type hints, no context managers
- New: snake_case, type hints, context managers, async iterators

## Migration Path

For existing code using the old handler:

```python
# Old code
handler = CCP4i2DjangoDbHandler()
job = handler.createJob("ctruncate", "My Job", str(parent_uuid))
handler.updateJobStatus(str(job.uuid), finishStatus=CPluginScript.SUCCEEDED)

# New code (minimal changes)
handler = AsyncDatabaseHandler(project_uuid)
job = await handler.create_job("ctruncate", "My Job", parent_uuid)
await handler.update_job_status(job.uuid, models.Job.Status.FINISHED)

# New code (recommended - automatic tracking)
async with handler.track_job(plugin):
    await plugin.execute()
```

## Integration with CPluginScript

The new handler is designed to work seamlessly with the existing database-related attributes in `CPluginScript`:

```python
class CPluginScript:
    # These attributes are automatically set by track_job()
    _dbHandler: Optional[AsyncDatabaseHandler]  # Database handler instance
    _dbProjectId: Optional[uuid.UUID]           # Project UUID
    _dbJobId: Optional[uuid.UUID]               # Job UUID in database
    _dbJobNumber: Optional[str]                 # Job number (e.g., "1.2.3")
```

Example integration:

```python
plugin = CPluginScript(taskName="ctruncate")

# Handler automatically populates these attributes
async with handler.track_job(plugin):
    print(f"Job UUID: {plugin._dbJobId}")
    print(f"Job Number: {plugin._dbJobNumber}")

    await plugin.execute()

    # Plugin can access its database info
    if plugin._dbJobId:
        print(f"Job saved in database: {plugin._dbJobId}")
```

## Error Handling

The new handler uses proper async exception handling:

```python
from django.db import IntegrityError

try:
    async with handler.track_job(plugin):
        await plugin.execute()
except IntegrityError as e:
    logger.error(f"Database integrity error: {e}")
    # Job status automatically updated to FAILED
except Exception as e:
    logger.error(f"Plugin execution failed: {e}")
    # Status preserved in database
```

## Testing

The async nature makes testing cleaner:

```python
import pytest

@pytest.mark.asyncio
async def test_job_tracking():
    handler = AsyncDatabaseHandler(project_uuid=test_project.uuid)
    plugin = CPluginScript(taskName="test_plugin")

    async with handler.track_job(plugin):
        # Simulate execution
        plugin._status = CPluginScript.RUNNING
        await asyncio.sleep(0.1)
        plugin._status = CPluginScript.SUCCEEDED

    # Verify database state
    job = await sync_to_async(models.Job.objects.get)(uuid=plugin._dbJobId)
    assert job.status == models.Job.Status.FINISHED
```
