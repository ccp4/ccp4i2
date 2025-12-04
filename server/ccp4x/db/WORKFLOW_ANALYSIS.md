# CCP4i2 Database-Backed Execution Workflow Analysis

## Overview

This document analyzes the current Django-based job tracking system and identifies opportunities to modernize it using the new Qt-free CData architecture.

## Current Workflow (run_job.py)

The existing job execution follows these stages:

1. **Database Setup** - Initialize `CCP4i2DjangoDbHandler`
2. **Application Instance** - Create Qt event loop (`QtCore.QEventLoop`)
3. **Plugin Retrieval** - Load plugin and set database attributes
4. **Set Output Paths** - Pre-assign output file paths based on job directory
5. **Import Input Files** - Copy external files to `CCP4_IMPORTED_FILES/`, create database records
6. **Execute Plugin** - Run `the_plugin.process()` in Qt event loop
7. **Glean Output Files** - Extract file metadata from container, create database records
8. **Status Update** - Mark job finished/failed

```python
def run_job(jobId: str):
    db_handler = setup_db_handler(new_job)
    application_inst = setup_application_instance(new_job)  # QtCore.QEventLoop
    the_plugin = retrieve_plugin(new_job, application_inst, db_handler)

    new_job.status = models.Job.Status.RUNNING
    new_job.save()

    set_output_file_names(the_plugin.container, projectId=str(new_job.project.uuid), ...)
    _import_files(new_job, the_plugin)

    executePlugin(the_plugin, new_job)  # Runs in Qt event loop

    # On finish signal: glean files, update status
```

## File Import Workflow (import_files.py)

**Purpose**: Copy external input files to project's `CCP4_IMPORTED_FILES/` directory and track in database

```python
def _process_input(theJob, plugin, input: CDataFile):
    # Check if file already has database ID
    if input.dbFileId is not None:
        theFile = models.File.objects.get(uuid=input.dbFileId)
    else:
        # Copy file to CCP4_IMPORTED_FILES
        sourceFilePath = Path(str(input.relPath)) / str(input.baseName)
        destFilePath = Path(project.directory) / "CCP4_IMPORTED_FILES" / sourceFilePath.name
        shutil.copyfile(sourceFilePath, destFilePath)

        # Extract metadata using LEGACY API
        file_mime_type = input.qualifiers("mimeTypeName")  # ❌ String-based access
        annotation = str(input.annotation)
        job_param_name = extract_from_first_bracketed(input.objectPath())

        # Create database record
        theFile = models.File(
            name=str(destFilePath.name),
            annotation=annotation,
            type=file_type_object,
            job=theJob,
            job_param_name=job_param_name,
            directory=2,  # IMPORT_DIR
        )
        theFile.save()

        # Update container to point to imported file
        input.dbFileId.set(theFile.uuid)
        input.relPath.set("CCP4_IMPORTED_FILES")
        input.baseName.set(destFilePath.name)

    # Create FileUse record (IN role)
    models.FileUse(file=theFile, job=theJob, role=models.FileUse.Role.IN, ...).save()
```

**Legacy Approach Issues**:
1. ❌ Uses `.qualifiers("mimeTypeName")` - string-based, no type safety
2. ❌ Uses `getattr(item, "annotation", None)` - fragile attribute access
3. ❌ Manually extracts `job_param_name` from `.objectPath()` string parsing
4. ❌ Direct attribute mutation: `.dbFileId.set()`, `.relPath.set()`

## File Gleaning Workflow (glean_job_files.py)

**Purpose**: Extract output file metadata from container after job execution and register in database

```python
def glean_job_files(jobId, container, roleList=[0, 1]):
    job = models.Job.objects.get(uuid=jobId)

    # Find all CDataFile objects in output container using recursive search
    outputs = find_objects(
        container.outputData,
        lambda a: isinstance(a, CDataFile),
        recursive=True
    )

    # Process each output file
    for item in outputs:
        if item.exists():
            create_new_file(job, item)

def create_new_file(job, item: CDataFile):
    # Extract metadata using LEGACY API
    file_type = item.qualifiers("mimeTypeName")  # ❌ String access
    sub_type = getattr(item, "subType", None)     # ❌ getattr fallback
    content = getattr(item, "contentFlag", None)   # ❌ getattr fallback
    annotation = getattr(item, "annotation", None) # ❌ getattr fallback
    name = str(getattr(item, "baseName", None))   # ❌ getattr fallback

    try:
        sub_type = int(sub_type)  # ❌ Manual type conversion
    except AttributeError:
        sub_type = None

    try:
        content = int(content)  # ❌ Manual type conversion
    except AttributeError:
        content = None

    job_param_name = extract_from_first_bracketed(item.objectPath())

    # Create database record
    the_file = models.File(
        name=name,
        annotation=annotation,
        type=file_type_object,
        sub_type=sub_type,
        content=content,
        job=job,
        job_param_name=job_param_name,
        directory=1,  # JOB_DIR
    )
    the_file.save()

    # Link back to container
    item.dbFileId.set(str(the_file.uuid))
```

**Legacy Approach Issues**:
1. ❌ Uses `.qualifiers("mimeTypeName")` instead of metadata system
2. ❌ Extensive use of `getattr()` with fallbacks - no type safety
3. ❌ Manual type conversions with try/except blocks
4. ❌ String parsing of `objectPath()` to extract parameter name
5. ❌ No use of CData's introspection capabilities

## Output File Path Setting (set_output_file_names.py)

**Purpose**: Pre-assign file paths for output files before execution

```python
def set_output_file_names(container, projectId, jobNumber, force=True):
    relPath = Path("CCP4_JOBS").joinpath(*[f"job_{n}" for n in jobNumber.split(".")])

    dataList = container.outputData.dataOrder()  # ❌ Legacy CContainer method
    for objectName in dataList:
        dobj = container.outputData.find(objectName)  # ❌ String-based find
        if isinstance(dobj, CDataFile) and (force or not dobj.isSet()):
            dobj.setOutputPath(jobName=jobName, projectId=projectId, relPath=str(relPath))
```

**Legacy Approach Issues**:
1. ❌ Uses `.dataOrder()` and `.find()` - legacy CContainer API
2. ❌ String-based object lookup
3. ❌ No use of hierarchical object traversal

## Opportunities with New CData System

### 1. **Type-Safe Metadata Access**

**Current (Legacy)**:
```python
file_type = item.qualifiers("mimeTypeName")  # Returns string or None
sub_type = getattr(item, "subType", None)
```

**New CData**:
```python
# Access through metadata system
qualifiers = item.get_merged_metadata('qualifiers')
file_type = qualifiers.get('mimeTypeName')

# Direct attribute access with type safety
if hasattr(item, 'subType'):
    sub_type = item.subType.value  # Type-safe access
```

### 2. **Hierarchical Traversal**

**Current (Legacy)**:
```python
# Recursive search with lambda predicates
outputs = find_objects(
    container.outputData,
    lambda a: isinstance(a, CDataFile),
    recursive=True
)
```

**New CData**:
```python
# Use hierarchy system
def find_all_files(parent_obj):
    """Recursively find all CDataFile objects in hierarchy"""
    files = []
    for child_name in parent_obj.childNames():
        child = getattr(parent_obj, child_name)
        if isinstance(child, CDataFile):
            files.append(child)
        # Recurse into nested objects
        if isinstance(child, HierarchicalObject):
            files.extend(find_all_files(child))
    return files

# Or use object_path() for debugging
for child in find_all_files(container.outputData):
    print(f"Found file: {child.object_path()}")  # e.g., "outputData.HKLOUT"
```

### 3. **Parameter Name Extraction**

**Current (Legacy)**:
```python
# String parsing with regex
def extract_from_first_bracketed(path: str) -> str:
    parts = path.split(".")
    for i, part in enumerate(parts):
        if re.search(r"\[.*\]", part):
            return ".".join(parts[i:])
    return parts[-1]

job_param_name = extract_from_first_bracketed(item.objectPath())
```

**New CData**:
```python
# Direct attribute access
job_param_name = item.name  # Object's name is its parameter name
# Or use object_path() and parse more reliably
full_path = item.object_path()  # e.g., "container.outputData.HKLOUT"
param_name = full_path.split('.')[-1]  # "HKLOUT"
```

### 4. **Value State Checking**

**Current (Legacy)**:
```python
if item.exists():
    process_file(item)
elif unSetMissingFiles:
    item.unSet()
```

**New CData (Same, but more robust)**:
```python
if item.exists():
    process_file(item)
elif unSetMissingFiles:
    item.unSet()

# But can also check value state
if item.getValueState() == ValueState.EXPLICITLY_SET:
    # User set this value, don't override
    pass
```

### 5. **Smart Attribute Inspection**

**Current (Legacy)**:
```python
# Manual type checking and fallbacks
sub_type = getattr(item, "subType", None)
try:
    sub_type = int(sub_type)
except (AttributeError, TypeError):
    sub_type = None
```

**New CData**:
```python
# Use metadata to know what attributes exist
attributes = item.get_merged_metadata('attributes')
if 'subType' in attributes:
    sub_type = item.subType.value if item.subType.isSet() else None
```

### 6. **Async-Compatible File Operations**

**Current (Legacy)**:
```python
# Synchronous file operations
shutil.copyfile(sourceFilePath, destFilePath)
theFile.save()
```

**New CData (Can be made async)**:
```python
# Async file operations
await async_copy_file(source_path, dest_path)
await sync_to_async(the_file.save)()
```

## Proposed Modern Workflow

### Stage 1: Async Job Execution

```python
async def run_job_async(job_uuid: uuid.UUID):
    """Modern async job execution with database tracking"""
    job = await sync_to_async(models.Job.objects.get)(uuid=job_uuid)

    # Create database handler
    db_handler = AsyncDatabaseHandler(project_uuid=job.project.uuid)

    # Create plugin instance
    plugin = await create_plugin_for_job(job)

    # Track job lifecycle automatically
    async with db_handler.track_job(plugin):
        # Import files (async)
        await import_input_files_async(job, plugin, db_handler)

        # Execute plugin (async)
        result = await plugin.execute()

        # Files automatically gleaned on completion by track_job context manager

    return result
```

### Stage 2: Modern File Import

```python
async def import_input_files_async(job, plugin, db_handler):
    """Import input files using new CData introspection"""
    # Find all CDataFile objects in inputData
    input_files = find_all_files(plugin.inputData)

    for file_obj in input_files:
        # Check if already in database
        if hasattr(file_obj, 'dbFileId') and file_obj.dbFileId.isSet():
            file_uuid = uuid.UUID(str(file_obj.dbFileId))
            await db_handler.register_input_file(
                job_uuid=job.uuid,
                file_uuid=file_uuid,
                param_name=file_obj.name
            )
        else:
            # Import external file
            await import_external_file(job, file_obj, db_handler)

async def import_external_file(job, file_obj, db_handler):
    """Copy external file to project directory and register"""
    # Extract metadata using new CData system
    qualifiers = file_obj.get_merged_metadata('qualifiers')
    file_type = qualifiers.get('mimeTypeName', 'unknown')
    annotation = str(file_obj.annotation) if hasattr(file_obj, 'annotation') else ""

    source_path = Path(str(file_obj.relPath)) / str(file_obj.baseName)
    dest_path = Path(job.project.directory) / "CCP4_IMPORTED_FILES" / source_path.name

    # Async file copy
    await async_copy_file(source_path, dest_path)

    # Register in database
    file_record = await db_handler.register_imported_file(
        job_uuid=job.uuid,
        file_path=dest_path,
        file_type=file_type,
        param_name=file_obj.name,
        annotation=annotation,
        source_path=source_path,
    )

    # Update container
    file_obj.dbFileId.set(str(file_record.uuid))
    file_obj.relPath.set("CCP4_IMPORTED_FILES")
    file_obj.baseName.set(dest_path.name)
```

### Stage 3: Modern File Gleaning

```python
async def glean_job_files_async(job_uuid, output_container):
    """Extract output files using new CData introspection"""
    files_created = []

    # Find all CDataFile objects
    output_files = find_all_files(output_container)

    for file_obj in output_files:
        # Check if file exists on disk
        if not file_obj.exists():
            continue

        # Extract metadata using new CData system
        qualifiers = file_obj.get_merged_metadata('qualifiers')
        attributes = file_obj.get_merged_metadata('attributes')

        file_info = {
            'file_path': Path(str(file_obj)),
            'file_type': qualifiers.get('mimeTypeName', 'unknown'),
            'param_name': file_obj.name,
            'annotation': qualifiers.get('guiLabel', ''),
        }

        # Extract optional attributes if they exist
        if 'subType' in attributes and hasattr(file_obj, 'subType'):
            if file_obj.subType.isSet():
                file_info['sub_type'] = file_obj.subType.value

        if 'contentFlag' in attributes and hasattr(file_obj, 'contentFlag'):
            if file_obj.contentFlag.isSet():
                file_info['content_flag'] = file_obj.contentFlag.value

        # Register in database
        file_record = await db_handler.register_output_file(
            job_uuid=job_uuid,
            **file_info
        )

        # Link back to container
        file_obj.dbFileId.set(str(file_record.uuid))

        files_created.append(file_record)

    return files_created
```

### Stage 4: Performance Indicator Gleaning

```python
async def glean_performance_indicators_async(output_container, job_uuid):
    """Extract KPIs using new CData introspection"""
    # Find all CPerformanceIndicator objects
    kpis = find_objects_by_type(output_container, CPerformanceIndicator)

    for kpi in kpis:
        # Use childNames() to iterate through KPI values
        for param_name in kpi.childNames():
            value_obj = getattr(kpi, param_name)

            # Use type information from CData
            if isinstance(value_obj, CFloat) and value_obj.isSet():
                await db_handler.register_job_float_value(
                    job_uuid=job_uuid,
                    key=param_name,
                    value=value_obj.value
                )
            elif isinstance(value_obj, CString) and value_obj.isSet():
                await db_handler.register_job_char_value(
                    job_uuid=job_uuid,
                    key=param_name,
                    value=str(value_obj)
                )
```

## Key Improvements Summary

| Aspect | Legacy Approach | New CData Approach |
|--------|----------------|-------------------|
| Metadata Access | `.qualifiers("key")` string access | `.get_merged_metadata('qualifiers')` type-safe |
| Attribute Access | `getattr()` with fallbacks | Direct access with `hasattr()` + type checks |
| Type Safety | Manual `int()` conversions with try/except | `.value` property returns correct type |
| Object Traversal | `find_objects()` with lambda | `.childNames()` + hierarchical traversal |
| Parameter Names | Regex parsing of `.objectPath()` | `.name` attribute or `object_path().split('.')` |
| Value State | Basic `.isSet()` check | `.getValueState()` with state tracking |
| Execution Model | Qt event loop (sync) | Async/await |
| Signal System | Qt signals | Modern Python Signal[T] with async |
| Status Tracking | Manual `updateJobStatus()` calls | Automatic via signal connections |
| Error Handling | Try/except around everything | Type-safe access reduces errors |

## Next Steps

1. **Create modern utility functions**:
   - `find_all_files(container)` - Use CData hierarchy traversal
   - `extract_file_metadata(file_obj)` - Use metadata system
   - `find_objects_by_type(container, type_class)` - Type-safe search

2. **Extend AsyncDatabaseHandler**:
   - `register_imported_file()` - For input file tracking
   - `register_job_float_value()` / `register_job_char_value()` - For KPIs
   - Better integration with CData's signal system

3. **Create async job runner**:
   - Replace Qt event loop with async/await
   - Use `track_job()` context manager for automatic tracking
   - Support nested job execution with proper hierarchy

4. **Update CDataFile classes**:
   - Ensure all have `mimeTypeName` in qualifiers
   - Standardize `contentFlag` and `subType` attributes
   - Add helper methods for database interaction

5. **Add tests**:
   - Test modern file gleaning with new CData
   - Test async job execution
   - Test hierarchical job creation
   - Integration tests with real Django database
