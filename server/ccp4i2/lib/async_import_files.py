"""
Async file import operations using modern CData introspection.

This module replaces the legacy import_files.py with async/await operations
and uses the new CData metadata system instead of fragile string-based access.
"""

import asyncio
import logging
import shutil
import uuid
from pathlib import Path
from typing import Optional

from asgiref.sync import sync_to_async

# Import CData utilities for legacy field name mapping
from .cdata_utils import get_file_type_from_class

logger = logging.getLogger(f"ccp4x:{__name__}")


async def import_input_files_async(job, plugin, db_handler):
    """
    Import input files for a job using modern async operations.

    This function:
    1. Finds all CDataFile objects in inputData using hierarchical traversal
    2. For files already in database, creates FileUse records
    3. For external files, copies to CCP4_IMPORTED_FILES and registers
    4. Updates container objects with new file locations

    Args:
        job: Django Job model instance
        plugin: CPluginScript instance
        db_handler: AsyncDatabaseHandler instance

    Returns:
        Number of files imported

    Example:
        >>> await import_input_files_async(job, plugin, db_handler)
        3  # Imported 3 files
    """
    files_imported = 0

    # Find all input files using modern hierarchical traversal
    # plugin.container.inputData contains the input file objects
    input_data = plugin.container.inputData if hasattr(plugin.container, 'inputData') else plugin.container
    input_files = input_data.find_all_files()
    logger.info(f"Found {len(input_files)} input file objects")

    for file_obj in input_files:
        try:
            # Check if file already has database ID (previously registered)
            if hasattr(file_obj, 'dbFileId'):
                db_file_id_obj = file_obj.dbFileId
                if hasattr(db_file_id_obj, 'isSet') and db_file_id_obj.isSet():
                    file_uuid_str = str(db_file_id_obj)
                    if len(file_uuid_str.strip()) > 0:
                        # File already in database, just create FileUse
                        await db_handler.register_input_file(
                            job_uuid=job.uuid,
                            file_uuid=uuid.UUID(file_uuid_str),
                            param_name=file_obj.objectName(),
                        )
                        logger.info(f"Registered existing file for {file_obj.objectName()}")
                        files_imported += 1
                        continue

            # Check if file has a baseName set (external file to import)
            if hasattr(file_obj, 'baseName'):
                base_name_obj = file_obj.baseName
                if hasattr(base_name_obj, 'isSet') and base_name_obj.isSet():
                    base_name = str(base_name_obj).strip()
                    if len(base_name) > 0:
                        # Import external file
                        await import_external_file_async(job, file_obj, db_handler)
                        files_imported += 1

        except Exception as e:
            logger.exception(f"Error importing file {file_obj.objectName()}: {e}")

    # Save updated parameters to input_params.xml
    # This preserves the dbFileId, relPath, baseName changes made during import
    if files_imported > 0:
        input_params_file = job.directory / "input_params.xml"

        # DEBUG: Check what's in the container before saving
        logger.info(f"[DEBUG import_input_files_async] About to save input_params.xml for job {job.number}")
        logger.info(f"[DEBUG import_input_files_async] Plugin container type: {type(plugin.container).__name__}")
        logger.info(f"[DEBUG import_input_files_async] Has inputData: {hasattr(plugin.container, 'inputData')}")
        if hasattr(plugin.container, 'inputData'):
            logger.info(f"[DEBUG import_input_files_async] inputData type: {type(plugin.container.inputData).__name__}")
            input_data_children = [c.objectName() for c in plugin.container.inputData.children() if hasattr(c, 'objectName')]
            logger.info(f"[DEBUG import_input_files_async] inputData children: {input_data_children}")

        logger.info(f"Saving updated parameters to {input_params_file}")
        error = await sync_to_async(plugin.saveDataToXml)(str(input_params_file))
        if error and hasattr(error, 'hasError') and error.hasError():
            logger.error(f"Failed to save parameters: {error}")
        else:
            logger.info("Successfully saved updated input_params.xml with imported file paths")

    logger.info(f"Imported {files_imported} input files")
    return files_imported


async def import_external_file_async(job, file_obj, db_handler):
    """
    Import an external file to CCP4_IMPORTED_FILES using modern CData introspection.

    This function implements checksum-based deduplication for imported files:
    - If a file with the same checksum and type already exists, reuse it
    - Otherwise, copy the file and create new File/FileImport records

    Args:
        job: Django Job model instance
        file_obj: CDataFile object
        db_handler: AsyncDatabaseHandler instance
    """
    # Extract metadata using modern CData system
    core_metadata = file_obj.to_metadata_dict()

    # Map to legacy field names for database operations
    metadata = {
        'name': file_obj.objectName(),
        'file_type': get_file_type_from_class(file_obj),
        'content_flag': core_metadata.get('contentFlag'),
        'sub_type': core_metadata.get('subType'),
        'annotation': core_metadata.get('annotation'),
        'gui_label': file_obj.get_qualifier('guiLabel', ''),
    }

    # Get source file path
    source_path = await get_source_file_path(job, file_obj)
    if source_path is None or not source_path.exists():
        logger.warning(f"Source file does not exist: {source_path}")
        return

    # Calculate checksum from SOURCE file (before copying)
    # This allows us to check for duplicates without copying first
    checksum = None
    try:
        import hashlib
        def compute_checksum_sync(path):
            """Synchronous checksum computation to be wrapped with sync_to_async."""
            md5_hash = hashlib.md5()
            with open(path, 'rb') as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    md5_hash.update(chunk)
            return md5_hash.hexdigest()
        checksum = await sync_to_async(compute_checksum_sync)(source_path)
        logger.debug(f"Calculated checksum for {source_path.name}: {checksum[:8]}...")
    except Exception as e:
        logger.warning(f"Could not calculate checksum for {source_path}: {e}")

    # Check if file with same checksum and type already exists
    existing_file = None
    if checksum and hasattr(db_handler, 'find_imported_file_by_checksum'):
        existing_file = await db_handler.find_imported_file_by_checksum(
            checksum=checksum,
            file_type=metadata['file_type']
        )

    if existing_file:
        # Reuse existing file - no copy needed!
        logger.info(
            f"Deduplicating import: reusing existing file {existing_file.name} "
            f"(checksum {checksum[:8]}...) instead of copying {source_path.name}"
        )
        file_record = existing_file
        dest_path = existing_file.path

        # Create FileUse record to track that this job uses the existing file
        await db_handler.register_input_file(
            job_uuid=job.uuid,
            file_uuid=file_record.uuid,
            param_name=metadata['name'],
        )
    else:
        # No duplicate found - proceed with normal import
        # Determine destination path in CCP4_IMPORTED_FILES
        import_dir = Path(job.project.directory) / "CCP4_IMPORTED_FILES"

        # Check if source file is already in CCP4_IMPORTED_FILES
        # This happens when files are pre-processed (e.g., MTZ splitting)
        # In this case, we can skip the copy step
        try:
            source_resolved = source_path.resolve()
            import_dir_resolved = import_dir.resolve()
            if source_resolved.parent == import_dir_resolved:
                # Source is already in CCP4_IMPORTED_FILES - no copy needed
                dest_path = source_resolved
                logger.info(
                    f"Source file {source_path.name} already in CCP4_IMPORTED_FILES, "
                    f"skipping copy"
                )
            else:
                # Source is external - need to copy
                dest_path = import_dir / source_path.name

                # Ensure unique filename
                dest_path = await ensure_unique_path(dest_path)

                # Ensure import directory exists
                await sync_to_async(import_dir.mkdir)(parents=True, exist_ok=True)

                # Copy file asynchronously
                await async_copy_file(source_path, dest_path)
                logger.info(f"Copied {source_path} to {dest_path}")
        except (OSError, ValueError) as e:
            # If path resolution fails, fall back to copying
            logger.debug(f"Path resolution failed, proceeding with copy: {e}")
            dest_path = import_dir / source_path.name
            dest_path = await ensure_unique_path(dest_path)
            await sync_to_async(import_dir.mkdir)(parents=True, exist_ok=True)
            await async_copy_file(source_path, dest_path)
            logger.info(f"Copied {source_path} to {dest_path}")

        # Register in database
        file_record = await db_handler.register_imported_file(
            job_uuid=job.uuid,
            file_path=dest_path,
            file_type=metadata['file_type'],
            param_name=metadata['name'],
            source_path=source_path,
            annotation=metadata.get('annotation', f"Imported from {source_path.name}"),
            checksum=checksum,
        )

    # Update container to reflect new location
    if hasattr(file_obj, 'dbFileId'):
        file_obj.dbFileId.set(str(file_record.uuid))

    if hasattr(file_obj, 'relPath'):
        file_obj.relPath.set("CCP4_IMPORTED_FILES")

    if hasattr(file_obj, 'baseName'):
        file_obj.baseName.set(dest_path.name)

    if hasattr(file_obj, 'project'):
        file_obj.project.set(str(job.project.uuid))

    # Reload file with new location
    if hasattr(file_obj, 'loadFile'):
        await sync_to_async(file_obj.loadFile)()

    if hasattr(file_obj, 'setContentFlag'):
        await sync_to_async(file_obj.setContentFlag)()

    logger.info(f"Successfully imported {file_obj.objectName()}")


async def get_source_file_path(job, file_obj) -> Optional[Path]:
    """
    Get the source file path from a CDataFile object.

    Args:
        job: Django Job model instance
        file_obj: CDataFile object

    Returns:
        Path to source file, or None if cannot be determined
    """
    base_name = None
    rel_path = None

    # Extract baseName
    if hasattr(file_obj, 'baseName'):
        base_name_obj = file_obj.baseName
        if hasattr(base_name_obj, 'isSet') and base_name_obj.isSet():
            base_name = str(base_name_obj).strip()

    # Extract relPath
    if hasattr(file_obj, 'relPath'):
        rel_path_obj = file_obj.relPath
        if hasattr(rel_path_obj, 'isSet') and rel_path_obj.isSet():
            rel_path = str(rel_path_obj).strip()

    if not base_name:
        return None

    # Try different path constructions
    if rel_path:
        # Try relPath / baseName
        source_path = Path(rel_path) / base_name
        if source_path.exists():
            return source_path

        # Try project_dir / relPath / baseName
        source_path = Path(job.project.directory) / rel_path / base_name
        if source_path.exists():
            return source_path

    # Try just baseName (current directory or absolute)
    source_path = Path(base_name)
    if source_path.exists():
        return source_path

    return None


async def ensure_unique_path(dest_path: Path) -> Path:
    """
    Ensure destination path is unique by appending _1, _2, etc if needed.

    Args:
        dest_path: Desired destination path

    Returns:
        Unique path that doesn't exist
    """
    def path_exists(p):
        return p.exists()

    if not await sync_to_async(path_exists)(dest_path):
        return dest_path

    # Path exists, find unique name
    counter = 1
    while True:
        stem = dest_path.stem
        suffix = dest_path.suffix
        parent = dest_path.parent

        new_path = parent / f"{stem}_{counter}{suffix}"
        if not await sync_to_async(path_exists)(new_path):
            return new_path

        counter += 1

        # Safety check
        if counter > 1000:
            raise RuntimeError(f"Could not find unique path for {dest_path}")


async def async_copy_file(source: Path, dest: Path):
    """
    Copy a file asynchronously using asyncio.

    Args:
        source: Source file path
        dest: Destination file path
    """
    # Use sync_to_async to avoid blocking
    await sync_to_async(shutil.copyfile)(str(source), str(dest))


async def save_params_after_import(plugin, job):
    """
    Save plugin parameters after importing files.

    This updates the job's parameter file to reflect the new file locations.

    Args:
        plugin: CPluginScript instance
        job: Django Job model instance
    """
    # Import save_params_for_job from modern utilities
    from .utils.parameters.save_params import save_params_for_job

    # DEBUG: Check what's in the container before saving
    logger.info(f"[DEBUG save_params_after_import] About to save input_params.xml for job {job.number}")
    logger.info(f"[DEBUG save_params_after_import] Plugin container type: {type(plugin.container).__name__}")
    logger.info(f"[DEBUG save_params_after_import] Has inputData: {hasattr(plugin.container, 'inputData')}")
    if hasattr(plugin.container, 'inputData'):
        logger.info(f"[DEBUG save_params_after_import] inputData type: {type(plugin.container.inputData).__name__}")
        logger.info(f"[DEBUG save_params_after_import] inputData children: {[c.objectName() for c in plugin.container.inputData.children() if hasattr(c, 'objectName')]}")

    # Run synchronously (it's a quick operation)
    await sync_to_async(save_params_for_job)(plugin, job, mode="JOB_INPUT")
    logger.info(f"Saved parameters for job {job.number}")


# Legacy compatibility wrapper
async def import_files_async_compat(job, plugin):
    """
    Legacy compatibility wrapper that doesn't require db_handler.

    This is for gradual migration from the old import_files.py.

    Args:
        job: Django Job model instance
        plugin: CPluginScript instance

    Returns:
        Updated plugin instance
    """
    from ..db.async_db_handler import AsyncDatabaseHandler

    # Create temporary handler
    handler = AsyncDatabaseHandler(project_uuid=job.project.uuid)

    # Import files
    await import_input_files_async(job, plugin, handler)

    # Save parameters
    await save_params_after_import(plugin, job)

    return plugin
