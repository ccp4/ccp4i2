"""
Modern async database handler for CCP4i2 plugin execution tracking.

This handler provides a clean API for integrating CPluginScript with the Django database,
using modern Python patterns like async/await, context managers, and type hints.
"""

import datetime
import logging
import uuid
from pathlib import Path
from typing import Optional, Dict, Any, List
from contextlib import asynccontextmanager

from asgiref.sync import sync_to_async
from django.db import transaction
from django.utils import timezone

from core.CCP4PluginScript import CPluginScript

# Import using Django's registered app name to avoid app registry errors
from ccp4x.db import models

logger = logging.getLogger(__name__)


class AsyncDatabaseHandler:
    """
    Modern async database handler for CCP4i2 plugin execution.

    This handler provides:
    - Async/await interface for non-blocking database operations
    - Automatic job lifecycle tracking
    - Signal integration for status updates
    - Context manager support for transactions
    - Type-safe API with modern Python patterns

    Example Usage:
        db_handler = AsyncDatabaseHandler(project_uuid)

        # Create job in database
        job = await db_handler.create_job(
            task_name="ctruncate",
            title="Convert intensities to amplitudes",
            parent_job_uuid=parent_uuid
        )

        # Track job lifecycle automatically
        async with db_handler.track_job(plugin_instance):
            await plugin_instance.execute()

        # Job status, files, and metadata are automatically updated
    """

    def __init__(self, project_uuid: uuid.UUID):
        """
        Initialize database handler for a specific project.

        Args:
            project_uuid: UUID of the project to track jobs for
        """
        self.project_uuid = project_uuid
        self._project: Optional[models.Project] = None

    @property
    def projectId(self) -> uuid.UUID:
        """Legacy compatibility: Get project UUID as projectId."""
        return self.project_uuid

    @property
    def projectName(self) -> str:
        """Legacy compatibility: Get project name (requires sync access)."""
        if self._project is None:
            # Synchronously fetch project if not cached
            # This is for legacy compatibility only - async code should use get_project()
            import asyncio
            try:
                loop = asyncio.get_event_loop()
                if loop.is_running():
                    # In async context, can't block - return a default
                    return "unknown"
                else:
                    # Safe to run sync
                    return loop.run_until_complete(self.get_project()).name
            except RuntimeError:
                # No event loop - use sync
                from asgiref.sync import async_to_sync
                return async_to_sync(self.get_project)().name
        return self._project.name

    @property
    def jobId(self) -> Optional[uuid.UUID]:
        """Legacy compatibility: Get current job ID."""
        # This would need to be tracked separately if needed
        return None

    @property
    def jobNumber(self) -> Optional[str]:
        """Legacy compatibility: Get current job number."""
        # This would need to be tracked separately if needed
        return None

    async def get_project(self) -> models.Project:
        """Get the project instance, cached after first retrieval."""
        if self._project is None:
            self._project = await sync_to_async(
                models.Project.objects.get
            )(uuid=self.project_uuid)
        return self._project

    async def create_job(
        self,
        task_name: str,
        title: Optional[str] = None,
        parent_job_uuid: Optional[uuid.UUID] = None,
        job_number: Optional[str] = None,
    ) -> models.Job:
        """
        Create a new job in the database.

        Args:
            task_name: Name of the plugin/task (e.g., "ctruncate")
            title: Human-readable job title
            parent_job_uuid: UUID of parent job for nested execution
            job_number: Optional explicit job number (e.g., "1" or "1.2")

        Returns:
            Created Job instance

        Example:
            job = await handler.create_job(
                task_name="refmac",
                title="Refinement round 1",
                parent_job_uuid=parent.uuid
            )
        """
        project = await self.get_project()

        @sync_to_async
        def _create_job():
            with transaction.atomic():
                # Determine parent job if specified
                parent_job = None
                if parent_job_uuid:
                    parent_job = models.Job.objects.get(uuid=parent_job_uuid)

                # Auto-generate job number if not provided
                if job_number is None:
                    if parent_job:
                        # Get highest child number for this parent
                        siblings = models.Job.objects.filter(parent=parent_job)
                        if siblings.exists():
                            max_num = max(
                                int(s.number.split('.')[-1])
                                for s in siblings
                            )
                            child_num = max_num + 1
                        else:
                            child_num = 1
                        computed_number = f"{parent_job.number}.{child_num}"
                    else:
                        # Top-level job - compute next job number from existing jobs
                        # Don't rely on project.last_job_number which can get out of sync
                        existing_jobs = models.Job.objects.filter(
                            project=project,
                            parent__isnull=True  # Only top-level jobs
                        )
                        if existing_jobs.exists():
                            # Find the maximum job number among top-level jobs
                            max_num = max(
                                int(j.number.split('.')[0])  # Handle "1", "2", etc.
                                for j in existing_jobs
                            )
                            next_num = max_num + 1
                        else:
                            next_num = 1

                        computed_number = str(next_num)

                        # Update last_job_number to keep it in sync (for legacy compat)
                        project.last_job_number = next_num
                        project.save()
                else:
                    computed_number = job_number

                # Create the job
                job = models.Job.objects.create(
                    project=project,
                    parent=parent_job,
                    number=computed_number,
                    title=title or task_name,
                    task_name=task_name,
                    status=models.Job.Status.PENDING,
                )

                return job

        return await _create_job()

    def createSubJob(
        self,
        taskName: str,
        parentJobId: str,
        jobNumber: str,
        jobTitle: Optional[str] = None,
    ) -> str:
        """
        Create a sub-job in the database (synchronous wrapper for makePluginObject).

        This method provides a synchronous interface for CPluginScript.makePluginObject()
        to create child job records in the database. It wraps the async create_job method.

        Args:
            taskName: Name of the plugin/task for the sub-job (e.g., "refmac")
            parentJobId: UUID of the parent job (as string, with or without hyphens)
            jobNumber: Sub-job number relative to parent (e.g., "1", "2")
            jobTitle: Optional human-readable title for the sub-job

        Returns:
            str: UUID of the created job (without hyphens)

        Example:
            sub_job_id = handler.createSubJob(
                taskName="refmac",
                parentJobId="abc123...",
                jobNumber="1"
            )
        """
        from asgiref.sync import async_to_sync

        # Normalize parent job UUID
        if '-' not in parentJobId and len(parentJobId) == 32:
            # Add hyphens to UUID if missing
            parent_uuid = uuid.UUID(parentJobId)
        else:
            parent_uuid = uuid.UUID(parentJobId)

        # Create job using async method (sync wrapper)
        job = async_to_sync(self.create_job)(
            task_name=taskName,
            title=jobTitle or f"{taskName} sub-job {jobNumber}",
            parent_job_uuid=parent_uuid,
            # Let create_job auto-generate the compound job number (e.g., "1.1")
            job_number=None,
        )

        logger.debug(f"[createSubJob] Created sub-job {job.number} for task {taskName}")
        return str(job.uuid).replace("-", "")

    async def update_job_status(
        self,
        job_uuid: uuid.UUID,
        status: int,
        finish_time: Optional[datetime.datetime] = None,
    ) -> None:
        """
        Update job status in the database.

        Args:
            job_uuid: UUID of job to update
            status: New status (models.Job.Status enum value)
            finish_time: Optional finish time (auto-set for FINISHED status)
        """
        @sync_to_async
        def _update():
            with transaction.atomic():
                # job_uuid is already a UUID object (see parameter type hint)
                job = models.Job.objects.get(uuid=job_uuid)
                job.status = status

                if status == models.Job.Status.FINISHED:
                    job.finish_time = finish_time or timezone.now()

                job.save()

        await _update()

    def updateJobStatus(
        self,
        jobId=None,
        status=None,
        finishStatus=None,
        container=None,
        dbOutputData=None,
    ):
        """
        Legacy sync wrapper for updating job status.

        This method provides compatibility with the legacy CPluginScript interface
        that calls updateJobStatus synchronously after job completion.

        Args:
            jobId: Job UUID as string
            status: Direct status value (optional)
            finishStatus: Plugin finish status (0=SUCCEEDED, 1=FAILED, etc.)
            container: Job container for gleaning output files
            dbOutputData: Deprecated, ignored

        Returns:
            CPluginScript.SUCCEEDED
        """
        from asgiref.sync import async_to_sync

        logger.debug(
            "updateJobStatus (sync wrapper): jobId=%s, status=%s, finishStatus=%s",
            jobId, status, finishStatus
        )

        if jobId is None:
            logger.warning("updateJobStatus called with no jobId")
            return CPluginScript.SUCCEEDED

        try:
            job_uuid = uuid.UUID(jobId) if isinstance(jobId, str) else jobId

            # Convert finishStatus to database status
            if status is None and finishStatus is not None:
                # Map plugin finish status to Job.Status
                # CPluginScript: SUCCEEDED=0, FAILED=1, MARK_TO_DELETE=2, etc.
                if finishStatus == 0:
                    status = models.Job.Status.FINISHED
                elif finishStatus in (1, 2):
                    status = models.Job.Status.FAILED
                else:
                    status = models.Job.Status.UNSATISFACTORY

            if status is not None:
                async_to_sync(self.update_job_status)(job_uuid, status)

                # Glean files if job finished successfully
                if status == models.Job.Status.FINISHED and container is not None:
                    try:
                        async_to_sync(self.glean_job_files)(job_uuid, container)
                    except Exception as e:
                        logger.warning("Failed to glean job files: %s", e)

        except Exception as e:
            logger.error("Error in updateJobStatus: %s", e, exc_info=True)

        return CPluginScript.SUCCEEDED

    async def register_output_file(
        self,
        job_uuid: uuid.UUID,
        file_path: Path,
        file_type: str,
        param_name: str,
        content_flag: Optional[int] = None,
        sub_type: Optional[int] = None,
        annotation: str = "",
    ) -> models.File:
        """
        Register an output file for a job.

        Args:
            job_uuid: UUID of the job that created this file
            file_path: Path to the file (relative to job directory)
            file_type: File type name (e.g., "hklin", "xyzout")
            param_name: Parameter name in plugin (e.g., "HKLOUT")
            content_flag: Content flag from CCP4 data objects
            sub_type: Sub-type from CCP4 data objects
            annotation: Human-readable description

        Returns:
            Created File instance
        """
        @sync_to_async
        def _register():
            with transaction.atomic():
                job = models.Job.objects.get(uuid=job_uuid)

                # Get or create file type
                file_type_obj, _ = models.FileType.objects.get_or_create(
                    name=file_type,
                    defaults={"description": f"File type: {file_type}"}
                )

                # Get or create file record - prevent duplicates
                # Use job + job_param_name as unique identifier since a job
                # shouldn't have two output files with the same param_name
                file_obj, created = models.File.objects.get_or_create(
                    job=job,
                    job_param_name=param_name,
                    defaults={
                        "name": file_path.name,
                        "directory": models.File.Directory.JOB_DIR,
                        "type": file_type_obj,
                        "sub_type": sub_type,
                        "content": content_flag,
                        "annotation": annotation,
                    }
                )

                if not created:
                    # File already exists - update metadata if needed
                    logger.debug(f"File record already exists for {param_name} in job {job_uuid}")
                    # Update fields that might have changed
                    file_obj.name = file_path.name
                    file_obj.type = file_type_obj
                    file_obj.sub_type = sub_type
                    file_obj.content = content_flag
                    if annotation:
                        file_obj.annotation = annotation
                    file_obj.save()

                # Get or create FileUse record - also prevent duplicates
                models.FileUse.objects.get_or_create(
                    file=file_obj,
                    job=job,
                    role=models.FileUse.Role.OUT,
                    job_param_name=param_name,
                )

                return file_obj

        return await _register()

    async def register_input_file(
        self,
        job_uuid: uuid.UUID,
        file_uuid: uuid.UUID,
        param_name: str,
    ) -> None:
        """
        Register that a job uses an existing file as input.

        Args:
            job_uuid: UUID of the job using this file
            file_uuid: UUID of the file being used
            param_name: Parameter name in plugin (e.g., "HKLIN")
        """
        @sync_to_async
        def _register():
            with transaction.atomic():
                job = models.Job.objects.get(uuid=job_uuid)
                file_obj = models.File.objects.get(uuid=file_uuid)

                # Check if FileUse already exists
                if not models.FileUse.objects.filter(
                    file=file_obj,
                    job=job,
                    role=models.FileUse.Role.IN,
                    job_param_name=param_name,
                ).exists():
                    models.FileUse.objects.create(
                        file=file_obj,
                        job=job,
                        role=models.FileUse.Role.IN,
                        job_param_name=param_name,
                    )

        await _register()

    async def find_imported_file_by_checksum(
        self,
        checksum: str,
        file_type: str,
    ) -> Optional[models.File]:
        """
        Find an existing imported file with matching checksum and type.

        This is used to deduplicate imported files - if the same file
        has already been imported, we can reuse it instead of copying again.

        IMPORTANT: This ONLY applies to imported files (FileImport records).
        Generated output files are never deduplicated, even if checksums match.

        Args:
            checksum: MD5 checksum of the file
            file_type: File type name (e.g., "application/CCP4-mtz")

        Returns:
            Existing File instance if found, None otherwise
        """
        @sync_to_async
        def _find():
            try:
                # Look for FileImport with matching checksum
                file_import = models.FileImport.objects.filter(
                    checksum=checksum,
                    file__type__name=file_type,
                    file__directory=models.File.Directory.IMPORT_DIR,
                ).select_related('file', 'file__type', 'file__job__project').first()

                if file_import:
                    # Verify file still exists on disk
                    file_path = file_import.file.path
                    if file_path.exists():
                        logger.info(
                            f"Found existing imported file with checksum {checksum[:8]}... "
                            f"at {file_path}"
                        )
                        return file_import.file
                    else:
                        logger.warning(
                            f"Found FileImport with checksum {checksum[:8]}... "
                            f"but file doesn't exist: {file_path}"
                        )

                return None

            except Exception as e:
                logger.debug(f"Error finding file by checksum: {e}")
                return None

        return await _find()

    async def register_imported_file(
        self,
        job_uuid: uuid.UUID,
        file_path: Path,
        file_type: str,
        param_name: str,
        source_path: Path,
        annotation: str = "",
        checksum: Optional[str] = None,
    ) -> models.File:
        """
        Register an imported file that was copied to CCP4_IMPORTED_FILES.

        Args:
            job_uuid: UUID of the job that imported this file
            file_path: Path to the imported file in CCP4_IMPORTED_FILES
            file_type: File type name (e.g., "application/CCP4-mtz")
            param_name: Parameter name in plugin
            source_path: Original source path before import
            annotation: Human-readable description
            checksum: Optional file checksum

        Returns:
            Created File instance
        """
        @sync_to_async
        def _register():
            with transaction.atomic():
                job = models.Job.objects.get(uuid=job_uuid)

                # Get or create file type
                file_type_obj, _ = models.FileType.objects.get_or_create(
                    name=file_type,
                    defaults={"description": f"File type: {file_type}"}
                )

                # Create file record
                file_obj = models.File.objects.create(
                    name=file_path.name,
                    directory=models.File.Directory.IMPORT_DIR,
                    type=file_type_obj,
                    annotation=annotation,
                    job=job,
                    job_param_name=param_name,
                )

                # Create FileImport record
                models.FileImport.objects.create(
                    file=file_obj,
                    name=str(source_path),
                    checksum=checksum or "",
                    last_modified=timezone.now(),
                )

                # Create FileUse record
                models.FileUse.objects.create(
                    file=file_obj,
                    job=job,
                    role=models.FileUse.Role.IN,
                    job_param_name=param_name,
                )

                return file_obj

        return await _register()

    async def register_job_float_value(
        self,
        job_uuid: uuid.UUID,
        key: str,
        value: float,
        description: Optional[str] = None,
    ) -> None:
        """
        Register a float KPI value for a job.

        Args:
            job_uuid: UUID of the job
            key: KPI key name
            value: Float value
            description: Optional description of the KPI
        """
        @sync_to_async
        def _register():
            with transaction.atomic():
                job = models.Job.objects.get(uuid=job_uuid)

                # Get or create key
                job_value_key, _ = models.JobValueKey.objects.get_or_create(
                    name=key,
                    defaults={"description": description or key}
                )

                # Create or update value
                models.JobFloatValue.objects.update_or_create(
                    job=job,
                    key=job_value_key,
                    defaults={"value": value}
                )

        await _register()

    async def register_job_char_value(
        self,
        job_uuid: uuid.UUID,
        key: str,
        value: str,
        description: Optional[str] = None,
    ) -> None:
        """
        Register a string KPI value for a job.

        Args:
            job_uuid: UUID of the job
            key: KPI key name
            value: String value
            description: Optional description of the KPI
        """
        @sync_to_async
        def _register():
            with transaction.atomic():
                job = models.Job.objects.get(uuid=job_uuid)

                # Get or create key
                job_value_key, _ = models.JobValueKey.objects.get_or_create(
                    name=key,
                    defaults={"description": description or key}
                )

                # Create or update value
                models.JobCharValue.objects.update_or_create(
                    job=job,
                    key=job_value_key,
                    defaults={"value": value}
                )

        await _register()

    async def glean_job_files(
        self,
        job_uuid: uuid.UUID,
        container,
        plugin=None,
    ) -> List[models.File]:
        """
        Extract file information from a job's output container using modern CData utilities.

        This method inspects the job's output container and registers all
        output files in the database.

        Args:
            job_uuid: UUID of the job
            container: CDataContainer with output data
            plugin: Optional CPluginScript instance to provide dbHandler context

        Returns:
            List of created File instances
        """
        files_created = []

        # Use modern hierarchical traversal to find all files
        output_files = container.find_all_files()
        logger.info(f"Found {len(output_files)} files in output container")
        logger.debug(f"[DEBUG glean_job_files] Found {len(output_files)} files in output container")

        for file_obj in output_files:
            # Set temporary plugin reference so file can access dbHandler
            # This is needed because parent chain may be broken during gleaning
            if plugin:
                file_obj._temp_plugin_ref = plugin

            try:
                logger.debug(f"[DEBUG glean_job_files] Processing {file_obj.objectName()}:")
                pass  # DEBUG: print(f"  isSet(): {file_obj.isSet() if hasattr(file_obj, 'isSet') else 'N/A'}")
                pass  # DEBUG: print(f"  exists(): {file_obj.exists() if hasattr(file_obj, 'exists') else 'N/A'}")
                pass  # DEBUG: print(f"  getFullPath(): {file_obj.getFullPath() if hasattr(file_obj, 'getFullPath') else 'N/A'}")

                # Check if file is set (baseName has been assigned)
                if not hasattr(file_obj, 'isSet') or not file_obj.isSet():
                    logger.debug(f"Skipping unset file: {file_obj.objectName()}")
                    pass  # DEBUG: print(f"  -> Skipping: not set")
                    continue

                # Check if file exists on disk
                if not hasattr(file_obj, 'exists') or not file_obj.exists():
                    logger.debug(f"Skipping non-existent file: {file_obj.objectName()} (path: {file_obj.getFullPath() if hasattr(file_obj, 'getFullPath') else 'N/A'})")
                    pass  # DEBUG: print(f"  -> Skipping: doesn't exist")
                    continue
            finally:
                # Clean up temporary reference
                if hasattr(file_obj, '_temp_plugin_ref'):
                    delattr(file_obj, '_temp_plugin_ref')

            pass  # DEBUG: print(f"  -> Will glean this file")

            # Extract metadata using modern CData system
            try:
                # Use core method to extract metadata
                core_metadata = file_obj.to_metadata_dict()

                # Map to legacy field names for database registration
                from ..lib.cdata_utils import get_file_type_from_class

                # Extract param_name from objectPath() for proper list item naming
                # e.g., "outputData.HKLOUT[0]" -> "HKLOUT[0]"
                # e.g., "outputData.FREEROUT" -> "FREEROUT"
                full_path = file_obj.objectPath()

                # IMPORTANT: Only process files that are actually in outputData
                # Pipeline code often assigns inputData references to variables that
                # then get traversed. Skip any files not in the outputData hierarchy.
                if 'outputData' not in full_path:
                    logger.debug(f"Skipping file not in outputData: {full_path}")
                    continue

                if '.outputData.' in full_path:
                    param_name = full_path.split('.outputData.', 1)[1]
                elif 'outputData.' in full_path:
                    param_name = full_path.split('outputData.', 1)[1]
                else:
                    # Fallback to objectName() if path doesn't contain outputData
                    param_name = file_obj.objectName()

                metadata = {
                    'name': param_name,
                    'file_type': get_file_type_from_class(file_obj),
                    'content_flag': core_metadata.get('contentFlag'),
                    'sub_type': core_metadata.get('subType'),
                    'annotation': core_metadata.get('annotation'),
                    'gui_label': file_obj.get_qualifier('guiLabel', ''),
                }

                logger.debug(f"[DEBUG glean_job_files] Extracted metadata for {file_obj.objectName()}:")
                pass  # DEBUG: print(f"  file_type (mimeTypeName): {metadata['file_type']}")
                pass  # DEBUG: print(f"  content_flag: {metadata.get('content_flag', 'NOT SET')}")
                pass  # DEBUG: print(f"  sub_type: {metadata.get('sub_type', 'NOT SET')}")
                pass  # DEBUG: print(f"  gui_label: {metadata.get('gui_label', 'NOT SET')}")

                # Get full file path
                file_path = Path(str(file_obj))

                # Handle CPdbDataFile: rename mmcif files to .cif extension
                from core.CCP4ModelData import CPdbDataFile
                if isinstance(file_obj, CPdbDataFile):
                    print(f"[GLEAN DEBUG] Found CPdbDataFile: {file_obj.objectName()}")
                    print(f"[GLEAN DEBUG]   Current file_path: {file_path}")
                    print(f"[GLEAN DEBUG]   Current suffix: {file_path.suffix}")
                    logger.debug(f"[GLEAN DEBUG] Found CPdbDataFile: {file_obj.objectName()}")
                    logger.debug(f"[GLEAN DEBUG]   Current file_path: {file_path}")
                    logger.debug(f"[GLEAN DEBUG]   Current suffix: {file_path.suffix}")
                    try:
                        print(f"[GLEAN DEBUG]   About to call isMMCIF()...")
                        is_mmcif = file_obj.isMMCIF()
                        print(f"[GLEAN DEBUG]   isMMCIF() returned: {is_mmcif}")
                        logger.debug(f"[GLEAN DEBUG]   isMMCIF() returned: {is_mmcif}")

                        if is_mmcif:
                            print(f"[GLEAN DEBUG]   File IS mmCIF format")
                            # Check if file needs renaming (wrong extension for mmCIF)
                            # ONLY rename .pdb files to .cif to avoid mislabeled PDB files that are actually mmCIF
                            # Do NOT rename .mmcif or .ent files - those are valid extensions
                            current_suffix = file_path.suffix.lower()
                            logger.debug(f"[GLEAN DEBUG]   File is mmCIF, checking suffix: {current_suffix}")

                            if current_suffix == '.pdb':
                                # New path with .cif extension (to fix mislabeled .pdb files that are actually mmCIF)
                                new_file_path = file_path.with_suffix('.cif')
                                logger.debug(f"[GLEAN DEBUG]   Will rename .pdb to .cif: {new_file_path}")

                                # 1. Rename file on disk
                                import shutil
                                shutil.move(str(file_path), str(new_file_path))
                                logger.info(f"✅ Renamed mmCIF file from {file_path.name} to {new_file_path.name}")

                                # 2. Update file_obj baseName (update the CDataFile metadata)
                                file_obj.baseName.set(new_file_path.name)

                                # 3. Update file_path variable for database registration
                                file_path = new_file_path
                            else:
                                logger.debug(f"[GLEAN DEBUG]   Suffix {current_suffix} is already correct (.cif, .mmcif, .ent are all valid), no rename needed")
                    except Exception as e:
                        logger.warning(f"Failed to rename mmCIF file {file_path.name}: {e}")
                        import traceback
                        logger.warning(traceback.format_exc())
                        # Continue with original path if rename fails

                # Register the file in database
                file_record = await self.register_output_file(
                    job_uuid=job_uuid,
                    file_path=file_path,
                    file_type=metadata['file_type'],
                    param_name=metadata['name'],
                    content_flag=metadata.get('content_flag'),
                    sub_type=metadata.get('sub_type'),
                    annotation=metadata.get('annotation', metadata.get('gui_label', '')),
                )

                pass  # DEBUG: print(f"  -> Registered as: {file_record.type.name if file_record.type else 'NO TYPE'}")

                # Link back to container
                if hasattr(file_obj, 'dbFileId'):
                    file_obj.dbFileId.set(str(file_record.uuid))

                files_created.append(file_record)

            except Exception as e:
                logger.exception(f"Error gleaning file {file_obj.objectName()}: {e}")

        return files_created

    async def glean_performance_indicators(
        self,
        job_uuid: uuid.UUID,
        container,
    ) -> int:
        """
        Extract performance indicators (KPIs) from output container.

        Args:
            job_uuid: UUID of the job
            container: CDataContainer with output data

        Returns:
            Number of KPIs extracted
        """
        # Import here to avoid circular dependency
        from ..lib.cdata_utils import extract_kpi_values
        # Will need to import CPerformanceIndicator type
        try:
            from core.CCP4PerformanceData import CPerformanceIndicator
        except ImportError:
            try:
                from core.CCP4PerformanceData import CPerformanceIndicator
            except ImportError:
                logger.warning("Could not import CPerformanceIndicator")
                return 0

        count = 0

        # Find all performance indicator objects
        kpis = container.find_children_by_type(CPerformanceIndicator)

        for kpi in kpis:
            try:
                # Extract all KPI values
                values = extract_kpi_values(kpi)

                # Register each value in database
                for key, value in values.items():
                    if isinstance(value, float):
                        await self.register_job_float_value(
                            job_uuid=job_uuid,
                            key=key,
                            value=value,
                        )
                        count += 1
                    elif isinstance(value, str) and len(value) > 0:
                        await self.register_job_char_value(
                            job_uuid=job_uuid,
                            key=key,
                            value=value,
                        )
                        count += 1

            except Exception as e:
                logger.exception(f"Error gleaning KPIs from {kpi.object_path()}: {e}")

        return count

    @asynccontextmanager
    async def track_job(self, plugin: CPluginScript):
        """
        Context manager for automatic job lifecycle tracking.

        This context manager:
        1. Creates job record in database (if not exists)
        2. Connects to plugin signals to track status changes
        3. Updates database on status transitions
        4. Registers output files on completion

        Usage:
            async with db_handler.track_job(plugin_instance):
                await plugin_instance.execute()
            # Job is automatically tracked and files registered

        Args:
            plugin: CPluginScript instance to track
        """
        # Create job if doesn't exist in database
        if plugin.get_db_job_id() is None:
            parent_uuid = None
            if plugin.parent() and hasattr(plugin.parent(), 'get_db_job_id'):
                parent_uuid = plugin.parent().get_db_job_id()

            job = await self.create_job(
                task_name=plugin.TASKNAME or "unknown",
                title=plugin.name,
                parent_job_uuid=parent_uuid,
            )
            plugin.set_db_job_id(job.uuid)
            plugin.set_db_job_number(job.number)

            # Set project ID for database-aware file handling
            plugin._dbProjectId = self.project_uuid

            # Attach this handler to the plugin so CDataFile can find it
            plugin._db_handler = self

            # Set work directory to job directory for proper file paths
            # This ensures output files are created in the correct job directory
            from pathlib import Path
            job_dir_sync = await sync_to_async(lambda: job.directory)()
            plugin.workDirectory = Path(job_dir_sync)

        # Mark plugin as being tracked by track_job context manager
        # This prevents double-gleaning since track_job handles gleaning for top-level jobs
        # IMPORTANT: This must be set OUTSIDE the if block, because the plugin may already
        # have a job ID set via setDbData() before entering track_job
        plugin._tracked_by_track_job = True

        job_uuid = plugin.get_db_job_id()

        try:
            # Mark as running
            await self.update_job_status(job_uuid, models.Job.Status.RUNNING)

            # Execute the plugin
            yield plugin

            # After execution, update job status based on plugin status
            plugin_status = plugin.get_status()
            status_map = {
                CPluginScript.SUCCEEDED: models.Job.Status.FINISHED,
                CPluginScript.FAILED: models.Job.Status.FAILED,
            }

            if plugin_status in status_map:
                db_status = status_map[plugin_status]
                await self.update_job_status(job_uuid, db_status)
                logger.info(f"Job {job_uuid} status updated to {db_status}")

            # After execution, glean output files and KPIs if finished successfully
            logger.debug(f"[DEBUG track_job] plugin_status = {plugin_status}, SUCCEEDED = {CPluginScript.SUCCEEDED}")
            if plugin_status == CPluginScript.SUCCEEDED:
                logger.debug(f"[DEBUG track_job] Status is SUCCEEDED, gleaning files...")
                output_container = plugin.container.outputData if hasattr(plugin.container, 'outputData') else None
                logger.debug(f"[DEBUG track_job] output_container = {output_container}")
                logger.debug(f"[DEBUG track_job] output_container is not None = {output_container is not None}")
                if output_container is not None:
                    # Pass plugin so file objects can access dbHandler during gleaning
                    files_gleaned = await self.glean_job_files(job_uuid, output_container, plugin=plugin)
                    logger.info(f"Gleaned {len(files_gleaned)} output files")
                    logger.debug(f"[DEBUG track_job] Gleaned {len(files_gleaned)} output files")

                    kpis_gleaned = await self.glean_performance_indicators(job_uuid, output_container)
                    logger.info(f"Gleaned {kpis_gleaned} performance indicators")
                    logger.debug(f"[DEBUG track_job] Gleaned {kpis_gleaned} performance indicators")

                    # Save params.xml with updated dbFileId values
                    if len(files_gleaned) > 0:
                        logger.debug(f"[DEBUG track_job] Saving params.xml after gleaning...")
                        from ..lib.utils.parameters.save_params import save_params_for_job
                        job = await sync_to_async(models.Job.objects.get)(uuid=job_uuid)
                        await sync_to_async(save_params_for_job)(plugin, job, mode="PARAMS")
                        logger.info(f"Saved params.xml with gleaned file IDs")
                        logger.debug(f"[DEBUG track_job] Saved params.xml with {len(files_gleaned)} file IDs")
                else:
                    logger.debug(f"[DEBUG track_job] No output container found!")
            else:
                logger.debug(f"[DEBUG track_job] Status is NOT SUCCEEDED, skipping gleaning")

        finally:
            # Cleanup if needed
            pass

    def run_subjob(self, plugin: CPluginScript, parent_job_uuid: uuid.UUID) -> int:
        """
        Run a subjob plugin with full database tracking (synchronous wrapper).

        This method is intended for pipeline code that creates subjobs via
        makePluginObject(). It handles:
        - Creating the subjob record in database with proper parent FK
        - Setting up database context on the plugin
        - Running process()
        - Gleaning output files and KPIs
        - Updating job status

        Unlike track_job (which is async), this method can be called from
        synchronous pipeline code.

        Args:
            plugin: The subjob plugin instance (already initialized with container)
            parent_job_uuid: UUID of the parent job

        Returns:
            Plugin status code (CPluginScript.SUCCEEDED, FAILED, etc.)
        """
        from asgiref.sync import async_to_sync

        async def _run_subjob_async():
            # Get parent job to determine job number
            parent_job = await sync_to_async(models.Job.objects.get)(uuid=parent_job_uuid)

            # Count existing subjobs to determine the next job number
            existing_subjobs = await sync_to_async(
                lambda: models.Job.objects.filter(parent=parent_job).count()
            )()
            subjob_index = existing_subjobs + 1
            job_number = f"{parent_job.number}.{subjob_index}"

            # Create subjob record in database
            task_name = plugin.TASKNAME or type(plugin).__name__
            job = await self.create_job(
                task_name=task_name,
                title=plugin.name if hasattr(plugin, 'name') else task_name,
                parent_job_uuid=parent_job_uuid,
            )

            # Set database context on plugin
            plugin.set_db_job_id(job.uuid)
            plugin.set_db_job_number(job.number)
            plugin._dbProjectId = self.project_uuid
            plugin._db_handler = self
            plugin._tracked_by_track_job = True  # Mark so process() doesn't double-glean

            logger.info(f"Created subjob {job.number} ({task_name}) under parent {parent_job.number}")

            try:
                # Mark as running
                await self.update_job_status(job.uuid, models.Job.Status.RUNNING)

                # Execute the plugin synchronously (this is what pipelines expect)
                status = plugin.process()

                # Update job status based on result
                status_map = {
                    CPluginScript.SUCCEEDED: models.Job.Status.FINISHED,
                    CPluginScript.FAILED: models.Job.Status.FAILED,
                }
                if status in status_map:
                    await self.update_job_status(job.uuid, status_map[status])

                # Glean output files and KPIs if succeeded
                if status == CPluginScript.SUCCEEDED:
                    output_container = plugin.container.outputData if hasattr(plugin.container, 'outputData') else None
                    if output_container is not None:
                        files_gleaned = await self.glean_job_files(job.uuid, output_container, plugin=plugin)
                        logger.info(f"Subjob {job.number}: Gleaned {len(files_gleaned)} output files")

                        kpis_gleaned = await self.glean_performance_indicators(job.uuid, output_container)
                        logger.info(f"Subjob {job.number}: Gleaned {kpis_gleaned} performance indicators")

                        # Save params.xml with updated dbFileId values
                        if len(files_gleaned) > 0:
                            from ..lib.utils.parameters.save_params import save_params_for_job
                            await sync_to_async(save_params_for_job)(plugin, job, mode="PARAMS")

                return status

            except Exception as e:
                logger.error(f"Subjob {job.number} failed with exception: {e}")
                await self.update_job_status(job.uuid, models.Job.Status.FAILED)
                raise

        return async_to_sync(_run_subjob_async)()

    async def find_file_by_path(
        self,
        file_path: str,
        job_number: str,
        filename: str
    ) -> Optional[models.File]:
        """
        Find an existing file record in the database by job and filename.

        This method is called by CDataFile.setFullPath() to update file attributes
        from the database when a file path is set in a database-aware context.

        Args:
            file_path: Full absolute path to the file
            job_number: Job number (e.g., "1.2")
            filename: Base filename

        Returns:
            File instance if found, None otherwise
        """
        try:
            return await sync_to_async(models.File.objects.get)(
                job__project__uuid=self.project_uuid,
                job__number=job_number,
                name=filename
            )
        except models.File.DoesNotExist:
            return None
        except Exception as e:
            logger.debug(f"Error finding file by path: {e}")
            return None

    async def get_file_path(self, file_uuid: uuid.UUID) -> Optional[str]:
        """
        Get the full file path from a database file UUID.

        This method is called by CDataFile.getFullPath() when dbFileId is set,
        allowing the file path to be retrieved from the database.

        Args:
            file_uuid: UUID of the file record

        Returns:
            Full path to the file, or None if not found
        """
        try:
            file_record = await sync_to_async(models.File.objects.get)(uuid=file_uuid)
            return str(file_record.path) if file_record.path else None
        except models.File.DoesNotExist:
            return None
        except Exception as e:
            logger.debug(f"Error getting file path: {e}")
            return None

    def find_file_by_path_sync(
        self,
        file_path: str,
        job_number: str,
        filename: str
    ) -> Optional[Dict[str, Any]]:
        """
        Synchronous version of find_file_by_path for use by CDataFile.

        Returns a dict with file attributes instead of a model instance
        to avoid exposing Django models to core code.

        Returns:
            Dict with keys: uuid, name, path, relative_path, or None if not found
        """
        try:
            file_record = models.File.objects.get(
                job__project__uuid=self.project_uuid,
                job__number=job_number,
                name=filename
            )
            return {
                'uuid': str(file_record.uuid),
                'name': file_record.name,
                'path': str(file_record.path) if file_record.path else None,
                'relative_path': file_record.relative_path,
            }
        except models.File.DoesNotExist:
            return None
        except Exception as e:
            logger.debug(f"Error finding file by path: {e}")
            return None

    def get_file_path_sync(self, file_uuid: uuid.UUID) -> Optional[str]:
        """
        Synchronous version of get_file_path for use by CDataFile.

        Args:
            file_uuid: UUID of the file record

        Returns:
            Full path to the file, or None if not found or if path doesn't exist
        """
        from pathlib import Path

        try:
            file_record = models.File.objects.get(uuid=file_uuid)
            if not file_record.path:
                return None

            # Check if the file actually exists at the recorded path
            # This prevents returning stale paths from previous test runs or moved files
            file_path = Path(file_record.path)
            if not file_path.exists():
                logger.debug(f"File path in database doesn't exist: {file_record.path}, returning None")
                return None

            return str(file_record.path)
        except models.File.DoesNotExist:
            return None
        except Exception as e:
            logger.debug(f"Error getting file path: {e}")
            return None

    def parse_file_path(self, file_path: str) -> Optional[Dict[str, str]]:
        """
        Parse a full file path into its components: project UUID, job number, filename.

        Used by CDataFile.setFullPath() to extract database context from file paths.

        Args:
            file_path: Full absolute path to the file

        Returns:
            Dict with keys: project_id, rel_path, filename, or None if parsing fails
        """
        file_path_path = Path(file_path)
        project_id_dir_tuples = models.Project.objects.all().values_list('uuid', 'directory')
        # Use next() with default value to avoid StopIteration if no match found
        project_id_dir = next(
            (project_id_dir_tuple for project_id_dir_tuple
             in project_id_dir_tuples if
             str(file_path).startswith(str(project_id_dir_tuple[1]))),
            None  # Return None if no matching project directory found
        )
        if not project_id_dir:
            return None
        project_id = str(project_id_dir[0])
        project_dir = project_id_dir[1]

        if Path(project_dir) / "CCP4_IMPORTED_FILES" is file_path_path.parent:
            # Imported file - no job number
            filename = file_path_path.name
            return {
                'project_id': project_id,
                'rel_path': 'CCP4_IMPORTED_FILES',
                'filename': filename,
            }
        
        job_number_dirs = [(job.number, job.directory) for job in models.Job.objects.filter(
            project__uuid=uuid.UUID(project_id)
        )]
        # Use next() with default value to avoid StopIteration if no match found
        job_number_dir = next(
            (job_number_dir for job_number_dir in job_number_dirs
             if str(file_path_path.parent) == str(job_number_dir[1])),
            None  # Return None if no matching job directory found
        )
        if not job_number_dir:
            return None
        filename = file_path_path.name
        rel_path = str(file_path_path.parent.relative_to(project_dir))
        return {
            'project_id': project_id,
            'rel_path': rel_path,
            'filename': filename,
        }

    def getProjectDirectory(self, project_id: str) -> Optional[str]:
        """
        Get the directory path for a project by its UUID.

        Used by CDataFile.getFullPath() to construct full file paths from
        project/relPath/baseName components.

        This is a synchronous method that safely handles both sync and async contexts.

        Args:
            project_id: UUID of the project (as string)

        Returns:
            Full path to the project directory, or None if not found
        """
        try:
            import uuid as uuid_module
            import asyncio

            # Check if we're in an async context
            try:
                asyncio.get_running_loop()
                is_async = True
            except RuntimeError:
                is_async = False

            project_uuid = uuid_module.UUID(str(project_id))

            if is_async:
                # We're in async context - need to run query in a thread pool
                def _get_project():
                    try:
                        proj = models.Project.objects.get(uuid=project_uuid)
                        return str(proj.directory) if proj.directory else None
                    except models.Project.DoesNotExist:
                        return None

                import concurrent.futures
                with concurrent.futures.ThreadPoolExecutor() as executor:
                    future = executor.submit(_get_project)
                    return future.result(timeout=5.0)
            else:
                # Synchronous context - safe to query directly
                project = models.Project.objects.get(uuid=project_uuid)
                return str(project.directory) if project.directory else None

        except models.Project.DoesNotExist:
            logger.debug(f"Project not found: {project_id}")
            return None
        except Exception as e:
            logger.debug(f"Error getting project directory: {e}")
            return None


def plugin_status_to_job_status(finish_status: int) -> int:
    """
    Convert CPluginScript finish status to Job.Status enum value.

    Args:
        finish_status: CPluginScript status code

    Returns:
        models.Job.Status enum value
    """
    status_map = {
        CPluginScript.SUCCEEDED: models.Job.Status.FINISHED,
        CPluginScript.FAILED: models.Job.Status.FAILED,
        CPluginScript.INTERRUPTED: models.Job.Status.INTERRUPTED,
        CPluginScript.MARK_TO_DELETE: models.Job.Status.TO_DELETE,
        CPluginScript.UNSATISFACTORY: models.Job.Status.UNSATISFACTORY,
    }
    return status_map.get(finish_status, models.Job.Status.FAILED)
