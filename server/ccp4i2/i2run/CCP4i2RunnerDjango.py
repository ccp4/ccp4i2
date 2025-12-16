"""
CCP4i2RunnerDjango - Tidied version with optimized Django ORM usage.

Improvements:
- Eliminated redundant database queries
- Used Django ORM best practices (.select_related, .exists(), etc.)
- Cleaner Path operations
- Better error messages
"""

import logging
from pathlib import Path
from django.conf import settings
from django.utils.text import slugify
from ..cli.i2run.CCP4i2RunnerBase import CCP4i2RunnerBase
from .i2run_components import PluginPopulator
from ..db import models
from ..api import serializers
from ..lib.utils.jobs.create import create_job
from ..lib.utils.parameters.save_params import save_params_for_job

logger = logging.getLogger(f"ccp4i2:{__name__}")


class CCP4i2RunnerDjango(CCP4i2RunnerBase):
    """
    Django-backed i2run implementation.

    Uses Django ORM for project/job management and AsyncDatabaseHandler
    for file operations.
    """

    def __init__(self, the_args=None, command_line=None, parser=None, parent=None):
        """
        Initialize Django runner.

        Note: Validation of args is handled by parent class.
        """
        super().__init__(
            the_args=the_args,
            command_line=command_line,
            parser=parser,
            parent=parent
        )

    def fileForFileUse(
        self,
        projectName=None,
        projectId=None,
        task_name=None,
        jobIndex=None,
        jobParamName=None,
        paramIndex=-1,
    ):
        """
        Resolve file reference from fileUse syntax.

        Returns dictionary with file metadata for database-aware file operations.
        """
        # Validate required parameters
        if projectName is None and projectId is None:
            raise ValueError("Must provide either projectName or projectId")
        if jobParamName is None:
            raise ValueError("jobParamName is required")
        if jobIndex is None:
            raise ValueError("jobIndex is required")

        # Get project (single query)
        if projectId:
            theProject = models.Project.objects.get(uuid=projectId)
        else:
            theProject = models.Project.objects.get(name=projectName)

        # Get jobs with task_name filter if provided (single query with select_related)
        jobs_query = models.Job.objects.filter(project=theProject).select_related('project')
        if task_name is not None:
            jobs_query = jobs_query.filter(task_name=task_name)

        job_count = jobs_query.count()
        logger.info(f"Filtered job list has {job_count} jobs")

        if job_count <= jobIndex:
            raise ValueError(
                f'Requested job index {jobIndex} but only {job_count} jobs with '
                f'task_name="{task_name}" exist in project "{theProject.name}"'
            )

        # Get the specific job (use slicing instead of list conversion)
        theJob = jobs_query[jobIndex]
        logger.info(f"Selected job: {theJob}")

        # Try output files first (single query with select_related)
        output_files = models.File.objects.filter(
            job=theJob,
            job_param_name=jobParamName
        ).select_related('job__project')

        if output_files.exists():
            theFile = output_files[paramIndex]
        else:
            # Try input files (optimized with values_list for ID extraction)
            input_file_ids = models.FileUse.objects.filter(
                job=theJob,
                job_param_name=jobParamName
            ).values_list('file_id', flat=True)

            if not input_file_ids:
                raise ValueError(
                    f"Could not find file for job {theJob} with "
                    f"job_param_name={jobParamName} and index={paramIndex}"
                )

            input_files = models.File.objects.filter(
                id__in=input_file_ids
            ).select_related('job__project')

            if input_files.count() <= paramIndex:
                raise ValueError(
                    f"File index {paramIndex} out of range for "
                    f"job_param_name={jobParamName} (found {input_files.count()} files)"
                )

            theFile = input_files[paramIndex]

        # Build file metadata dictionary
        fileDict = {
            "project": str(theFile.job.project.uuid).replace("-", ""),
            "baseName": theFile.name,
            "dbFileId": str(theFile.uuid).replace("-", ""),
        }

        # Determine relPath based on file directory type
        if theFile.directory == models.File.Directory.IMPORT_DIR:
            fileDict["relPath"] = "CCP4_IMPORTED_FILES"
        else:
            # For job files, extract relative path from job directory
            job_dir = Path(theFile.job.directory)
            if "CCP4_JOBS" in job_dir.parts:
                # Find CCP4_JOBS and take everything after it
                parts = job_dir.parts
                ccp4_jobs_index = parts.index("CCP4_JOBS")
                fileDict["relPath"] = str(Path(*parts[ccp4_jobs_index:]))
            else:
                # Fallback
                fileDict["relPath"] = str(job_dir.name)

        logger.info(f"Resolved file: {fileDict}")
        return fileDict

    def projectWithName(self, projectName, projectPath=None):
        """
        Get or create project by name.

        Returns project UUID.
        """
        logger.info(f"Getting/creating project: {projectName}")

        # Try to get existing project
        try:
            project = models.Project.objects.get(name=projectName)
            logger.info(f"Found existing project: {project.uuid}")
            return project.uuid
        except models.Project.DoesNotExist:
            pass

        # Create new project
        if projectPath is None:
            projectPath = settings.CCP4I2_PROJECTS_DIR / slugify(projectName)

        # Use serializer for validation
        serializer = serializers.ProjectSerializer(
            data={
                "name": projectName,
                "directory": str(Path(projectPath).resolve()),
            }
        )

        if not serializer.is_valid():
            raise ValueError(f"Invalid project data: {serializer.errors}")

        # Save creates and saves in one step
        newProject = serializer.save()

        logger.info(
            f'Created project "{newProject.name}" at '
            f'{Path(newProject.directory).resolve()} with uuid {newProject.uuid}'
        )
        return newProject.uuid

    def projectJobWithTask(self, projectId, task_name=None):
        """
        Create a new job for the given task in the project.

        Returns job UUID (without hyphens).
        """
        logger.info(f"Creating job for task '{task_name}' in project {projectId}")

        created_job_uuid = create_job(projectId=str(projectId), taskName=task_name)

        logger.info(f"Created job {created_job_uuid} for task '{task_name}'")
        return created_job_uuid.replace("-", "")

    def pluginWithArgs(self, parsed_args, workDirectory=None, jobId=None):
        """
        Create plugin instance and populate with command-line arguments.

        Uses PluginPopulator component for clean separation of concerns.
        Attaches Django database handler for database-aware file operations.
        """
        # Get job and work directory in one query (with select_related for project)
        job = None
        if jobId is not None:
            job = models.Job.objects.select_related('project').get(uuid=jobId)
            if workDirectory is None:
                workDirectory = job.directory

        logger.info(f"Work directory: {workDirectory}")

        # Create plugin instance
        from ccp4i2.core.CCP4TaskManager import TASKMANAGER
        thePlugin = TASKMANAGER().get_plugin_class(
            parsed_args.task_name
        )(jobId=jobId, workDirectory=workDirectory, parent=None)

        self.work_directory = workDirectory

        # Attach Django database handler for database-aware operations
        if job is not None:
            from ..db.async_db_handler import AsyncDatabaseHandler
            thePlugin._dbHandler = AsyncDatabaseHandler(project_uuid=str(job.project.uuid))
            thePlugin._dbProjectId = str(job.project.uuid).replace("-", "")
            logger.debug(f"Attached dbHandler with project {job.project.uuid}")

        # Populate plugin using PluginPopulator component
        allKeywords = self.keywordsOfTaskName(parsed_args.task_name)
        PluginPopulator.populate(thePlugin, parsed_args, allKeywords)

        # DEBUG: Check UNMERGEDFILES state after population
        if hasattr(thePlugin.container, 'inputData') and hasattr(thePlugin.container.inputData, 'UNMERGEDFILES'):
            uf = thePlugin.container.inputData.UNMERGEDFILES
            print(f"\n[DEBUG pluginWithArgs] After population:")
            print(f"[DEBUG pluginWithArgs]   UNMERGEDFILES type: {type(uf).__name__}")
            print(f"[DEBUG pluginWithArgs]   UNMERGEDFILES length: {len(uf)}")
            print(f"[DEBUG pluginWithArgs]   UNMERGEDFILES isSet(): {uf.isSet()}")
            if len(uf) > 0:
                print(f"[DEBUG pluginWithArgs]   First item type: {type(uf[0]).__name__}")
                print(f"[DEBUG pluginWithArgs]   First item isSet(): {uf[0].isSet()}")
                if hasattr(uf[0], 'file'):
                    print(f"[DEBUG pluginWithArgs]   First item.file: {uf[0].file}")
                    print(f"[DEBUG pluginWithArgs]   First item.file.isSet(): {uf[0].file.isSet()}")
                    print(f"[DEBUG pluginWithArgs]   First item.file.baseName: {uf[0].file.baseName}")
            print()

        # Save params to database (creates input_params.xml)
        if job is not None:
            print(f"[DEBUG pluginWithArgs] About to call save_params_for_job, exclude_unset=True")
            save_params_for_job(thePlugin, the_job=job, exclude_unset=True)
            print(f"[DEBUG pluginWithArgs] save_params_for_job returned")

        return thePlugin

    def execute(self):
        """
        Execute the job after validation.

        Returns (jobId, exit_code) tuple.
        """
        thePlugin = self.getPlugin(arguments_parsed=True)

        if self.jobId is None:
            raise ValueError("jobId not set - cannot execute")
        if self.projectId is None:
            raise ValueError("projectId not set - cannot execute")

        # Save params.xml (job output metadata)
        #thePlugin.saveParams()
        print(self.jobId)

        from ccp4i2.core.base_object.error_reporting import Severity

        # Step 1: Run validity() - allows plugins to adjust qualifiers for embedded wrappers
        # (e.g., servalcat_pipe sets allowUndefined on metalCoordWrapper.inputData.XYZIN)
        # Note: validity() errors are informational only for now - we don't block execution
        # This allows pipeline developers time to implement proper validity checks
        validity_error = thePlugin.validity()

        print("\n" + "=" * 80)
        print("VALIDITY CHECK")
        print("=" * 80)
        if validity_error and validity_error.maxSeverity() >= Severity.WARNING:
            print(validity_error.report(severity_threshold=Severity.WARNING))
            if validity_error.maxSeverity() >= Severity.ERROR:
                print("\n⚠️  Validity errors found (continuing anyway)")
                logger.warning(
                    "Job %s validity check has errors (non-blocking): %s",
                    self.jobId,
                    validity_error.report(severity_threshold=Severity.ERROR)
                )
        else:
            print("✓ No validity issues")
        print("=" * 80)

        # Step 2: Run checkInputData() - validates input file paths (this IS a hard failure)
        input_error = thePlugin.checkInputData()

        print("\n" + "=" * 80)
        print("INPUT DATA CHECK")
        print("=" * 80)
        if input_error and input_error.maxSeverity() >= Severity.WARNING:
            print(input_error.report(severity_threshold=Severity.WARNING))
        else:
            print("✓ All input files valid")
        print("=" * 80)

        if input_error and input_error.maxSeverity() >= Severity.ERROR:
            print("\n❌ INPUT DATA CHECK FAILED - Cannot execute job\n")
            logger.error(
                "Job %s validation failed: %s",
                self.jobId,
                input_error.report(severity_threshold=Severity.ERROR)
            )
            return self.jobId, 1  # Failure

        # Save input parameters to input_params.xml so async_run_job can load them
        # This is essential because async_run_job creates a fresh plugin instance
        # and needs to populate it with the command-line arguments we just processed
        from ..db import models
        job = models.Job.objects.get(uuid=self.jobId)
        save_params_for_job(thePlugin, the_job=job, mode="JOB_INPUT", exclude_unset=True)
        logger.info(f"Saved input parameters to input_params.xml for async runner to load")

        # Execute job using async runner
        from asgiref.sync import async_to_sync
        from ccp4i2.lib.async_run_job import run_job_async

        try:
            async_to_sync(run_job_async)(self.jobId)
            logger.info("Job %s executed successfully", self.jobId)
            return self.jobId, 0  # Success
        except Exception as e:
            logger.error("Job %s execution failed: %s", self.jobId, str(e))
            return self.jobId, 1  # Failure
