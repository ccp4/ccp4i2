"""
CCP4X Job Management API ViewSet

This module provides a comprehensive REST API for managing CCP4 crystallographic computing jobs
within the CCP4i2 Django application. It handles job lifecycle operations including creation,
execution, monitoring, cloning, and export functionality.

The JobViewSet follows Django REST Framework patterns and integrates with the CCP4 task management
system to provide a robust interface for computational crystallography workflows.

Classes:
    JobViewSet: Main REST API viewset for job operations

Dependencies:
    - Django REST Framework for API structure
    - CCP4i2 core libraries for crystallographic computing
    - XML processing for job parameters and reports
    - Subprocess management for job execution

Author: CCP4i2 Development Team
License: CCP4 License
Version: Compatible with CCP4i2 and Django 4.2+

    - Supports distributed file storage
    - Handles container-based job execution
    - Includes proper error handling for cloud environments
"""

# Add these imports at the top of the file (after existing imports)
import logging
import datetime
import json
import os
from xml.etree import ElementTree as ET
from pytz import timezone
from django.http import Http404
from rest_framework.response import Response
from rest_framework.parsers import MultiPartParser, JSONParser
from ccp4i2.core import CCP4TaskManager
from ccp4i2.core.CCP4Container import CContainer
from ccp4i2.core import CCP4ErrorHandling
from ..lib.utils.jobs.i2run import i2run_for_job
from ..lib.utils.parameters.load_xml import load_nested_xml
# validate_container no longer used - validation/ endpoint now uses unified validate_job utility
from ..lib.utils.files.digest import digest_param_file
from ..lib.utils.containers.validate import getEtree  # Still used for error handling in other endpoints
from ..lib.utils.parameters.set_input_by_context import set_input_by_context_job
from ..lib.utils.jobs.preview import preview_job
import tempfile
from pathlib import Path
from django.http import FileResponse
from ..db.export_project import export_project_to_zip
import subprocess
from xml.etree import ElementTree as ET
from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.viewsets import ModelViewSet
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework import status
from rest_framework.permissions import IsAuthenticated
from . import serializers
from ..db import models
from ..lib.utils.reporting.i2_report import generate_job_report
from ..lib.utils.jobs.clone import clone_job  # Modern clone utility with Result pattern
from ..lib.utils.navigation.dependencies import find_dependent_jobs
from ..lib.utils.navigation.dependencies import delete_job_and_dependents
# Modern utilities
from ..lib.utils.plugins.get_plugin import get_job_plugin
from ..lib.utils.containers.json_encoder import CCP4i2JsonEncoder

# Modern imports - all now using modern utilities
from ..lib.utils.files.upload_param import upload_file_param
from ..lib.utils.containers.get_container import get_job_container
from ..lib.utils.helpers.object_method import object_method
from ..lib.utils.navigation.what_next import get_what_next
from django.http import JsonResponse
from django.conf import settings
from django.utils.text import slugify

# Uniform API response helpers
from ..lib.response import api_success, api_error

logger = logging.getLogger(f"ccp4i2:{__name__}")


class JobViewSet(ModelViewSet):
    """
    Django REST Framework ViewSet for managing CCP4 crystallographic computing jobs.

    This ViewSet provides comprehensive job management functionality including:
    - Standard CRUD operations (list, create, retrieve, update, delete)
    - Job execution and monitoring
    - Parameter management and validation
    - File handling and export
    - Job cloning and dependency tracking
    - XML report generation


    Attributes:
        queryset (QuerySet): All Job model instances
        serializer_class (JobSerializer): Serializer for Job model
        parser_classes (list): Supported content parsers for requests
        filterset_fields (list): Fields available for filtering

    API Endpoints:
        Standard REST endpoints:
        - GET /api/jobs/ - List all jobs
        - POST /api/jobs/ - Create new job
        - GET /api/jobs/{id}/ - Retrieve specific job
        - PUT /api/jobs/{id}/ - Update job
        - DELETE /api/jobs/{id}/ - Delete job and dependencies

        Custom action endpoints:
        - GET /api/jobs/{id}/what_next/ - Get suggested next steps
        - POST /api/jobs/{id}/set_context_job/ - Set context for job input
        - POST /api/jobs/{id}/object_method/ - Execute method on job object
        - GET /api/jobs/{id}/params_xml/ - Get job parameters as XML
        - GET /api/jobs/{id}/report_xml/ - Get job report as XML
        - GET /api/jobs/{id}/dependent_jobs/ - Get jobs that depend on this job
        - POST /api/jobs/{id}/clone/ - Clone an existing job
        - POST /api/jobs/{id}/run/ - Execute a job
        - GET /api/jobs/{id}/container/ - Get job container data
        - GET /api/jobs/{id}/diagnostic_xml/ - Get diagnostic information
        - GET /api/jobs/{id}/digest/ - Get file digest information
        - GET /api/jobs/{id}/i2run_command/ - Get command line for job execution
        - GET /api/jobs/{id}/digest_param_file/ - Digest specific parameter file
        - GET /api/jobs/{id}/def_xml/ - Get job definition XML
        - GET /api/jobs/{id}/validation/ - Validate job parameters
        - POST /api/jobs/{id}/set_parameter/ - Set job parameter
        - POST /api/jobs/{id}/upload_file_param/ - Upload file parameter
        - POST /api/jobs/{id}/preview/ - Preview job with external viewer
        - GET /api/jobs/{id}/files/ - Get files associated with job
        - GET /api/jobs/{id}/export_job/ - Export job as ZIP archive

        - Jobs execute in container environments
        - Subprocess execution adapted for container constraints
        - Error handling includes cloud-specific scenarios
        - Export functionality supports blob storage integration

    Security Notes:
        - All endpoints require appropriate authentication
        - File operations validate paths to prevent directory traversal
        - Subprocess execution uses secure parameter validation
        - Sensitive data is masked in logs

    Example Usage:
        ```python
        # Get job report
        GET /api/jobs/123/report_xml/

        # Clone a job
        POST /api/jobs/123/clone/

        # Run a job
        POST /api/jobs/123/run/

        # Export job data
        GET /api/jobs/123/export_job/
        ```
    """

    queryset = models.Job.objects.all()
    serializer_class = serializers.JobSerializer
    parser_classes = [FormParser, MultiPartParser, JSONParser]
    filterset_fields = ["project"]
    permission_classes = [IsAuthenticated]

    def destroy(self, request, *args, **kwargs):
        """
        Delete a job and all its dependent jobs.

        This method safely removes a job and all jobs that depend on it,
        maintaining data integrity across the project hierarchy.

        Args:
            request (Request): HTTP request object
            *args: Variable length argument list
            **kwargs: Arbitrary keyword arguments

        Returns:
            Response: Success response with status message

        Raises:
            Http404: If the job is not found

        Example:
            DELETE /api/jobs/123/

            - Ensures proper cleanup of container-based job artifacts
        """
        try:
            instance = self.get_object()
            logger.warning("Deleting job %s", instance)
            delete_job_and_dependents(instance)
            # Note: Adding response body to prevent JavaScript network error
            # when response has no body content
            return api_success({"deleted": True})
        except Http404:
            return api_error("Job not found", status=404)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def what_next(self, request, pk=None):
        """
        Get suggested next steps for workflow continuation.

        Analyzes the current job state and provides recommendations for
        subsequent analysis steps based on CCP4 crystallographic workflows.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            JsonResponse: Dictionary containing suggested next steps

        Response Format:
            {
                "status": "Success",
                "suggestions": [
                    {
                        "task_name": "refmac5",
                        "description": "Refinement with Refmac5",
                        "priority": "high"
                    }
                ]
            }

        Example:
            GET /api/jobs/123/what_next/
        """
        try:
            job = models.Job.objects.get(id=pk)
            what_next_data = get_what_next(job)
            return api_success(what_next_data)
        except ValueError as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=400)
        except models.Job.DoesNotExist as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)
        except Exception as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def set_context_job(self, request, pk=None):
        """
        Set context job for automatic input parameter configuration.

        Links the current job to a context job, automatically configuring
        input parameters based on the context job's output files.

        Args:
            request (Request): HTTP request with JSON body containing:
                - context_job_uuid (str): UUID of the context job
            pk (int): Primary key of the target job

        Returns:
            JsonResponse: Updated job data or error message

        Request Body:
            {
                "context_job_uuid": "550e8400-e29b-41d4-a716-446655440000"
            }

        Response Format:
            {
                "status": "Success",
                "new_job": {<serialized job data>}
            }

        Example:
            POST /api/jobs/123/set_context_job/
            {
                "context_job_uuid": "550e8400-e29b-41d4-a716-446655440000"
            }
        """
        try:
            job = models.Job.objects.get(id=pk)
            form_data = json.loads(request.body.decode("utf-8"))
            context_job_uuid = form_data["context_job_uuid"]
            set_input_by_context_job(str(job.uuid), context_job_uuid)
            serializer = serializers.JobSerializer(job)
            return api_success({"new_job": serializer.data})
        except CCP4ErrorHandling.CException as err:
            error_tree = getEtree(err)
            ET.indent(error_tree, " ")
            return api_error(ET.tostring(error_tree).decode("utf-8"), status=400)
        except models.Job.DoesNotExist as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)
        except Exception as err:
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def object_method(self, request, pk=None):
        """
        Execute a method on a job object within the CCP4 framework.

        Provides dynamic access to CCP4 job object methods for advanced
        parameter manipulation and data extraction.

        Args:
            request (Request): HTTP request with JSON body containing:
                - object_path (str): Path to the target object
                - method_name (str): Name of method to execute
                - args (list, optional): Method arguments
                - kwargs (dict, optional): Method keyword arguments
            pk (int): Primary key of the job

        Returns:
            JsonResponse: Method execution result or error

        Request Body:
            {
                "object_path": "inputData.XYZIN",
                "method_name": "getFileName",
                "args": [],
                "kwargs": {}
            }

        Response Format:
            {
                "status": "Success",
                "result": <method return value>
            }

        Example:
            POST /api/jobs/123/object_method/

        Security Notes:
            - Method execution is sandboxed within CCP4 framework
            - Only CCP4-defined methods are accessible
        """
        form_data = json.loads(request.body.decode("utf-8"))
        job = models.Job.objects.get(id=pk)
        object_path = form_data["object_path"]
        method_name = form_data["method_name"]
        args = form_data.get("args", [])
        kwargs = form_data.get("kwargs", {})
        try:
            result = object_method(job, object_path, method_name, args, kwargs)
            logger.debug("result %s", result)
            return api_success({"result": result})
        except CCP4ErrorHandling.CException as err:
            error_tree = getEtree(err)
            logger.debug("error_tree %s", error_tree)
            ET.indent(error_tree, " ")
            return api_error(ET.tostring(error_tree).decode("utf-8"), status=400)

    @action(
        detail=True,
        methods=["get", "put"],
        serializer_class=serializers.JobSerializer,
    )
    def params_xml(self, request, pk=None):
        """
        Retrieve or update job parameters as XML document.

        GET: Returns the job's parameter configuration in XML format, either from
        params.xml (for completed jobs) or input_params.xml (for pending jobs).

        PUT: Updates the job's input_params.xml with the provided XML content.
        Only allowed for jobs in PENDING status.

        Uses the unified utility from ccp4i2.lib.utils.jobs.reports for
        consistent behavior with CLI commands.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: XML content or error message

        Response Format (GET):
            {
                "status": "Success",
                "xml": "<xml>...</xml>"
            }

        Response Format (PUT):
            {
                "status": "Success",
                "message": "Parameters saved successfully"
            }

        File Priority (GET):
            1. params.xml (for running/completed jobs)
            2. input_params.xml (fallback for pending jobs)

        Example:
            GET /api/jobs/123/params_xml/
            PUT /api/jobs/123/params_xml/  (body: {"xml": "<CCP4i2_body>...</CCP4i2_body>"})

        Architecture:
            - Uses get_job_params_xml() utility for GET
            - Direct file write for PUT (validation via CPluginScript planned)
            - Shared with get_job_report --type params command
            - Consistent file path logic and error handling
        """
        try:
            the_job = models.Job.objects.get(id=pk)

            if request.method == "PUT":
                return self._put_params_xml(request, the_job)

            # GET request - use unified utility
            from ..lib.utils.jobs.reports import get_job_params_xml

            result = get_job_params_xml(the_job)

            if result.success:
                return api_success({"xml": result.data})
            else:
                return api_error(result.error, status=404)

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(f"Job not found: {str(err)}", status=404)
        except Exception as err:
            logger.exception("Unexpected error getting params XML for job %s", pk, exc_info=err)
            return api_error(f"Unexpected error: {str(err)}", status=500)

    def _put_params_xml(self, request, job):
        """
        Handle PUT request to update job's input_params.xml.

        Only allowed for jobs in PENDING status. Validates the XML is well-formed
        before saving.

        Args:
            request: HTTP request with JSON body containing 'xml' key
            job: Job model instance

        Returns:
            Response with success/error status
        """
        import xml.etree.ElementTree as ET

        # Check job status - only allow editing pending jobs
        if job.status not in [models.Job.Status.UNKNOWN, models.Job.Status.PENDING]:
            return api_error(
                f"Cannot modify parameters on job with status '{job.get_status_display()}'. "
                f"Only PENDING jobs can be edited.",
                status=400,
                details={
                    "job_id": str(job.uuid),
                    "job_status": job.status,
                    "allowed_statuses": ["UNKNOWN", "PENDING"]
                }
            )

        try:
            # Parse request body
            body = json.loads(request.body.decode("utf-8"))
            xml_content = body.get("xml")

            if not xml_content:
                return api_error("Missing 'xml' field in request body", status=400)

            # Validate XML is well-formed
            try:
                ET.fromstring(xml_content)
            except ET.ParseError as e:
                return api_error(
                    f"Invalid XML: {str(e)}",
                    status=400,
                    details={"parse_error": str(e)}
                )

            # Write to input_params.xml
            input_params_file = job.directory / "input_params.xml"
            logger.info("Saving params XML to %s for job %s", input_params_file, job.uuid)

            with open(input_params_file, "w", encoding="utf-8") as f:
                f.write(xml_content)

            logger.info("Successfully saved params XML to %s", input_params_file)
            return api_success({
                "message": "Parameters saved successfully",
                "file": str(input_params_file)
            })

        except json.JSONDecodeError as e:
            return api_error(f"Invalid JSON in request body: {str(e)}", status=400)
        except IOError as e:
            logger.exception("Failed to write params XML for job %s", job.uuid)
            return api_error(f"Failed to write file: {str(e)}", status=500)
        except Exception as e:
            logger.exception("Unexpected error saving params XML for job %s", job.uuid)
            return api_error(f"Unexpected error: {str(e)}", status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def report_xml(self, request, pk=None):
        """
        Generate and retrieve XML report for a job.

        Creates a comprehensive XML report containing job results, statistics,
        and analysis outcomes. Reports are cached for performance and regenerated
        as needed based on job status.

        Uses the unified utility from ccp4i2.lib.utils.jobs.reports for
        consistent behavior with CLI commands.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: XML report content with cache headers

        Response Format:
            {
                "status": "Success",
                "xml": b"<report>...</report>"
            }

        Report Content:
            - Job execution statistics
            - Result quality metrics
            - Output file summaries
            - Error/warning messages
            - Crystallographic statistics (if applicable)

        Caching Strategy:
            - Generated reports cached to report_xml.xml
            - Cache invalidated on job status changes
            - No-cache headers for real-time updates

        Example:
            GET /api/jobs/123/report_xml/

        Architecture:
            - Uses get_job_report_xml() utility
            - Shared with get_job_report --type report command
            - Consistent caching logic and report generation
        """
        try:
            the_job = models.Job.objects.get(id=pk)

            # Use unified utility
            from ..lib.utils.jobs.reports import get_job_report_xml

            # Check if regeneration requested
            regenerate = request.GET.get('regenerate', 'false').lower() == 'true'

            result = get_job_report_xml(the_job, regenerate=regenerate)

            if result.success:
                response = api_success({"xml": result.data})
                # Add no-cache headers for real-time updates
                response["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
                response["Pragma"] = "no-cache"
                response["Expires"] = "0"
                return response
            else:
                error_status = 404 if "not found" in result.error.lower() else 500
                return api_error(result.error, status=error_status)

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(f"Job not found: {str(err)}", status=404)
        except Exception as err:
            logger.exception("Unexpected error getting report XML for job %s", pk, exc_info=err)
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def regenerate_report(self, request, pk=None):
        """
        Force regeneration of the job's report XML.

        Deletes any cached report and regenerates it from scratch. This is useful
        when report generation logic has been updated or when the cached report
        may be stale or corrupted.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: Success status with regenerated report XML

        Response Format:
            {
                "status": "Success",
                "data": {
                    "regenerated": true,
                    "xml": "<report>...</report>"
                }
            }

        Example:
            POST /api/jobs/123/regenerate_report/
        """
        try:
            the_job = models.Job.objects.get(id=pk)

            # Use unified utility with regenerate=True
            from ..lib.utils.jobs.reports import get_job_report_xml

            result = get_job_report_xml(the_job, regenerate=True)

            if result.success:
                response = api_success({
                    "regenerated": True,
                    "xml": result.data
                })
                # Add no-cache headers
                response["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
                response["Pragma"] = "no-cache"
                response["Expires"] = "0"
                return response
            else:
                error_status = 404 if "not found" in result.error.lower() else 500
                return api_error(result.error, status=error_status)

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(f"Job not found: {str(err)}", status=404)
        except Exception as err:
            logger.exception("Unexpected error regenerating report for job %s", pk, exc_info=err)
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def dependent_jobs(self, request, pk=None):
        """
        Retrieve jobs that depend on the specified job.

        Returns a list of jobs that use outputs from the current job as inputs,
        enabling dependency tracking and impact analysis for job modifications.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: Serialized list of dependent jobs

        Response Format:
            [
                {
                    "id": 456,
                    "title": "Refinement Job",
                    "task_name": "refmac5",
                    "status": "FINISHED",
                    ...
                }
            ]

        Use Cases:
            - Impact analysis before job deletion
            - Workflow visualization
            - Dependency validation

        Example:
            GET /api/jobs/123/dependent_jobs/
        """
        try:
            the_job = models.Job.objects.get(id=pk)
            dependent_jobs = find_dependent_jobs(the_job)
            serializer = serializers.JobSerializer(dependent_jobs, many=True)
            # DRF standard for list endpoints - return array directly
            return Response(serializer.data)
        except (ValueError, models.Job.DoesNotExist) as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def clone(self, request, pk=None):
        """
        Create a copy of an existing job with identical parameters.

        Clones a job including all input parameters, allowing users to rerun
        analyses with modifications or create variations of successful workflows.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job to clone

        Returns:
            Response: Serialized data of the newly created job

        Response Format:
            {
                "id": 789,
                "uuid": "new-uuid-here",
                "title": "Copy of Original Job",
                "task_name": "same_task",
                "status": "PENDING",
                ...
            }

        Cloning Process:
            1. Copy job configuration and parameters
            2. Generate new UUID and job number
            3. Reset status to PENDING
            4. Preserve input file references
            5. Clear output data

        Example:
            POST /api/jobs/123/clone/

            - Preserves file references without duplication
        """
        try:
            old_job_id = models.Job.objects.get(id=pk).uuid
            result = clone_job(old_job_id)

            if result.success:
                serializer = serializers.JobSerializer(result.data)
                # DRF standard for create - return serializer data directly
                return Response(serializer.data, status=201)
            else:
                return result.to_response(error_status=400)

        except models.Job.DoesNotExist as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def run(self, request, pk=None):
        """
        Execute a job using environment-appropriate backend.

        Automatically adapts to deployment context:
        - Local Mode: Executes job via subprocess (laptop/development)
        - Azure Mode: Queues job via Service Bus (container apps)

        The execution mode is determined from environment variables.
        See ccp4i2.lib.context_dependent_run for implementation details.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job to execute

        Returns:
            Response: Updated job data with appropriate status

        Example:
            POST /api/jobs/123/run/

        Environment Variables:
            EXECUTION_MODE: Explicit mode ('local' or 'azure')
            SERVICE_BUS_CONNECTION_STRING: Azure connection (implies azure)
            CCP4: Path to CCP4 installation (for local mode)
        """
        try:
            from ..lib.utils.jobs.context_run import run_job_context_aware

            job = models.Job.objects.get(id=pk)

            # Execute job using context-aware backend
            result = run_job_context_aware(job)

            if result["success"]:
                serializer = serializers.JobSerializer(result["data"])
                # DRF standard - return serializer data directly
                return Response(serializer.data)
            else:
                return api_error(result["error"], status=result["status"])

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)
        except Exception as err:
            logger.exception("Unexpected error running job %s", pk, exc_info=err)
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def run_local(self, request, pk=None):
        """
        Execute a job locally regardless of environment configuration.

        Forces local execution even in Azure environments, useful for:
        - Tasks requiring direct filesystem access
        - Interactive or GUI-based tasks
        - Tasks with specific local dependencies
        - Development and debugging scenarios

        This endpoint bypasses the automatic environment detection and
        always executes jobs via subprocess on the current machine.

        Args:
            request (Request): HTTP request object with optional JSON body:
                - synchronous (bool): If True, blocks until job completes.
                    Default is False (returns immediately after starting).
            pk (int): Primary key of the job to execute

        Returns:
            Response: Updated job data with appropriate status.
            When synchronous=True, returns the final job state after completion.

        Request Body (optional):
            {
                "synchronous": true
            }

        Example:
            # Asynchronous (default) - returns immediately
            POST /api/jobs/123/run_local/

            # Synchronous - blocks until job completes
            POST /api/jobs/123/run_local/
            {"synchronous": true}

        Note:
            - Requires CCP4 installation on local machine
            - May not work in pure container environments
            - Frontend can selectively use this for specific task types
            - Synchronous mode is useful for scripting and automation
        """
        try:
            from ..lib.utils.jobs.context_run import run_job_context_aware

            job = models.Job.objects.get(id=pk)

            # Parse synchronous parameter from request body
            synchronous = False
            if request.body:
                try:
                    body_data = json.loads(request.body.decode("utf-8"))
                    synchronous = body_data.get("synchronous", False)
                except (json.JSONDecodeError, UnicodeDecodeError):
                    pass  # Use default (False) if body is invalid

            logger.info(
                "Forcing local execution for job %s (uuid=%s, task=%s, synchronous=%s)",
                job.id,
                job.uuid,
                job.task_name,
                synchronous,
            )

            # Execute job with forced local mode
            result = run_job_context_aware(job, force_local=True, synchronous=synchronous)

            if result["success"]:
                serializer = serializers.JobSerializer(result["data"])
                # DRF standard - return serializer data directly
                return Response(serializer.data)
            else:
                return api_error(result["error"], status=result["status"])

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)
        except Exception as err:
            logger.exception(
                "Unexpected error running job locally %s", pk, exc_info=err
            )
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def container(self, request, pk=None):
        """
        Retrieve job container data as JSON.

        Returns the CCP4 container object data in JSON format, providing
        access to the complete job configuration and runtime information.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: JSON representation of job container

        Response Format:
            {
                "status": "Success",
                "result": {
                    "inputData": {...},
                    "outputData": {...},
                    "controlParameters": {...}
                }
            }

        Container Data:
            - Input data definitions and values
            - Output data specifications
            - Control parameters and settings
            - Task-specific configuration

        Example:
            GET /api/jobs/123/container/
        """
        try:
            the_job = models.Job.objects.get(id=pk)

            # Modern approach: Use CPluginScript architecture
            # get_job_plugin automatically creates a dbHandler for file path resolution
            plugin = get_job_plugin(the_job)

            # Serialize container to JSON using modern encoder
            container_json = json.dumps(
                plugin.container,
                cls=CCP4i2JsonEncoder
            )

            return api_success({"result": json.loads(container_json)})

        except models.Job.DoesNotExist as err:
            logger.exception("Job %s not found", pk, exc_info=err)
            return api_error("Job not found", status=404)
        except FileNotFoundError as err:
            logger.exception(
                "Failed to find parameters for job %s",
                pk,
                exc_info=err,
            )
            return api_error("Job parameters not found", status=404)
        except Exception as err:
            logger.exception(
                "Failed to get container for job %s",
                pk,
                exc_info=err,
            )
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def diagnostic_xml(self, request, pk=None):
        """
        Retrieve diagnostic information as XML.

        Returns detailed diagnostic data generated during job execution,
        including error messages, warnings, and debugging information.

        Uses the unified utility from ccp4i2.lib.utils.jobs.reports for
        consistent behavior with CLI commands.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: XML diagnostic content or error message

        Response Format:
            {
                "status": "Success",
                "xml": "<diagnostic>...</diagnostic>"
            }

        Diagnostic Content:
            - Execution errors and warnings
            - Resource usage statistics
            - Input validation results
            - Program-specific debug information

        Example:
            GET /api/jobs/123/diagnostic_xml/

        Note:
            Diagnostic files are created during job execution and may not
            exist for jobs that haven't run or failed early.

        Architecture:
            - Uses get_job_diagnostic_xml() utility
            - Shared with get_job_report --type diagnostic command
            - Consistent error handling
        """
        try:
            the_job = models.Job.objects.get(id=pk)

            # Use unified utility
            from ..lib.utils.jobs.reports import get_job_diagnostic_xml

            result = get_job_diagnostic_xml(the_job)

            if result.success:
                return api_success({"xml": result.data})
            else:
                return api_error(result.error, status=404)

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(f"Job not found: {str(err)}", status=404)
        except Exception as err:
            logger.exception("Unexpected error getting diagnostic XML for job %s", pk, exc_info=err)
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def digest(self, request, pk=None):
        """
        Generate digest information for a specific object path.

        Creates a summary digest of data at the specified object path,
        providing quick access to key information without full data retrieval.

        Args:
            request (Request): HTTP request with query parameters:
                - object_path (str): Path to the object to digest
            pk (int): Primary key of the job

        Returns:
            Response: Digest information dictionary

        Query Parameters:
            - object_path: CCP4 object path (e.g., "inputData.XYZIN")

        Example:
            GET /api/jobs/123/digest/?object_path=inputData.XYZIN
        """
        try:
            the_job = models.Job.objects.get(id=pk)
            object_path = request.GET.get("object_path", "")
            # Strip trailing slash if present (but don't strip last character otherwise)
            if object_path.endswith("/"):
                object_path = object_path[:-1]
            logger.info("Digesting file %s", object_path)
            response_dict = digest_param_file(the_job, object_path)
            return api_success(response_dict)
        except (ValueError, models.Job.DoesNotExist) as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=400)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def i2run_command(self, request, pk=None):
        """
        Generate command line equivalent for job execution.

        Returns the i2run command that would execute this job from the
        command line, useful for debugging and external execution.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: Command string or error message

        Response Format:
            {
                "status": "Success",
                "command": "i2run task_name -p project_path ..."
            }

        Example:
            GET /api/jobs/123/i2run_command/
        """
        try:
            the_job = models.Job.objects.get(id=pk)
            response_string = i2run_for_job(the_job)
            return api_success({"command": response_string})
        except (ValueError, models.Job.DoesNotExist) as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=400)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def digest_param_file(self, request, pk=None):
        """
        Generate digest for a specific parameter file.

        Creates a digest summary of a file parameter, automatically determining
        whether it's an input or output file based on job associations.

        Args:
            request (Request): HTTP request with query parameters:
                - job_param_name (str): Name of the job parameter
            pk (int): Primary key of the job

        Returns:
            Response: File digest information

        Response Format:
            {
                "status": "Success",
                "digest": {
                    "file_type": "PDB",
                    "content_summary": "...",
                    "statistics": {...}
                }
            }

        Example:
            GET /api/jobs/123/digest_param_file/?job_param_name=XYZIN
        """
        job_param_name = request.GET.get("job_param_name")
        object_path = job_param_name[:-1]
        try:
            the_job = models.Job.objects.get(id=pk)
            the_file = models.File.objects.get(
                job=the_job, job_param_name=job_param_name[:-1]
            )
            imports = models.FileImport.objects.filter(file=the_file)
            is_import = imports.count() > 0
            if is_import:
                object_path = f"{the_job.task_name}.inputData.{job_param_name[:-1]}"
            else:
                object_path = f"{the_job.task_name}.outputData.{job_param_name[:-1]}"
            response_dict = digest_param_file(the_job, object_path)
            return api_success({"digest": response_dict})
        except (ValueError, models.Job.DoesNotExist) as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=400)
        except Exception as err:
            logging.exception(
                "Failed to digest file %s %s", pk, object_path, exc_info=err
            )
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def def_xml(self, request, pk=None):
        """
        Retrieve task definition XML for the job.

        Returns the CCP4 task definition XML that defines the structure
        and parameters for the job's computational task.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: Unpacked XML task definition

        Response Format:
            {
                "status": "Success",
                "xml": b"<task_definition>...</task_definition>"
            }

        XML Content:
            - Task parameter definitions
            - Input/output specifications
            - Validation rules
            - GUI layout information

        Example:
            GET /api/jobs/123/def_xml/
        """
        try:
            the_job = models.Job.objects.get(id=pk)
            def_xml_path = CCP4TaskManager.TASKMANAGER().locate_def_xml(
                task_name=the_job.task_name
            )
            with open(def_xml_path, "r") as def_xml_file:
                def_xml = def_xml_file.read()
                packedXML = ET.fromstring(def_xml)
                unpackedXML = load_nested_xml(packedXML)
                ET.indent(unpackedXML, " ")
                return api_success({"xml": ET.tostring(unpackedXML).decode("utf-8")})
        except (ValueError, models.Job.DoesNotExist) as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=400)
        except FileNotFoundError as err:
            logger.exception(
                "Failed to find file %s",
                def_xml_path,
                exc_info=err,
            )
            return api_error(str(err), status=404)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def validation(self, request, pk=None):
        """
        Validate job parameters and return error report.

        Performs comprehensive validation of job parameters against
        task requirements and returns detailed error information in XML format.

        Uses the unified CPluginScript architecture for proper container
        hierarchy and consistent behavior with CLI validation commands.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: XML validation error report

        Response Format:
            {
                "status": "Success",
                "xml": b"<error_report>...</error_report>"
            }

        Validation Checks:
            - Required parameter validation
            - Data type verification
            - File existence checks
            - Parameter dependency validation
            - Value range verification

        Error Report Structure:
            - Parameter-specific errors
            - Severity levels (error, warning, info)
            - Suggested corrections
            - Dependency conflict details

        Example:
            GET /api/jobs/123/validation/

        Architecture:
            - Uses CPluginScript + dbHandler for proper context
            - Shared with validate_job management command
            - Uses new CErrorReport.getErrors() API
            - Consistent Result[T] pattern for error handling
        """
        try:
            the_job = models.Job.objects.get(id=pk)

            # Use unified utility (CPluginScript architecture)
            from ..lib.utils.jobs.validate import validate_job

            result = validate_job(the_job)

            if result.success:
                error_etree = result.data
                # Remove stack traces from output (too verbose for API)
                stack_elements = error_etree.findall(".//stack")
                for stack_element in stack_elements:
                    parent = error_etree.find(f".//*[stack='{stack_element}']/..")
                    if parent is not None:
                        parent.remove(stack_element)

                ET.indent(error_etree, " ")
                return api_success({"xml": ET.tostring(error_etree).decode("utf-8")})
            else:
                return api_error(result.error, status=400, details=result.error_details)

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(f"Job not found: {str(err)}", status=404)
        except Exception as err:
            logger.exception(
                "Unexpected error validating job %s", pk, exc_info=err
            )
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def set_parameter(self, request, pk=None):
        """
        Set a parameter value in the job's input configuration.

        Updates a specific parameter in the job's input_params.xml file,
        allowing dynamic modification of job configuration before execution.

        Uses the unified CPluginScript architecture for proper database
        synchronization and consistent behavior with CLI commands.

        Args:
            request (Request): HTTP request with JSON body containing:
                - object_path (str): Path to the parameter to modify
                - value (any): New value for the parameter
            pk (int): Primary key of the job

        Returns:
            JsonResponse: Status and updated parameter information

        Request Body:
            {
                "object_path": "inputData.XYZIN.fileName",
                "value": "/path/to/new/file.pdb"
            }

        Response Format:
            {
                "status": "Success",
                "data": {
                    "object_path": "inputData.XYZIN.fileName",
                    "value_set": true,
                    "message": "Parameter set successfully"
                }
            }

        Example:
            POST /api/jobs/123/set_parameter/

        Security Notes:
            - Parameter paths validated against task definition
            - File paths checked for security violations
            - Type checking enforced for parameter values

        Architecture:
            - Uses CPluginScript + dbHandler for proper DB sync
            - Shared with set_job_parameter management command
            - Consistent Result[T] pattern for error handling
        """
        try:
            form_data = json.loads(request.body.decode("utf-8"))
            job = models.Job.objects.get(id=pk)
            object_path = form_data["object_path"]
            value = form_data["value"]

            # Use unified utility (CPluginScript architecture)
            from ..lib.utils.parameters.set_param import set_parameter as set_job_param

            result = set_job_param(job, object_path, value)

            if result.success:
                return api_success(result.data)
            else:
                return api_error(result.error, status=400, details=result.error_details)

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(f"Job not found: {str(err)}", status=404)
        except Exception as err:
            logger.exception("Unexpected error setting parameter for job %s", pk, exc_info=err)
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def get_parameter(self, request, pk=None):
        """
        Get a parameter value from the job's configuration.

        Retrieves a specific parameter from the job's current configuration,
        allowing inspection of job setup before or after execution.

        Uses the unified CPluginScript architecture for consistent behavior
        with set_parameter and CLI commands.

        Args:
            request (Request): HTTP request with query parameter:
                - object_path (str): Path to the parameter to retrieve
            pk (int): Primary key of the job

        Returns:
            JsonResponse: Status and parameter information

        Query Parameters:
            ?object_path=inputData.XYZIN

        Response Format:
            {
                "status": "Success",
                "data": {
                    "path": "inputData.XYZIN",
                    "value": "/path/to/file.pdb",
                    "object_type": "CDataFile",
                    "file_path": "/path/to/file.pdb",
                    "db_file_id": "uuid-string",
                    "base_name": "file.pdb"
                }
            }

        Example:
            GET /api/jobs/123/get_parameter/?object_path=inputData.NCYCLES

        Architecture:
            - Uses CPluginScript for proper object hierarchy
            - Shared with get_job_parameter management command
            - Consistent Result[T] pattern for error handling
        """
        try:
            job = models.Job.objects.get(id=pk)
            object_path = request.GET.get("object_path")

            if not object_path:
                return api_error("Missing required parameter: object_path", status=400)

            # Use unified utility (CPluginScript architecture)
            from ..lib.utils.parameters.get_param import get_parameter as get_job_param

            result = get_job_param(job, object_path)

            if result.success:
                return api_success(result.data)
            else:
                return api_error(result.error, status=400, details=result.error_details)

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(f"Job not found: {str(err)}", status=404)
        except Exception as err:
            logger.exception("Unexpected error getting parameter for job %s", pk, exc_info=err)
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def upload_file_param(self, request, pk=None):
        """
        Upload and set a file parameter for the job.

        Handles file upload and associates it with a specific job parameter,
        managing file storage and parameter configuration automatically.

        Args:
            request (Request): Multipart HTTP request with file upload
            pk (int): Primary key of the job

        Returns:
            JsonResponse: Status and updated parameter information

        Request Format:
            - Multipart form data with file attachment
            - Additional fields for parameter configuration

        Response Format:
            {
                "status": "Success",
                "updated_item": {
                    "parameter_name": "XYZIN",
                    "file_path": "/uploaded/file/path",
                    "file_size": 12345
                }
            }

        Example:
            POST /api/jobs/123/upload_file_param/
            Content-Type: multipart/form-data

            - Secure file handling with validation
            - Atomic upload operations
        """
        job = models.Job.objects.get(id=pk)
        try:
            result = upload_file_param(job, request)
            return api_success({"updated_item": result})
        except CCP4ErrorHandling.CException as err:
            error_tree = getEtree(err)
            ET.indent(error_tree, " ")
            return api_error(ET.tostring(error_tree).decode("utf-8"), status=400)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.JobSerializer,
    )
    def preview(self, request, pk=None):
        """
        Launch external viewer for job results preview.

        Opens job results in an external molecular viewer application
        for interactive analysis and visualization.

        Args:
            request (Request): HTTP request with viewer specification
            pk (int): Primary key of the job

        Returns:
            Response: Success status or error message

        Request Data:
            - viewer (str): Name of the viewer application

        Supported Viewers:
            - CCP4mg: Molecular graphics program
            - Coot: Electron density fitting tool
            - PyMOL: Molecular visualization system

        Example:
            POST /api/jobs/123/preview/
            {
                "viewer": "ccp4mg"
            }

        Note:
            This endpoint triggers external process launch and requires
            appropriate viewer software installation.
        """
        try:
            the_job = models.Job.objects.get(id=pk)
            the_viewer = request.data.get("viewer")
            preview_job(the_viewer, str(the_job.directory))
            return api_success({"previewed": True})
        except models.File.DoesNotExist as err:
            logging.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def files(self, request, pk=None):
        """
        Retrieve files associated with a specific job.

        Returns a complete list of input and output files for the job,
        including file metadata and access information.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            Response: Serialized list of file objects

        Response Format:
            [
                {
                    "id": 789,
                    "uuid": "file-uuid",
                    "filename": "structure.pdb",
                    "file_type": "PDB",
                    "size": 12345,
                    "job_param_name": "XYZIN",
                    "is_input": true,
                    "created_date": "2024-01-01T12:00:00Z"
                }
            ]

        File Categories:
            - Input files: User-provided data
            - Output files: Generated results
            - Intermediate files: Processing artifacts
            - Log files: Execution records

        Side Effects:
            - Updates project last_access timestamp
            - Enables project activity tracking

        Example:
            GET /api/jobs/123/files/
        """
        job = models.Job.objects.get(pk=pk)
        serializer = serializers.FileSerializer(
            models.File.objects.filter(job=job), many=True
        )
        job.project.last_access = datetime.datetime.now(tz=timezone("UTC"))
        job.project.save()
        return Response(serializer.data)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def export_job(self, request, pk=None):
        """
        Export job and its dependencies as a ZIP archive.

        Creates a downloadable ZIP archive containing the job data,
        associated files, and dependency information for backup,
        sharing, or migration purposes.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job to export

        Returns:
            FileResponse: ZIP archive download or error response

        Archive Contents:
            - Job configuration and parameters
            - Input and output files
            - Reports and logs
            - Dependency information
            - Project metadata

        Response Headers:
            - Content-Type: application/zip
            - Content-Disposition: attachment
            - X-Export-Info: JSON metadata about export

        Export Metadata:
            - Project information
            - Job details
            - Export timestamp
            - File inventory

        Example:
            GET /api/jobs/123/export_job/

        Download:
            - Filename: {project_name}_job_{number}_{title}.zip
            - Streaming download for large archives
            - Automatic cleanup of temporary files

            - Streaming download optimized for bandwidth
            - Temporary file cleanup handled automatically

        Error Handling:
            - Validates job existence
            - Handles file access permissions
            - Manages storage space constraints
            - Provides detailed error messages
        """
        try:
            # Get the job and infer the project
            job = models.Job.objects.get(id=pk)
            project = job.project

            # Create job selection set (just this job number as string)
            job_selection = {str(job.number)}

            # Create temporary file for the export
            with tempfile.NamedTemporaryFile(suffix=".zip", delete=False) as temp_file:
                temp_path = Path(temp_file.name)

            try:
                # Export the project with job selection
                result_path = export_project_to_zip(
                    project=project, output_path=temp_path, job_selection=job_selection
                )

                # Generate a meaningful filename
                safe_project_name = "".join(
                    c for c in project.name if c.isalnum() or c in "._-"
                )
                safe_job_title = "".join(
                    c
                    for c in (job.title or job.task_name or f"job_{job.number}")
                    if c.isalnum() or c in "._-"
                )
                filename = f"{safe_project_name}_job_{job.number}_{safe_job_title}.zip"

                # Create file response
                def file_iterator():
                    with open(result_path, "rb") as f:
                        yield from f
                    # Clean up temp file after reading
                    try:
                        os.unlink(result_path)
                    except OSError:
                        pass

                response = FileResponse(
                    file_iterator(),
                    as_attachment=True,
                    filename=filename,
                    content_type="application/zip",
                )

                # Add custom headers with export info
                response["X-Export-Info"] = json.dumps(
                    {
                        "project_name": project.name,
                        "project_uuid": str(project.uuid),
                        "job_number": job.number,
                        "job_title": job.title or job.task_name,
                        "job_uuid": str(job.uuid),
                        "export_type": "job_selection",
                    }
                )

                logger.info(
                    f"Exported project archive for job {job.number} in project {project.name}"
                )
                return response

            except Exception as e:
                # Clean up temp file on error
                try:
                    os.unlink(temp_path)
                except OSError:
                    pass

                logger.exception("Export failed for job %s", pk, exc_info=e)
                return api_error(f"Export failed: {str(e)}", status=500)

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)
        except Exception as err:
            logger.exception(
                "Unexpected error during export for job %s", pk, exc_info=err
            )
            return api_error(f"Unexpected error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def export_job_file_menu(self, request, pk=None):
        """
        Retrieve file menu data for job export functionality.

        Gets the task-specific file menu configuration from the CCP4 Task Manager,
        which provides information about exportable files and their categories
        for the job's task type.

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the job

        Returns:
            JsonResponse: File menu configuration data or error message

        Response Format:
            {
                "status": "Success",
                "result": {
                    "menu_items": [...],
                    "file_categories": [...],
                    "export_options": [...]
                }
            }

        Menu Data Content:
            - Available file types for export
            - File category groupings
            - Export format options
            - Task-specific file handling rules

        Example:
            GET /api/jobs/123/export_job_file_menu/

            - Optimized for distributed computing environments

        Error Handling:
            - Validates job existence
            - Handles CCP4 Task Manager initialization errors
            - Provides detailed error messages for debugging
        """
        try:
            # Retrieve the job object
            job = models.Job.objects.get(id=pk)
            task_name = job.task_name

            logger.debug(
                "Retrieving file menu for job %s with task_name: %s", pk, task_name
            )

            # Get the task manager instance and retrieve menu data
            task_manager = CCP4TaskManager.TASKMANAGER()

            # TODO: Modernize exportJobFiles - legacy method not available in modern TaskManager
            # For now, return empty menu structure
            if hasattr(task_manager, 'exportJobFiles'):
                menu_result = task_manager.exportJobFiles(
                    taskName=task_name, jobId=job.uuid, mode="menu"
                )
            else:
                # Return stub menu structure for modern code
                logger.warning(
                    "exportJobFiles not available in TaskManager - returning empty menu"
                )
                menu_result = []

            # The task manager should return menu data - exact structure depends on CCP4 implementation
            # This might need adjustment based on the actual return type from TASKMANAGER
            menu_result = menu_result if menu_result else []

            logger.debug(
                "Successfully retrieved file menu for job %s, task %s", pk, task_name
            )

            return api_success({"result": menu_result})

        except models.Job.DoesNotExist as err:
            logger.exception("Failed to retrieve job with id %s", pk, exc_info=err)
            return api_error(f"Job not found: {str(err)}", status=404)

        except Exception as err:
            logger.exception(
                "Failed to get file menu for job %s with task_name %s",
                pk,
                (
                    getattr(job, "task_name", "unknown")
                    if "job" in locals()
                    else "unknown"
                ),
                exc_info=err,
            )
            return api_error(f"Task manager error: {str(err)}", status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
    )
    def export_job_file(self, request, pk=None):
        """
        Export a specific file from a job using CCP4 Task Manager.

        Retrieves and exports a specific file associated with the job using the
        CCP4 Task Manager's exportJobFiles functionality. The exported file is
        returned as a downloadable response.

        Args:
            request (Request): HTTP request object with query parameters:
                - mode (str): Export mode parameter for Task Manager
            pk (int): Primary key of the job

        Returns:
            FileResponse: Downloaded file or error response

        Query Parameters:
            - mode: Export mode (e.g., "pdb", "mtz", "log", etc.)

        Response:
            - Success: File download with appropriate headers
            - Error: JSON response with error details

        Example:
            GET /api/jobs/123/export_job_file/?mode=pdb
        """
        print("In export job file")
        print(f"Request: pk={pk}, mode={request.GET.get('mode')}")
        # Get the export mode from query parameters
        export_mode = request.GET.get("mode")
        if not export_mode:
            print("Missing export mode parameter")
            return api_error("Missing required 'mode' parameter", status=400)

        print(f"About to call utility function with pk={pk}, export_mode={export_mode}")
        # Import the utility function locally to avoid unused import lint error
        from ..lib.utils.jobs.export import export_job_file

        # Use the refactored utility function
        file_response, error_response = export_job_file(pk, export_mode)
        print(
            f"Utility function returned: file_response={file_response is not None}, error_response={error_response is not None}"
        )

        # Return either the file response or error response
        if error_response:
            print(f"Returning error response: {error_response}")
            return error_response
        print(f"Returning file response: {file_response}")
        return file_response
