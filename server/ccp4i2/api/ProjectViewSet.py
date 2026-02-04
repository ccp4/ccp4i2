import logging
import datetime
import json
import pathlib
import os
import subprocess
from asgiref.sync import async_to_sync
from pytz import timezone
from django.http import Http404
from django.http import FileResponse
from django.core.management import call_command
from rest_framework.response import Response
from rest_framework import status
from rest_framework.parsers import MultiPartParser, JSONParser

# Modern utilities
from ..lib.async_create_job import create_job_async

# Modern utilities
from ..lib.utils.navigation.list_project import list_project
from ..lib.utils.navigation.task_tree import get_task_tree
from ..lib.utils.files.preview import preview_file
from ..lib.utils.files.resolve_fileuse import resolve_fileuse, is_fileuse_pattern

from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.viewsets import ModelViewSet
from rest_framework.decorators import action
from rest_framework.response import Response
from rest_framework.permissions import IsAuthenticated
from django.db.models import Prefetch
from . import serializers
from ..db import models
from ..lib.utils.navigation.dependencies import delete_job_and_dependents
from django.http import JsonResponse
from django.conf import settings
from django.utils.text import slugify
from ..lib.response import api_success, api_error

logger = logging.getLogger(f"ccp4i2:{__name__}")


class ProjectViewSet(ModelViewSet):
    """
    ProjectViewSet

    This viewset provides various actions to interact with the `Project` model and its related entities.
    It includes endpoints for importing projects, retrieving associated files, jobs, and metadata, as well
    as handling project-specific operations like directory listing, file retrieval, and task creation.

    Actions:
        - import_project: Handles the upload and import of project files (ZIP format).
        - files: Retrieves a list of files associated with a specific project.
        - jobs: Retrieves a list of jobs associated with a specific project.
        - job_float_values: Retrieves all `JobFloatValue` instances associated with a specific project.
        - job_char_values: Retrieves job characteristic values for a specific project.
        - tags: Retrieves tags associated with a specific project.
        - directory: Retrieves the directory listing of a specific project.
        - project_file: Retrieves a specific file from the project's directory.
        - preview_file: Previews a specific file from the project's directory using a specified viewer.
        - task_tree: Retrieves the task tree structure (not directly tied to a specific project).
        - create_task: Creates a new task (job) within a specific project.

    Attributes:
        - queryset: The queryset used to retrieve `Project` objects.
        - serializer_class: The default serializer class for the viewset.
        - parser_classes: The parsers allowed for handling incoming data.

    Notes:
        - The `import_project` action ensures secure file handling and validates file types.
        - Several actions update the `last_access` timestamp of the project to track usage.
        - Directory traversal attacks are mitigated in file-related actions by validating file paths.
        - Logging is used extensively to capture errors and important events.

    Security:
        - All endpoints require authentication via Azure AD middleware.
        - The IsAuthenticated permission class provides an additional layer of defense.
    """

    queryset = models.Project.objects.all()
    serializer_class = serializers.ProjectSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    permission_classes = [IsAuthenticated]

    def get_queryset(self):
        """Optimize queryset based on action."""
        queryset = super().get_queryset()
        if self.action == 'list':
            # List view: order by most recent first, prefetch tags
            return queryset.prefetch_related('tags').order_by('-last_access')
        else:
            # Detail view: include tags
            return queryset.prefetch_related('tags')

    def get_serializer_class(self):
        """Use lightweight serializer for list view."""
        if self.action == 'list':
            return serializers.ProjectListSerializer
        return serializers.ProjectSerializer

    def destroy(self, request, *args, **kwargs):
        try:
            instance = self.get_object()
            logger.warning("Deleting project %s", instance)
            while (
                models.Job.objects.filter(project=instance, parent__isnull=True).count()
                > 0
            ):
                last_job = models.Job.objects.filter(
                    project=instance, parent__isnull=True
                ).last()
                logger.warning("Deleting job %s", last_job)
                delete_job_and_dependents(last_job)

            # Attempt some security by using a defined list of subdirectories for deletion
            for subdir in [
                "CCP4_COOT",
                "CCP4_IMPORTED_FILES",
                "CCP4_JOBS",
                "CCP4_PROJECT_FILES",
                "CCP4_TMP",
            ]:
                subdir_path = os.path.join(instance.directory, subdir)
                if os.path.exists(subdir_path):
                    logger.warning("Deleting subdirectory %s", subdir_path)
                    # Definitely need to
                    for root, dirs, files in os.walk(
                        subdir_path, topdown=False, followlinks=False
                    ):
                        for file in files:
                            file_path = os.path.join(root, file)
                        for directory in dirs:
                            dir_path = os.path.join(root, directory)
                            try:
                                os.rmdir(dir_path)
                            except Exception as e:
                                logger.warning(
                                    "Failed to delete directory %s: %s", dir_path, e
                                )
                            except Exception as e:
                                logger.warning(
                                    "Failed to delete file %s: %s", file_path, e
                                )
                    try:
                        os.rmdir(subdir_path)
                    except Exception as e:
                        logger.warning(
                            "Failed to delete directory %s: %s", subdir_path, e
                        )
            # Attempt to delete any special files in the project directory
            for special_file in [".DS_Store"]:
                special_file_path = os.path.join(instance.directory, special_file)
                if os.path.exists(special_file_path):
                    logger.warning("Deleting special file %s", special_file_path)
                    try:
                        os.remove(special_file_path)
                    except Exception as e:
                        logger.warning(
                            "Failed to delete file %s: %s", special_file_path, e
                        )
            # Attempt to delete the main project directory
            try:
                os.rmdir(instance.directory)
            except Exception as e:
                logger.warning(
                    "Failed to delete project directory %s: %s", instance.directory, e
                )
            instance.delete()
            logger.warning("Deleted project %s", instance)

            # Note I am adding a bit of body to the response because of an odd
            # javascript feature which presents as network error if no body in response.
            return api_success({"deleted": True})
        except Http404:
            return Http404("Job not found")

    @action(
        detail=False,
        methods=["post"],
        parser_classes=[
            JSONParser,
            FormParser,
            MultiPartParser,
        ],  # Allow JSON and form data
    )
    def import_project(self, request):
        uploaded_files = request.FILES.getlist("files")
        if not uploaded_files:
            return api_error("No file provided", status=400)
        # Define a secure, platform-independent storage path
        secure_storage_dir = pathlib.Path(settings.MEDIA_ROOT) / "uploaded_files"
        secure_storage_dir.mkdir(parents=True, exist_ok=True)

        # Save the file to the secure storage directory
        for uploaded_file in uploaded_files:
            if not uploaded_file.name.endswith(".zip"):
                return api_error("Invalid file type", status=400)
            # Ensure the filename is safe
            # if not uploaded_file.name.isalnum():
            #    return Response(
            #        {"status": "Failed", "reason": "Invalid file name"},
            #        status=status.HTTP_400_BAD_REQUEST,
            #    )
            file_path = secure_storage_dir / slugify(uploaded_file.name)
            with open(file_path, "wb") as destination:
                for chunk in uploaded_file.chunks():
                    destination.write(chunk)
            try:
                call_command("import_ccp4_project_zip", str(file_path), "--detach")
            except Exception as e:
                logger.exception(
                    "Failed to import project from %s", file_path, exc_info=e
                )
                return api_error(str(e), status=500)
            logger.warning("File uploaded and saved to %s", file_path)

        return api_success({"imported": True})

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def files(self, request, pk=None):
        """
        Retrieve a list of files associated with a specific project.

        Parameters:

            request (Request): The HTTP request object.
            pk (int, optional): The primary key of the project.

        Returns:
            Response: A Response object containing serialized file data.
        """

        project = models.Project.objects.get(pk=pk)
        serializer = serializers.FileSerializer(
            models.File.objects.filter(job__project=project), many=True
        )
        project.last_access = datetime.datetime.now(tz=timezone("UTC"))
        project.save()
        return Response(serializer.data)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.FileSerializer,
    )
    def file_uses(self, request, pk=None):
        """
        Retrieve a list of file_uses associated with a specific project.

        Parameters:

            request (Request): The HTTP request object.
            pk (int, optional): The primary key of the project.

        Returns:
            Response: A Response object containing serialized file_use data.
        """

        project = models.Project.objects.get(pk=pk)
        serializer = serializers.FileUseSerializer(
            models.FileUse.objects.filter(job__project=project), many=True
        )
        project.last_access = datetime.datetime.now(tz=timezone("UTC"))
        project.save()
        return Response(serializer.data)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobSerializer,
        parser_classes=[
            JSONParser,
            FormParser,
            MultiPartParser,
        ],  # Allow JSON and form data
    )
    def jobs(self, request, pk=None):
        """
        Retrieve a list of jobs associated with a specific project.
        Args:
            request (Request): The HTTP request object.
            pk (int, optional): The primary key of the project.
        Returns:
            Response: A Response object containing serialized job data.
        """

        project = models.Project.objects.get(pk=pk)
        serializer = serializers.JobSerializer(
            models.Job.objects.filter(project=project), many=True
        )
        project.last_access = datetime.datetime.now(tz=timezone("UTC"))
        project.save()
        return Response(serializer.data)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobFloatValueSerializer,
    )
    def job_float_values(self, request, pk=None):
        """
        Retrieve all JobFloatValue instances associated with a specific project.
        Args:
            request (Request): The HTTP request object.
            pk (int, optional): The primary key of the project.
        Returns:
            Response: A Response object containing serialized data of JobFloatValue instances.
        """

        project = models.Project.objects.get(pk=pk)
        serializer = serializers.JobFloatValueSerializer(
            models.JobFloatValue.objects.filter(job__project=project), many=True
        )
        project.last_access = datetime.datetime.now(tz=timezone("UTC"))
        project.save()
        return Response(serializer.data)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.JobCharValueSerializer,
    )
    def job_char_values(self, request, pk=None):
        """
        Retrieve job characteristic values for a specific project.
        Args:
            request (Request): The HTTP request object.
            pk (int, optional): The primary key of the project.
        Returns:
            Response: A Response object containing serialized job characteristic values.
        """

        project = models.Project.objects.get(pk=pk)
        serializer = serializers.JobCharValueSerializer(
            models.JobCharValue.objects.filter(job__project=project), many=True
        )
        project.last_access = datetime.datetime.now(tz=timezone("UTC"))
        project.save()
        return Response(serializer.data)

    @action(
        detail=True,
        methods=["get"],
    )
    def job_tree(self, request, pk=None):
        """
        Consolidated endpoint for job tree view.

        Returns hierarchical job data with embedded files and KPIs in a single request.
        This replaces 4 separate calls to: jobs/, files/, job_float_values/, job_char_values/

        The response is structured as a tree where each job includes:
        - All job fields (id, uuid, number, title, status, etc.)
        - files: Array of files associated with this job
        - kpis: Object with float_values and char_values for this job
        - children: Array of child jobs (recursively structured the same way)

        Jobs are sorted by job number (highest first) at each level.

        Args:
            request (Request): The HTTP request object.
            pk (int, optional): The primary key of the project.

        Returns:
            Response: A Response object containing the hierarchical job tree.
        """
        try:
            project = models.Project.objects.get(pk=pk)

            # Single optimized query with all related data prefetched
            # Note: related_name values from models.py:
            #   File.job -> related_name="files"
            #   JobFloatValue.job -> related_name="float_values"
            #   JobCharValue.job -> related_name="char_values"
            jobs = list(
                models.Job.objects.filter(project=project)
                .prefetch_related(
                    Prefetch(
                        "files",
                        queryset=models.File.objects.all(),
                    ),
                    # JobValueKey.name is the primary key, so we just need key_id (no select_related needed)
                    Prefetch(
                        "float_values",
                        queryset=models.JobFloatValue.objects.all(),
                    ),
                    Prefetch(
                        "char_values",
                        queryset=models.JobCharValue.objects.all(),
                    ),
                )
                .select_related("parent")
            )

            # Build lookup dictionaries for O(1) access
            jobs_by_id = {job.id: job for job in jobs}
            children_by_parent = {}
            for job in jobs:
                parent_id = job.parent_id
                if parent_id not in children_by_parent:
                    children_by_parent[parent_id] = []
                children_by_parent[parent_id].append(job)

            def parse_job_number(job_number: str) -> list:
                """Parse job number like '1.2.3' into [1, 2, 3] for sorting."""
                try:
                    return [int(part) for part in job_number.split(".")]
                except (ValueError, AttributeError):
                    return [0]

            def build_job_node(job):
                """Build a single job node with embedded files and KPIs."""
                # Serialize job data
                job_data = serializers.JobSerializer(job).data

                # Add embedded files (uses prefetched 'files' relation)
                job_data["files"] = serializers.FileSerializer(
                    job.files.all(), many=True
                ).data

                # Add embedded KPIs as a structured object
                # JobValueKey.name is the primary key, so kv.key_id gives us the string name directly
                float_vals = {kv.key_id: kv.value for kv in job.float_values.all()}
                char_vals = {kv.key_id: kv.value for kv in job.char_values.all()}

                # Debug logging
                if float_vals or char_vals:
                    logger.info(f"Job {job.number} KPIs: float={float_vals}, char={char_vals}")

                job_data["kpis"] = {
                    "float_values": float_vals,
                    "char_values": char_vals,
                }

                # Recursively build children
                child_jobs = children_by_parent.get(job.id, [])
                # Sort children by job number (highest first)
                child_jobs.sort(key=lambda j: parse_job_number(j.number), reverse=True)
                job_data["children"] = [build_job_node(child) for child in child_jobs]

                return job_data

            # Build tree starting from root jobs (parent=None)
            root_jobs = children_by_parent.get(None, [])
            # Sort root jobs by job number (highest first)
            root_jobs.sort(key=lambda j: parse_job_number(j.number), reverse=True)
            job_tree = [build_job_node(job) for job in root_jobs]

            # Update last_access only once
            project.last_access = datetime.datetime.now(tz=timezone("UTC"))
            project.save()

            return Response(
                {
                    "job_tree": job_tree,
                    "total_jobs": len(jobs),
                    "total_files": sum(len(job.files.all()) for job in jobs),
                }
            )

        except models.Project.DoesNotExist:
            return api_error("Project not found", status=404)
        except Exception as e:
            logger.exception("Failed to build job tree for project %s", pk, exc_info=e)
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["get", "post"],
        serializer_class=serializers.ProjectTagSerializer,
    )
    def tags(self, request, pk=None):
        """
        Retrieve tags for a specific project (GET) or add a tag to a project (POST).

        GET: Returns list of tags for the project
        POST: Expects 'tag_id' in request data to add tag to project
        """
        project = models.Project.objects.get(pk=pk)

        if request.method == "GET":
            # Return project tags
            project_tag_serializer = serializers.ProjectTagSerializer(
                project.tags, many=True
            )
            project.last_access = datetime.datetime.now(tz=timezone("UTC"))
            project.save()
            return Response(project_tag_serializer.data)

        elif request.method == "POST":
            # Add tag to project
            try:
                tag_id = request.data.get("tag_id")

                if not tag_id:
                    return Response(
                        {"error": "tag_id is required"},
                        status=status.HTTP_400_BAD_REQUEST,
                    )

                tag = models.ProjectTag.objects.get(pk=tag_id)
                project.tags.add(tag)
                project.last_access = datetime.datetime.now(tz=timezone("UTC"))
                project.save()

                return Response(
                    {
                        "status": "success",
                        "message": f"Tag '{tag.text}' added to project",
                    },
                    status=status.HTTP_200_OK,
                )

            except models.ProjectTag.DoesNotExist:
                return Response(
                    {"error": "Tag not found"},
                    status=status.HTTP_404_NOT_FOUND,
                )
            except Exception as e:
                logger.exception("Failed to add tag to project", exc_info=e)
                return Response(
                    {"error": str(e)},
                    status=status.HTTP_500_INTERNAL_SERVER_ERROR,
                )

    @action(
        detail=True,
        methods=["post"],
        url_path="add_tag_by_text",
        serializer_class=serializers.ProjectTagSerializer,
    )
    def add_tag_by_text(self, request, pk=None):
        """
        Add a tag to a project by text (get-or-create semantics).

        POST: Expects 'text' in request data. Creates the tag if it doesn't exist,
              then adds it to the project.

        Returns the tag that was added.
        """
        try:
            project = models.Project.objects.get(pk=pk)
            text = request.data.get("text")

            if not text:
                return Response(
                    {"error": "text is required"},
                    status=status.HTTP_400_BAD_REQUEST,
                )

            # Get or create the tag (with no parent)
            tag, created = models.ProjectTag.objects.get_or_create(
                text=text,
                parent=None,
            )

            # Add tag to project
            project.tags.add(tag)
            project.last_access = datetime.datetime.now(tz=timezone("UTC"))
            project.save()

            return Response(
                {
                    "status": "success",
                    "tag": serializers.ProjectTagSerializer(tag).data,
                    "created": created,
                    "message": f"Tag '{tag.text}' added to project",
                },
                status=status.HTTP_200_OK,
            )

        except models.Project.DoesNotExist:
            return Response(
                {"error": "Project not found"},
                status=status.HTTP_404_NOT_FOUND,
            )
        except Exception as e:
            logger.exception("Failed to add tag to project by text", exc_info=e)
            return Response(
                {"error": str(e)},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.ProjectSerializer,
    )
    def directory(self, request, pk=None):
        """
        Handles the request to retrieve the directory listing of a project.

        Args:
            request (HttpRequest): The HTTP request object.
            pk (int, optional): The primary key of the project to retrieve the directory listing for.

        Returns:
            JsonResponse: A JSON response containing the status and the directory listing of the project.
                          If successful, the response contains the directory listing in the "container" field.
                          If a TypeError occurs during JSON encoding, the response contains an error message.

        Raises:
            Project.DoesNotExist: If no project with the given primary key exists.
        """

        the_project = models.Project.objects.get(pk=pk)
        result = list_project(str(the_project.uuid))
        try:
            result_str = json.dumps(result, indent=2)
            return JsonResponse(
                json.loads(f'{{"status": "Success", "container": {result_str}}}')
            )
        except TypeError as err:
            logger.exception("Failed encoding listing of %s", exc_info=err)
            return JsonResponse(
                {status: "Failed", "container": {"Reason": "TypeError"}}
            )

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.ProjectSerializer,
    )
    def project_file(self, request, pk=None):
        """
        Retrieve a file from the specified project directory.
        This view handles GET requests to retrieve a file from a project's directory.
        It ensures that the requested file path is within the project's directory to
        prevent directory traversal attacks.

        Args:
            request (HttpRequest): The HTTP request object containing query parameters.
            pk (int, optional): The primary key of the project.

        Raises:
            Http404: If the requested file path is not within the project's directory.

        Returns:
            FileResponse: A response object containing the requested file.
        """
        the_project = models.Project.objects.get(pk=pk)
        file_path = request.GET.get("path")
        logger.info("File path %s", file_path)
        composite_path: pathlib.Path = pathlib.Path(the_project.directory) / file_path
        if pathlib.Path(the_project.directory) not in composite_path.resolve().parents:
            raise Http404("Unacceptable file")

        return FileResponse(open(composite_path, "rb"), filename=composite_path.name)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.ProjectSerializer,
    )
    def preview_file(self, request, pk=None):
        """
        Preview a file from the specified project directory, using the specified viewer.
        This view handles POST requests to preview a file from a project's directory.
        It ensures that the requested file path is within the project's directory to
        prevent directory traversal attacks.

        Args:
            request (HttpRequest): The HTTP request object containing query parameters.
            pk (int, optional): The primary key of the project.

        Raises:
            Http404: If the requested file path is not within the project's directory.

        Returns:
            JsonResponse: A JSON response indicating success or failure of the viewer launch.
        """
        the_project = models.Project.objects.get(pk=pk)
        file_path = request.data.get("path")
        viewer = request.data.get("viewer")
        composite_path = pathlib.Path(the_project.directory) / (
            file_path if not file_path.startswith("/") else file_path[1:]
        )
        if not composite_path.resolve().is_relative_to(
            pathlib.Path(the_project.directory)
        ):
            raise Http404("Unacceptable file - outside project directory")

        logger.info("Previewing file %s with viewer %s", composite_path, viewer)
        try:
            # call_command("preview_file", "-d", "-e", viewer, "-p", str(composite_path))
            preview_file(viewer, str(composite_path))
            return api_success({"previewed": True})
        except Exception as err:
            logger.exception(
                "Failed to preview file %s with viewer %s",
                composite_path,
                viewer,
                exc_info=err,
            )
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.ProjectSerializer,
    )
    def create_task(self, request, pk=None):
        """
        Create a new task/job in the project.

        Modern implementation using async_create_job (same as management command).

        Request body should contain:
            - task_name: Plugin name (required, e.g., "ctruncate", "refmac")
            - title: Job title (optional)
            - parent_job_uuid: Parent job UUID for nested jobs (optional)
            - input_params: Dict of input parameters (optional)
            - context_job_uuid: Context job UUID for input population (optional).
                If not provided, uses the most recent job in the project.
            - auto_context: Whether to auto-select context job (default: True).
                Set to False to create a job without context-based input population.

        Returns:
            JSON with status and new job details
        """
        try:
            the_project = models.Project.objects.get(pk=pk)

            # Parse request body
            body = json.loads(request.body.decode("utf-8"))
            task_name = body.get("task_name")

            if not task_name:
                return api_error("task_name is required", status=400)

            # Extract optional parameters
            title = body.get("title")
            parent_job_uuid = body.get("parent_job_uuid")
            input_params = body.get("input_params")
            auto_context = body.get("auto_context", True)

            # Convert context_job_uuid string to UUID if provided
            context_job_uuid_str = body.get("context_job_uuid")
            context_job_uuid = None
            if context_job_uuid_str:
                import uuid as uuid_module
                context_job_uuid = uuid_module.UUID(context_job_uuid_str)

            # Use modern async job creation (same as management command)
            # Use async_to_sync to properly handle Django database connections
            result = async_to_sync(create_job_async)(
                project_uuid=the_project.uuid,
                task_name=task_name,
                title=title,
                parent_job_uuid=parent_job_uuid,
                save_params=True,
                input_params=input_params,
                context_job_uuid=context_job_uuid,
                auto_context=auto_context,
            )

            # Get the created job from database
            new_job = models.Job.objects.get(uuid=result['job_uuid'])
            serializer = serializers.JobSerializer(new_job)

            return api_success({"new_job": serializer.data})

        except models.Project.DoesNotExist:
            logger.exception("Project %s not found", pk)
            return api_error("Project not found", status=404)
        except KeyError as e:
            logger.exception("Missing required field: %s", e)
            return api_error(f"Missing required field: {e}", status=400)
        except Exception as e:
            logger.exception("Failed to create task", exc_info=e)
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["post"],
        serializer_class=serializers.ProjectSerializer,
    )
    def export(self, request, pk=None):
        the_project = models.Project.objects.get(pk=pk)

        # Generate unique filepath based on project name, rooted in project.directory
        project_name = slugify(the_project.name or f"project_{the_project.id}")
        project_export = models.ProjectExport.objects.create(
            project=the_project, time=datetime.datetime.now(tz=timezone("UTC"))
        )
        project_export.save()
        timestamp = project_export.time.strftime("%Y%m%d_%H%M%S")
        export_file_name = f"{project_name}_export_{timestamp}.ccp4_project.zip"

        # Use CCP4_EXPORT_FILES directory (excluded from exports to prevent recursive inclusion)
        export_dir = os.path.join(the_project.directory, "CCP4_EXPORT_FILES")
        os.makedirs(export_dir, exist_ok=True)

        export_file_path = os.path.join(export_dir, export_file_name)

        # Ensure the export file path doesn't already exist (add counter if needed)
        counter = 1
        base_name = export_file_name
        while os.path.exists(export_file_path):
            name_without_ext = base_name.rsplit(".", 1)[0]
            export_file_name = f"{name_without_ext}_{counter}.ccp4_project.zip"
            export_file_path = os.path.join(export_dir, export_file_name)
            counter += 1

        # Create log file path with same base name but .export.log extension
        log_file_name = export_file_name.replace(".ccp4_project.zip", ".export.log")
        log_file_path = os.path.join(export_dir, log_file_name)

        # Start subprocess to run export_project management command in background
        try:
            with open(log_file_path, "w") as log_file:
                process = subprocess.Popen(
                    [
                        "ccp4-python",
                        "manage.py",
                        "export_project",
                        "-pi",
                        str(the_project.id),
                        "-o",
                        export_file_path,
                    ],
                    stdout=log_file,
                    stderr=subprocess.STDOUT,
                    start_new_session=True,
                )

            return JsonResponse(
                {
                    "status": "Success",
                    "export_file_name": export_file_name,
                    "log_file_name": log_file_name,
                    "process_id": process.pid,
                }
            )
        except Exception as e:
            logger.exception(
                "Failed to start export process for project %s",
                the_project.id,
                exc_info=e,
            )
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["get"],
        serializer_class=serializers.ProjectExportSerializer,
    )
    def exports(self, request, pk=None):
        """
        Retrieve a list of project exports for a specific project.

        This action returns all ProjectExport instances associated with the given project,
        ordered by creation time (most recent first). This allows users to see the history
        of exports and their status.

        Args:
            request (HttpRequest): The HTTP request object.
            pk (int, optional): The primary key of the project.

        Returns:
            Response: A Response object containing serialized ProjectExport data,
                     including export file names, creation times, and status information.
        """
        try:
            project = models.Project.objects.get(pk=pk)

            # Get all exports for this project, ordered by most recent first
            project_exports = models.ProjectExport.objects.filter(
                project=project
            ).order_by("-time")

            # Update project last access time
            project.last_access = datetime.datetime.now(tz=timezone("UTC"))
            project.save()

            # Serialize the export data
            serializer = serializers.ProjectExportSerializer(project_exports, many=True)

            return Response(
                serializer.data,
            )

        except models.Project.DoesNotExist:
            return api_error("Project not found", status=404)
        except Exception as e:
            logger.exception(
                "Failed to retrieve exports for project %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)

    @action(
        detail=True,
        methods=["delete"],
        url_path=r"tags/(?P<tag_id>\d+)",
    )
    def remove_tag(self, request, pk=None, tag_id=None):
        """
        Remove a tag from a specific project.
        """
        try:
            project = models.Project.objects.get(pk=pk)
            tag = models.ProjectTag.objects.get(pk=tag_id)
            project.tags.remove(tag)
            project.last_access = datetime.datetime.now(tz=timezone("UTC"))
            project.save()

            return Response(
                {
                    "status": "success",
                    "message": f"Tag '{tag.text}' removed from project",
                },
                status=status.HTTP_200_OK,
            )

        except models.Project.DoesNotExist:
            return Response(
                {"error": "Project not found"},
                status=status.HTTP_404_NOT_FOUND,
            )
        except models.ProjectTag.DoesNotExist:
            return Response(
                {"error": "Tag not found"},
                status=status.HTTP_404_NOT_FOUND,
            )
        except Exception as e:
            logger.exception("Failed to remove tag from project", exc_info=e)
            return Response(
                {"error": str(e)},
                status=status.HTTP_500_INTERNAL_SERVER_ERROR,
            )

    @action(
        detail=True,
        methods=["get", "post"],
        serializer_class=serializers.ProjectSerializer,
    )
    def resolve_fileuse(self, request, pk=None):
        """
        Resolve a fileUse reference to file metadata.

        FileUse syntax allows referencing files from previous jobs:
            [-1].XYZOUT[0]              - First XYZOUT from most recent job
            prosmart_refmac[-1].XYZOUT  - XYZOUT from most recent prosmart_refmac job
            refmac[-2].HKLOUT[0]        - HKLOUT from second-to-last refmac job

        The jobIndex can be negative (counting from end) or positive (from start).

        GET: Query parameter 'fileuse' contains the reference string
        POST: JSON body with 'fileuse' key

        Args:
            request (Request): HTTP request object
            pk (int): Primary key of the project

        Returns:
            JsonResponse: File metadata on success:
                {
                    "status": "Success",
                    "data": {
                        "project": "uuid-without-hyphens",
                        "baseName": "filename.pdb",
                        "dbFileId": "file-uuid-without-hyphens",
                        "relPath": "CCP4_JOBS/job_001",
                        "fullPath": "/full/path/to/file.pdb"
                    }
                }

        Example:
            GET /api/projects/123/resolve_fileuse/?fileuse=[-1].XYZOUT[0]
            POST /api/projects/123/resolve_fileuse/
                {"fileuse": "prosmart_refmac[-1].XYZOUT"}
        """
        try:
            project = models.Project.objects.get(pk=pk)

            # Get fileuse string from request
            if request.method == "GET":
                fileuse = request.GET.get("fileuse")
            else:
                body = json.loads(request.body.decode("utf-8"))
                fileuse = body.get("fileuse")

            if not fileuse:
                return api_error(
                    "Missing required parameter: 'fileuse'. "
                    "Example: [-1].XYZOUT[0] or task_name[-1].PARAM",
                    status=400
                )

            # Resolve the fileuse reference
            result = resolve_fileuse(project, fileuse)

            if result.success:
                return api_success(result.data)
            else:
                return api_error(result.error, status=400)

        except models.Project.DoesNotExist:
            return api_error("Project not found", status=404)
        except json.JSONDecodeError as e:
            return api_error(f"Invalid JSON in request body: {e}", status=400)
        except Exception as e:
            logger.exception(
                "Failed to resolve fileuse for project %s", pk, exc_info=e
            )
            return api_error(str(e), status=500)
