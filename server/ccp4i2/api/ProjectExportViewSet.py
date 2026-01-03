import logging
from django.http import FileResponse
from rest_framework.parsers import MultiPartParser, JSONParser, FormParser
from rest_framework.viewsets import ModelViewSet
from rest_framework.decorators import action
from rest_framework.response import Response
from . import serializers
from ..db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


class ProjectExportViewSet(ModelViewSet):
    queryset = models.ProjectExport.objects.all()
    serializer_class = serializers.ProjectExportSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]

    @action(
        detail=True,
        methods=["get"],
        permission_classes=[],
    )
    def download(self, _request, pk=None):
        export = self.get_object()
        project = export.project

        # Construct the expected file path based on the export creation logic
        from django.utils.text import slugify
        import os

        project_name = slugify(project.name or f"project_{project.id}")
        timestamp = export.time.strftime("%Y%m%d_%H%M%S")
        export_file_name = f"{project_name}_export_{timestamp}.ccp4_project.zip"
        export_file_path = os.path.join(
            project.directory, "CCP4_EXPORT_FILES", export_file_name
        )

        if os.path.exists(export_file_path):
            # FileResponse handles Content-Length and file closing automatically
            # Use as_attachment=True for proper download headers
            return FileResponse(
                open(export_file_path, "rb"),
                as_attachment=True,
                filename=export_file_name,
            )
        else:
            return Response({"error": "Export file not found"}, status=404)

    def destroy(self, request, *args, **kwargs):
        export = self.get_object()
        project = export.project

        print(f"In delete export for ProjectExport {export.id}")
        print(f"Project directory: {project.directory}")
        print(f"Export timestamp: {export.time}")

        # Construct the expected file paths based on the export creation logic
        from django.utils.text import slugify
        import os

        project_name = slugify(project.name or f"project_{project.id}")
        timestamp = export.time.strftime("%Y%m%d_%H%M%S")
        export_file_name = f"{project_name}_export_{timestamp}.ccp4_project.zip"
        export_dir = os.path.join(project.directory, "CCP4_EXPORT_FILES")
        export_file_path = os.path.join(export_dir, export_file_name)

        print(f"Constructed export file path: {export_file_path}")
        print(f"Export file exists: {os.path.exists(export_file_path)}")

        # Handle potential counter suffixes (same logic as export method)
        counter = 1
        base_name = export_file_name
        actual_export_path = export_file_path
        while (
            not os.path.exists(actual_export_path) and counter < 100
        ):  # Reasonable limit
            name_without_ext = base_name.rsplit(".", 1)[0]
            export_file_name = f"{name_without_ext}_{counter}.ccp4_project.zip"
            actual_export_path = os.path.join(export_dir, export_file_name)
            counter += 1

        print(f"Final export file path after counter check: {actual_export_path}")

        # Create log file path with same base name but .export.log extension
        log_file_name = export_file_name.replace(".ccp4_project.zip", ".export.log")
        log_file_path = os.path.join(export_dir, log_file_name)

        print(f"Log file path: {log_file_path}")
        print(f"Log file exists: {os.path.exists(log_file_path)}")

        # Delete the files if they exist
        files_deleted = []
        if os.path.exists(actual_export_path):
            try:
                os.remove(actual_export_path)
                files_deleted.append(actual_export_path)
                print(f"Successfully deleted export file: {actual_export_path}")
                logger.info("Deleted export file: %s", actual_export_path)
            except OSError as e:
                print(f"Failed to delete export file {actual_export_path}: {e}")
                logger.warning(
                    "Failed to delete export file %s: %s", actual_export_path, e
                )

        if os.path.exists(log_file_path):
            try:
                os.remove(log_file_path)
                files_deleted.append(log_file_path)
                print(f"Successfully deleted log file: {log_file_path}")
                logger.info("Deleted export log file: %s", log_file_path)
            except OSError as e:
                print(f"Failed to delete log file {log_file_path}: {e}")
                logger.warning(
                    "Failed to delete export log file %s: %s", log_file_path, e
                )

        print(f"Total files deleted: {len(files_deleted)}")
        if files_deleted:
            logger.info(
                "Deleted %d files for ProjectExport %s", len(files_deleted), export.id
            )

        # Call parent destroy method to delete the database record
        return super().destroy(request, *args, **kwargs)
