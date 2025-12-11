import logging
from django.http import FileResponse, JsonResponse
from rest_framework.response import Response
from rest_framework import status
from rest_framework.parsers import MultiPartParser, JSONParser

# Modern utilities
from ..lib.utils.files.digest import digest_file
from ..lib.utils.containers.json_encoder import CCP4i2JsonEncoder

# Modern utilities
from ..lib.utils.files.preview import preview_file
from xml.etree import ElementTree as ET
from rest_framework.parsers import FormParser, MultiPartParser
from rest_framework.viewsets import ModelViewSet
from rest_framework.decorators import action
from rest_framework.response import Response
from . import serializers
from ..db import models
from ..lib.response import api_success, api_error

logger = logging.getLogger(f"ccp4i2:{__name__}")


class FileViewSet(ModelViewSet):
    queryset = models.File.objects.all()
    serializer_class = serializers.FileSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]

    @action(
        detail=True,
        methods=["get"],
        permission_classes=[],
        serializer_class=serializers.FileSerializer,
    )
    def by_uuid(self, request, pk=None):
        try:
            the_file = models.File.objects.get(uuid=pk)
            serializer = serializers.FileSerializer(the_file, many=False)
            return Response(serializer.data)
        except models.File.DoesNotExist as err:
            logging.exception("Failed to retrieve file with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)

    @action(
        detail=True,
        methods=["get"],
        permission_classes=[],
        serializer_class=serializers.FileSerializer,
    )
    def download(self, request, pk=None):
        try:
            the_file = models.File.objects.get(id=pk)
            return FileResponse(open(the_file.path, "rb"), filename=the_file.name)
        except models.File.DoesNotExist as err:
            logging.exception("Failed to retrieve file with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)

    @action(
        detail=True,
        methods=["get"],
        permission_classes=[],
        serializer_class=serializers.FileSerializer,
    )
    def download_by_uuid(self, request, pk=None):
        try:
            the_file = models.File.objects.get(uuid=pk)
            return FileResponse(open(the_file.path, "rb"), filename=the_file.name)
        except models.File.DoesNotExist as err:
            logging.exception("Failed to retrieve file with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)

    @action(
        detail=True,
        methods=["get"],
        permission_classes=[],
        serializer_class=serializers.FileSerializer,
    )
    def digest(self, request, pk=None):
        """
        Get digest (summary) of file contents.

        Uses modern digest utility that:
        - Loads CDataFileContent subclass based on file type
        - Returns JSON dictionary representation
        - Tries to infer file type if not explicit

        Returns dict with file metadata like sequences, composition, cell, etc.
        """
        try:
            the_file = models.File.objects.get(id=pk)
            result = digest_file(the_file)

            # Check if digest returned an error
            if isinstance(result, dict) and result.get("status") == "Failed":
                return api_error(result.get("reason", "Digest failed"), status=400)

            return api_success(result)

        except models.File.DoesNotExist:
            logger.exception("File %s not found", pk)
            return api_error("File not found", status=404)
        except Exception as err:
            logger.exception("Failed to digest file %s", pk, exc_info=err)
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["get"],
        permission_classes=[],
        serializer_class=serializers.FileSerializer,
    )
    def digest_by_uuid(self, request, pk=None):
        """
        Get digest (summary) of file contents by UUID.

        Uses modern digest utility that:
        - Loads CDataFileContent subclass based on file type
        - Returns JSON dictionary representation
        - Tries to infer file type if not explicit

        Returns dict with file metadata like sequences, composition, cell, etc.
        """
        try:
            the_file = models.File.objects.get(uuid=pk)
            result = digest_file(the_file)

            # Check if digest returned an error
            if isinstance(result, dict) and result.get("status") == "Failed":
                return api_error(result.get("reason", "Digest failed"), status=400)

            return api_success(result)

        except models.File.DoesNotExist:
            logger.exception("File with UUID %s not found", pk)
            return api_error("File not found", status=404)
        except Exception as err:
            logger.exception("Failed to digest file with UUID %s", pk, exc_info=err)
            return api_error(str(err), status=500)

    @action(
        detail=True,
        methods=["post"],
        permission_classes=[],
        serializer_class=serializers.FileSerializer,
    )
    def preview(self, request, pk=None):
        try:
            the_file = models.File.objects.get(id=pk)
            the_viewer = request.data.get("viewer")
            preview_file(the_viewer, str(the_file.path))
            return api_success({"previewed": True})
        except models.File.DoesNotExist as err:
            logging.exception("Failed to retrieve file with id %s", pk, exc_info=err)
            return api_error(str(err), status=404)
