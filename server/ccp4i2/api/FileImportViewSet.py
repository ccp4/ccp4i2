import logging
from rest_framework.parsers import MultiPartParser, JSONParser, FormParser

from rest_framework.viewsets import ModelViewSet
from . import serializers
from ..db import models

logger = logging.getLogger(f"ccp4x:{__name__}")


class FileImportViewSet(ModelViewSet):
    queryset = models.FileImport.objects.all()
    serializer_class = serializers.FileImportSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]
