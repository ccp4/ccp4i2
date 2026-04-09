import logging
from rest_framework.parsers import MultiPartParser, JSONParser, FormParser
from rest_framework.permissions import IsAuthenticated
from rest_framework.viewsets import ModelViewSet
from . import serializers
from ..db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


class FileTypeViewSet(ModelViewSet):

    queryset = models.FileType.objects.all()
    serializer_class = serializers.FileTypeSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]
    permission_classes = [IsAuthenticated]
