import logging
from rest_framework.viewsets import ModelViewSet
from rest_framework.parsers import FormParser, MultiPartParser, JSONParser
from . import serializers
from ..db import models

logger = logging.getLogger(f"ccp4i2:{__name__}")


class ProjectTagViewSet(ModelViewSet):
    queryset = models.ProjectTag.objects.all()
    serializer_class = serializers.ProjectTagSerializer
    parser_classes = [JSONParser, FormParser, MultiPartParser]
