"""
URL configuration for Azure Extensions app.

This module registers the staged upload endpoints under /api/uploads/.
"""

from django.urls import path, include
from rest_framework import routers

from .views import StagedUploadViewSet

router = routers.DefaultRouter()
router.register("uploads", StagedUploadViewSet)

urlpatterns = [
    path("", include(router.urls)),
]
