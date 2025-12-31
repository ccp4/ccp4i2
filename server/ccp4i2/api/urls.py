from django.conf import settings
from django.conf.urls.static import static
from django.urls import include, path
from rest_framework import routers

from .ProjectExportViewSet import ProjectExportViewSet
from .ProjectViewSet import ProjectViewSet
from .ProjectTagViewSet import ProjectTagViewSet
from .FileViewSet import FileViewSet
from .JobViewSet import JobViewSet
from .FileTypeViewSet import FileTypeViewSet
from .FileImportViewSet import FileImportViewSet
from .FileUseViewSet import FileUseViewSet
from . import views

router = routers.DefaultRouter()
router.register("projects", ProjectViewSet)
router.register("projecttags", ProjectTagViewSet)
router.register("files", FileViewSet)
router.register("jobs", JobViewSet)
router.register("filetypes", FileTypeViewSet)
router.register("fileimports", FileImportViewSet)
router.register("fileuses", FileUseViewSet)
router.register("projectexports", ProjectExportViewSet)

urlpatterns = [
    path("", include(router.urls)),
    path("health/", views.health_check, name="health_check"),
    path("task_tree/", views.task_tree, name="task_tree"),
    path("active_jobs/", views.active_jobs, name="active_jobs"),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

# Include azure_extensions URLs if the app is installed
# This allows Azure-specific features (staged uploads) without polluting core ccp4i2
try:
    from django.apps import apps
    if apps.is_installed("azure_extensions"):
        from azure_extensions.urls import urlpatterns as azure_urls
        urlpatterns = urlpatterns + azure_urls
except ImportError:
    pass  # azure_extensions not available (non-Azure deployment)
