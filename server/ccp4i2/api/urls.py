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
from .StagedUploadViewSet import StagedUploadViewSet
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
router.register("uploads", StagedUploadViewSet)

urlpatterns = [
    path("", include(router.urls)),
    path("health/", views.health_check, name="health_check"),
    path("task_tree/", views.task_tree, name="task_tree"),
    path("active_jobs/", views.active_jobs, name="active_jobs"),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
