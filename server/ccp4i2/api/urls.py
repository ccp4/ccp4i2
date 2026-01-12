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

# Core API patterns (will be wrapped under /api/ccp4i2/ prefix)
_api_patterns = [
    path("", include(router.urls)),
    path("health/", views.health_check, name="health_check"),
    path("task_tree/", views.task_tree, name="task_tree"),
    path("active_jobs/", views.active_jobs, name="active_jobs"),
]

# Wrap all patterns under /api/ccp4i2/ for multi-app integration
# Frontend routes are under /ccp4i2/, API routes under /api/ccp4i2/
urlpatterns = [
    path("api/ccp4i2/", include(_api_patterns)),
] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

# Include azure_extensions URLs if the app is installed
# Azure extensions are part of ccp4i2, so they go under /api/ccp4i2/
if "azure_extensions" in settings.INSTALLED_APPS:
    try:
        from azure_extensions.urls import urlpatterns as azure_urls
        _api_patterns = _api_patterns + azure_urls
        # Rebuild urlpatterns with azure extensions included
        urlpatterns = [
            path("api/ccp4i2/", include(_api_patterns)),
        ] + static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
        print("Azure Extensions URLs loaded under /api/ccp4i2/")
    except ImportError as e:
        print(f"Warning: azure_extensions in INSTALLED_APPS but import failed: {e}")

# Include users URLs if the app is installed
if "users" in settings.INSTALLED_APPS:
    try:
        from users.urls import urlpatterns as users_urls
        urlpatterns = urlpatterns + [
            path("api/users/", include((users_urls, "users"))),
        ]
        print("Users URLs loaded at /api/users/")
    except ImportError as e:
        print(f"Warning: users in INSTALLED_APPS but import failed: {e}")

# Include compounds URLs if the app is enabled
# Compounds is a separate app, so it gets its own API namespace /api/compounds/
if getattr(settings, "COMPOUNDS_ENABLED", False):
    try:
        from compounds.urls import urlpatterns as compounds_urls
        urlpatterns = urlpatterns + [
            path("api/compounds/", include((compounds_urls, "compounds"))),
        ]
        print("Compounds URLs loaded at /api/compounds/")
    except ImportError as e:
        print(f"Warning: COMPOUNDS_ENABLED but import failed: {e}")
