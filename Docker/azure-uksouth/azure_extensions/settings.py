"""
Azure-specific Django settings.

This module extends the core ccp4i2 settings with Azure-specific configuration.
It should be used by setting DJANGO_SETTINGS_MODULE=azure_extensions.settings
in the Azure Container App environment.

This enables:
- The azure_extensions Django app (staged uploads, etc.)
- Azure Storage configuration for SAS URL uploads
- Azure Blob Storage for Django file uploads (STORAGES)
- Azure-specific URL routing
"""

import os
from pathlib import Path

# Import all settings from the core ccp4i2 settings
from ccp4i2.config.settings import *  # noqa: F401, F403

# Add azure_extensions and compounds apps to INSTALLED_APPS
INSTALLED_APPS = INSTALLED_APPS + [  # noqa: F405
    "storages",  # django-storages for Azure Blob Storage
    "azure_extensions",
    "users",
    "compounds.registry",
    "compounds.assays",
    "compounds.constructs",
    "reversion",
]

# Enable compounds URLs
COMPOUNDS_ENABLED = True

# Azure Storage configuration for staged uploads (large file uploads via SAS URL)
# These are used by the StagedUpload feature for files > 100MB
AZURE_STORAGE_ACCOUNT_NAME = os.environ.get("AZURE_STORAGE_ACCOUNT_NAME")
AZURE_STORAGE_ACCOUNT_KEY = os.environ.get("AZURE_STORAGE_ACCOUNT_KEY")
AZURE_STORAGE_CONNECTION_STRING = os.environ.get("AZURE_STORAGE_CONNECTION_STRING")

if AZURE_STORAGE_ACCOUNT_NAME:
    print(f"Azure Storage configured: {AZURE_STORAGE_ACCOUNT_NAME}")
elif AZURE_STORAGE_CONNECTION_STRING:
    print("Azure Storage configured via connection string")

# Django file storage configuration using Azure Blob Storage
# Uses django-storages[azure] for file uploads (cheaper than Azure Files for object storage)
# Authentication uses DefaultAzureCredential (Managed Identity in production)
if AZURE_STORAGE_ACCOUNT_NAME:
    # Remove deprecated STATICFILES_STORAGE (Django 4.2+ uses STORAGES instead)
    # This is imported from base settings but conflicts with the STORAGES dict
    globals().pop("STATICFILES_STORAGE", None)

    # Import DefaultAzureCredential for Managed Identity authentication
    # For User-Assigned Managed Identity, AZURE_CLIENT_ID env var must be set
    from azure.identity import DefaultAzureCredential

    STORAGES = {
        "default": {
            "BACKEND": "storages.backends.azure_storage.AzureStorage",
            "OPTIONS": {
                "account_name": AZURE_STORAGE_ACCOUNT_NAME,
                "azure_container": "django-uploads",
                # Use DefaultAzureCredential (Managed Identity) instead of account key
                # DefaultAzureCredential automatically uses AZURE_CLIENT_ID for User-Assigned MI
                "token_credential": DefaultAzureCredential(),
            },
        },
        "staticfiles": {
            "BACKEND": "django.contrib.staticfiles.storage.StaticFilesStorage",
        },
    }
    print("Django file storage configured: Azure Blob Storage (django-uploads container)")
else:
    # Local Docker Compose mode (no Azure Blob Storage)
    # Override MEDIA_ROOT to use a path inside the mounted /mnt/projects volume
    # This ensures uploaded files persist across container restarts
    PROJECTS_DIR = Path(os.environ.get("CCP4I2_PROJECTS_DIR", "/mnt/projects"))
    MEDIA_ROOT = PROJECTS_DIR / "media"  # noqa: F405
    MEDIA_ROOT.mkdir(parents=True, exist_ok=True)
    print(f"Django file storage configured: Local filesystem ({MEDIA_ROOT})")

# Platform admin emails (bootstrap admins from environment)
PLATFORM_ADMIN_EMAILS = os.environ.get("PLATFORM_ADMIN_EMAILS", "").split(",")
PLATFORM_ADMIN_EMAILS = [e.strip().lower() for e in PLATFORM_ADMIN_EMAILS if e.strip()]

print("Azure Extensions app enabled")
print("Compounds app enabled (registry, assays, constructs)")
if PLATFORM_ADMIN_EMAILS:
    print(f"Platform admins configured: {len(PLATFORM_ADMIN_EMAILS)} from environment")

# Logging configuration - log errors to stdout for Azure Container Apps
LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "verbose": {
            "format": "{levelname} {asctime} {module} {message}",
            "style": "{",
        },
    },
    "handlers": {
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "verbose",
        },
    },
    "root": {
        "handlers": ["console"],
        "level": "INFO",
    },
    "loggers": {
        "django": {
            "handlers": ["console"],
            "level": "INFO",
            "propagate": False,
        },
        "django.request": {
            "handlers": ["console"],
            "level": "ERROR",  # Log all request errors with full traceback
            "propagate": False,
        },
        "storages": {
            "handlers": ["console"],
            "level": "DEBUG",  # Debug logging for Azure Storage
            "propagate": False,
        },
        "azure": {
            "handlers": ["console"],
            "level": "DEBUG",  # Debug logging for Azure SDK
            "propagate": False,
        },
    },
}
