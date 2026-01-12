"""
Azure-specific Django settings.

This module extends the core ccp4i2 settings with Azure-specific configuration.
It should be used by setting DJANGO_SETTINGS_MODULE=azure_extensions.settings
in the Azure Container App environment.

This enables:
- The azure_extensions Django app (staged uploads, etc.)
- Azure Storage configuration for SAS URL uploads
- Azure-specific URL routing
"""

import os

# Import all settings from the core ccp4i2 settings
from ccp4i2.config.settings import *  # noqa: F401, F403

# Add azure_extensions and compounds apps to INSTALLED_APPS
INSTALLED_APPS = INSTALLED_APPS + [  # noqa: F405
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

# Platform admin emails (bootstrap admins from environment)
PLATFORM_ADMIN_EMAILS = os.environ.get("PLATFORM_ADMIN_EMAILS", "").split(",")
PLATFORM_ADMIN_EMAILS = [e.strip().lower() for e in PLATFORM_ADMIN_EMAILS if e.strip()]

print("Azure Extensions app enabled")
print("Compounds app enabled (registry, assays, constructs)")
if PLATFORM_ADMIN_EMAILS:
    print(f"Platform admins configured: {len(PLATFORM_ADMIN_EMAILS)} from environment")
