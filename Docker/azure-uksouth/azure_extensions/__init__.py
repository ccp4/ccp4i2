"""
Azure Extensions Django App.

This app provides Azure-specific functionality for CCP4i2:
- Staged uploads via SAS URLs for large files (bypassing HTTP body limits)
- Azure Blob Storage integration

This app is only installed in Azure deployments and should not be
included in the core CCP4i2 distribution.
"""

default_app_config = "azure_extensions.apps.AzureExtensionsConfig"
