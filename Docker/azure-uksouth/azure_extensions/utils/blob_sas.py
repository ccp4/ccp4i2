"""
Azure Blob Storage SAS URL generation for staged uploads.

This module provides utilities for generating SAS (Shared Access Signature) URLs
that allow clients to upload large files directly to Azure Blob Storage, bypassing
the HTTP body size limits of the API gateway.
"""

import logging
import os
from datetime import datetime, timedelta, timezone
from uuid import uuid4

from django.conf import settings

logger = logging.getLogger(f"azure_extensions:{__name__}")

# Container name for staged uploads
STAGING_CONTAINER = "staging-uploads"

# Default SAS token expiry (1 hour)
DEFAULT_EXPIRY_HOURS = 1


def get_blob_service_client():
    """
    Get Azure Blob Service Client using managed identity or connection string.

    In Azure Container Apps, uses DefaultAzureCredential (managed identity).
    For local development, falls back to connection string from settings.
    """
    try:
        from azure.storage.blob import BlobServiceClient
        from azure.identity import DefaultAzureCredential
    except ImportError:
        logger.error("azure-storage-blob or azure-identity not installed")
        raise ImportError(
            "Azure storage libraries required. Install with: "
            "pip install azure-storage-blob azure-identity"
        )

    # Try to get storage account name from settings
    storage_account_name = getattr(settings, "AZURE_STORAGE_ACCOUNT_NAME", None)

    if storage_account_name:
        # Use managed identity (recommended for Azure deployment)
        account_url = f"https://{storage_account_name}.blob.core.windows.net"
        try:
            # For User-Assigned Managed Identity, we need to specify the client ID
            # This can be set via AZURE_CLIENT_ID env var or passed explicitly
            managed_identity_client_id = os.environ.get("AZURE_CLIENT_ID")
            if managed_identity_client_id:
                logger.debug(f"Using managed identity with client ID: {managed_identity_client_id}")
                credential = DefaultAzureCredential(
                    managed_identity_client_id=managed_identity_client_id
                )
            else:
                logger.debug("Using DefaultAzureCredential without explicit client ID")
                credential = DefaultAzureCredential()
            client = BlobServiceClient(account_url=account_url, credential=credential)
            logger.debug(f"Using managed identity for storage account: {storage_account_name}")
            return client
        except Exception as e:
            logger.warning(f"Managed identity auth failed: {e}, trying connection string")

    # Fallback to connection string (for local development)
    connection_string = getattr(settings, "AZURE_STORAGE_CONNECTION_STRING", None)
    if connection_string:
        client = BlobServiceClient.from_connection_string(connection_string)
        logger.debug("Using connection string for blob storage")
        return client

    raise ValueError(
        "Azure storage not configured. Set AZURE_STORAGE_ACCOUNT_NAME or "
        "AZURE_STORAGE_CONNECTION_STRING in settings."
    )


def generate_upload_sas_url(
    filename: str,
    upload_type: str,
    expiry_hours: int = DEFAULT_EXPIRY_HOURS,
) -> tuple[str, str, datetime]:
    """
    Generate a SAS URL for uploading a file to blob storage.

    Args:
        filename: Original filename (used for extension, sanitized for blob name)
        upload_type: Type of upload ('project_import' or 'unmerged_data')
        expiry_hours: Hours until the SAS token expires

    Returns:
        Tuple of (sas_url, blob_path, expiry_datetime)

    Raises:
        ValueError: If storage is not configured
        Exception: If SAS generation fails
    """
    from azure.storage.blob import (
        BlobSasPermissions,
        generate_blob_sas,
    )

    # Generate unique blob name with original extension
    extension = filename.rsplit(".", 1)[-1] if "." in filename else "zip"
    unique_id = str(uuid4())
    blob_name = f"{upload_type}/{unique_id}.{extension}"

    # Calculate expiry time
    expiry_time = datetime.now(timezone.utc) + timedelta(hours=expiry_hours)

    # Get blob service client
    blob_service_client = get_blob_service_client()

    # Ensure container exists
    container_client = blob_service_client.get_container_client(STAGING_CONTAINER)
    try:
        container_client.get_container_properties()
    except Exception:
        logger.info(f"Creating staging container: {STAGING_CONTAINER}")
        container_client.create_container()

    # Get account details for SAS generation
    account_name = blob_service_client.account_name

    # For managed identity, we need to use user delegation key
    # For connection string/account key, we use account key directly
    try:
        # Try user delegation SAS (for managed identity)
        user_delegation_key = blob_service_client.get_user_delegation_key(
            key_start_time=datetime.now(timezone.utc) - timedelta(minutes=5),
            key_expiry_time=expiry_time,
        )

        sas_token = generate_blob_sas(
            account_name=account_name,
            container_name=STAGING_CONTAINER,
            blob_name=blob_name,
            user_delegation_key=user_delegation_key,
            permission=BlobSasPermissions(create=True, write=True),
            expiry=expiry_time,
        )
        logger.debug("Generated user delegation SAS token")

    except Exception as e:
        logger.debug(f"User delegation SAS failed ({e}), trying account key SAS")

        # Fallback to account key SAS (for connection string auth)
        # This requires the account key, which is available when using connection string
        account_key = getattr(settings, "AZURE_STORAGE_ACCOUNT_KEY", None)

        if not account_key:
            # Try to extract from connection string
            conn_str = getattr(settings, "AZURE_STORAGE_CONNECTION_STRING", "")
            if "AccountKey=" in conn_str:
                account_key = conn_str.split("AccountKey=")[1].split(";")[0]

        if not account_key:
            raise ValueError(
                "Cannot generate SAS: no account key available. "
                "Ensure managed identity has Storage Blob Data Contributor role, "
                "or provide AZURE_STORAGE_ACCOUNT_KEY in settings."
            )

        sas_token = generate_blob_sas(
            account_name=account_name,
            container_name=STAGING_CONTAINER,
            blob_name=blob_name,
            account_key=account_key,
            permission=BlobSasPermissions(create=True, write=True),
            expiry=expiry_time,
        )
        logger.debug("Generated account key SAS token")

    # Construct full URL
    sas_url = f"https://{account_name}.blob.core.windows.net/{STAGING_CONTAINER}/{blob_name}?{sas_token}"

    logger.info(f"Generated SAS URL for upload: {blob_name}, expires: {expiry_time}")

    return sas_url, blob_name, expiry_time


def verify_blob_exists(blob_path: str) -> bool:
    """
    Verify that a blob exists in the staging container.

    Args:
        blob_path: Path to the blob within the staging container

    Returns:
        True if blob exists, False otherwise
    """
    try:
        blob_service_client = get_blob_service_client()
        blob_client = blob_service_client.get_blob_client(
            container=STAGING_CONTAINER,
            blob=blob_path,
        )
        blob_client.get_blob_properties()
        return True
    except Exception as e:
        logger.debug(f"Blob not found or error: {blob_path} - {e}")
        return False


def get_blob_local_path(blob_path: str, download_dir: str = "/tmp/staged-uploads") -> str:
    """
    Download a blob to a local path for processing.

    Args:
        blob_path: Path to the blob within the staging container
        download_dir: Local directory to download to

    Returns:
        Local file path where blob was downloaded
    """
    blob_service_client = get_blob_service_client()
    blob_client = blob_service_client.get_blob_client(
        container=STAGING_CONTAINER,
        blob=blob_path,
    )

    # Create download directory
    os.makedirs(download_dir, exist_ok=True)

    # Extract filename from blob path
    filename = blob_path.split("/")[-1]
    local_path = os.path.join(download_dir, filename)

    logger.info(f"Downloading blob {blob_path} to {local_path}")

    with open(local_path, "wb") as f:
        download_stream = blob_client.download_blob()
        f.write(download_stream.readall())

    logger.info(f"Downloaded {blob_path} ({os.path.getsize(local_path)} bytes)")

    return local_path


def delete_blob(blob_path: str) -> bool:
    """
    Delete a blob from the staging container.

    Args:
        blob_path: Path to the blob within the staging container

    Returns:
        True if deleted successfully, False otherwise
    """
    try:
        blob_service_client = get_blob_service_client()
        blob_client = blob_service_client.get_blob_client(
            container=STAGING_CONTAINER,
            blob=blob_path,
        )
        blob_client.delete_blob()
        logger.info(f"Deleted blob: {blob_path}")
        return True
    except Exception as e:
        logger.warning(f"Failed to delete blob {blob_path}: {e}")
        return False
