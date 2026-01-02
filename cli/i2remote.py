#!/usr/bin/env ccp4-python
"""
i2remote - Remote CLI for CCP4i2

A command-line interface for CCP4i2 that mirrors the local `i2` CLI syntax
but communicates with a remote CCP4i2 server via HTTP API.

Uses ccp4-python for CCP4 data awareness (gemmi, clipper, etc.) to enable
inspection and validation of crystallographic data files.

This enables:
- Testing remote deployments with the same syntax as local
- Automated test suites that work against any CCP4i2 instance
- Scripted job submission and monitoring

Usage:
    i2remote login                        Authenticate with Azure AD (opens browser)
    i2remote login --device-code          Authenticate using device code (for headless)
    i2remote logout                       Clear saved credentials
    i2remote whoami                       Show current authentication status

    i2remote projects [list]              List all projects
    i2remote projects create <name>       Create a new project
    i2remote projects show <project>      Show project details
    i2remote projects tree <project>      Show project job tree

    i2remote jobs <project> [list]        List jobs in a project
    i2remote jobs create <project> <task> Create a new job
    i2remote jobs run <project> <job>     Run an existing job
    i2remote jobs tree <project> <job>    Show job file tree
    i2remote jobs clone <project> <job>   Clone a job
    i2remote jobs kpi <project> <job>     Show job KPIs
    i2remote jobs status <project> <job>  Show job status
    i2remote jobs wait <project> <job>    Wait for job completion
    i2remote jobs validate <project> <job> Validate output files (uses gemmi)
    i2remote jobs set-param <project> <job> <param> <value>  Set job parameter
    i2remote jobs get-param <project> <job> <param>  Get job parameter value
    i2remote jobs upload-param <project> <job> <param> <file> [col]  Upload file param

    i2remote files <project> <job> [list] List files in a job
    i2remote files cat <project> <job> <name>  Display file contents

    i2remote report <project> <job>       Get job report

    i2remote export job <project> <job>   Export a job to zip
    i2remote export project <project>     Export a project to zip
    i2remote import <zipfile>             Import a project from zip (small files)

    i2remote upload-request <file>        Get SAS URL for large file upload
    i2remote upload-complete <id>         Complete upload and trigger processing
    i2remote upload-status <id>           Check status of staged upload
    i2remote upload-list                  List all staged uploads

    i2remote config                       Show current configuration
    i2remote config set <key> <value>     Set configuration value

Environment Variables:
    CCP4I2_API_URL      Base URL for the CCP4i2 API (required)
    CCP4I2_API_TOKEN    Authentication token (if required)
    AZURE_CLIENT_ID     Azure AD client ID (for login command)
    AZURE_TENANT_ID     Azure AD tenant ID (for login command)

Configuration File:
    ~/.ccp4i2remote.json - Persistent configuration

Authentication Methods:
    1. Browser login (recommended for interactive use):
       i2remote login

    2. Device code flow (for SSH/headless environments):
       i2remote login --device-code

    3. Azure CLI (requires az cli installed):
       export CCP4I2_API_TOKEN=$(az account get-access-token \\
         --resource <client-id> --query accessToken -o tsv)

    4. Manual token:
       i2remote config set api_token <your-token>

Examples:
    # Authenticate and set server URL
    i2remote config set api_url https://myserver.azurecontainerapps.io/api/proxy
    i2remote login

    # List projects
    i2remote projects

    # Create and run a job
    i2remote jobs create myproject aimless
    i2remote jobs run myproject 1

    # Wait for completion and check results
    i2remote jobs wait myproject 1
    i2remote jobs kpi myproject 1

    # Set job parameters with file paths
    i2remote jobs set-param myproject 5 XYZIN fullPath=/path/to/file.pdb

    # Set job parameters using fileUse references (files from previous jobs)
    i2remote jobs set-param myproject 5 XYZIN '[-1].XYZOUT[0]'        # Last job's XYZOUT
    i2remote jobs set-param myproject 5 HKLIN 'refmac[-1].HKLOUT'     # Last refmac's HKLOUT
    i2remote jobs set-param myproject 5 XYZIN 'fileUse=[-2].XYZOUT'   # Explicit fileUse=
"""

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urljoin

try:
    import requests
except ImportError:
    print("Error: 'requests' package required. Install with: pip install requests", file=sys.stderr)
    sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# Configuration
# ─────────────────────────────────────────────────────────────────────────────

CONFIG_FILE = Path.home() / ".ccp4i2remote.json"

DEFAULT_CONFIG = {
    "api_url": None,
    "api_token": None,
    "refresh_token": None,
    "token_expires_at": None,
    "azure_client_id": None,
    "azure_tenant_id": None,
    "timeout": 30,
    "verify_ssl": True,
}


def load_config() -> dict:
    """Load configuration from file and environment."""
    config = DEFAULT_CONFIG.copy()

    # Load from config file
    if CONFIG_FILE.exists():
        try:
            with open(CONFIG_FILE) as f:
                file_config = json.load(f)
                config.update(file_config)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load config file: {e}", file=sys.stderr)

    # Override with environment variables
    if os.environ.get("CCP4I2_API_URL"):
        config["api_url"] = os.environ["CCP4I2_API_URL"]
    if os.environ.get("CCP4I2_API_TOKEN"):
        config["api_token"] = os.environ["CCP4I2_API_TOKEN"]
    if os.environ.get("AZURE_CLIENT_ID"):
        config["azure_client_id"] = os.environ["AZURE_CLIENT_ID"]
    if os.environ.get("AZURE_TENANT_ID"):
        config["azure_tenant_id"] = os.environ["AZURE_TENANT_ID"]

    return config


def save_config(config: dict):
    """Save configuration to file."""
    try:
        with open(CONFIG_FILE, "w") as f:
            json.dump(config, f, indent=2)
        os.chmod(CONFIG_FILE, 0o600)  # Secure permissions
    except IOError as e:
        print(f"Error: Could not save config file: {e}", file=sys.stderr)
        sys.exit(1)


CONFIG = load_config()


# ─────────────────────────────────────────────────────────────────────────────
# Azure AD Authentication
# ─────────────────────────────────────────────────────────────────────────────

def get_azure_config() -> tuple:
    """Get Azure AD client and tenant IDs from config or prompt."""
    client_id = CONFIG.get("azure_client_id")
    tenant_id = CONFIG.get("azure_tenant_id")

    if not client_id:
        print("Azure Client ID not configured.")
        print("Set via: i2remote config set azure_client_id <your-client-id>")
        print("Or env:  export AZURE_CLIENT_ID=<your-client-id>")
        sys.exit(1)

    if not tenant_id:
        print("Azure Tenant ID not configured.")
        print("Set via: i2remote config set azure_tenant_id <your-tenant-id>")
        print("Or env:  export AZURE_TENANT_ID=<your-tenant-id>")
        sys.exit(1)

    return client_id, tenant_id


def check_msal_available() -> bool:
    """Check if MSAL library is available."""
    try:
        import msal
        return True
    except ImportError:
        return False


def login_interactive(client_id: str, tenant_id: str) -> dict:
    """
    Perform interactive browser-based login.
    Opens a browser window for the user to authenticate.
    """
    try:
        import msal
    except ImportError:
        print("Error: MSAL library required for login.", file=sys.stderr)
        print("Install with: pip install msal", file=sys.stderr)
        sys.exit(1)

    authority = f"https://login.microsoftonline.com/{tenant_id}"
    # Use .default scope - works without explicit scope configuration in Azure AD
    # This requests all statically configured permissions for the app
    scopes = [f"{client_id}/.default"]

    app = msal.PublicClientApplication(
        client_id,
        authority=authority,
    )

    print("Opening browser for authentication...")
    print("If browser doesn't open, use: i2remote login --device-code")

    try:
        # Try interactive login with browser using localhost redirect
        # MSAL will start a local server to receive the auth code
        result = app.acquire_token_interactive(
            scopes=scopes,
            prompt="select_account",
        )
    except Exception as e:
        print(f"Browser login failed: {e}", file=sys.stderr)
        print("Try using: i2remote login --device-code", file=sys.stderr)
        sys.exit(1)

    if "error" in result:
        print(f"Authentication failed: {result.get('error_description', result.get('error'))}", file=sys.stderr)
        sys.exit(1)

    return result


def login_device_code(client_id: str, tenant_id: str) -> dict:
    """
    Perform device code flow login.
    User visits a URL and enters a code to authenticate.
    Works in headless/SSH environments.
    """
    try:
        import msal
    except ImportError:
        print("Error: MSAL library required for login.", file=sys.stderr)
        print("Install with: pip install msal", file=sys.stderr)
        sys.exit(1)

    authority = f"https://login.microsoftonline.com/{tenant_id}"
    # Use .default scope - works without explicit scope configuration in Azure AD
    scopes = [f"{client_id}/.default"]

    app = msal.PublicClientApplication(
        client_id,
        authority=authority,
    )

    # Initiate device code flow
    flow = app.initiate_device_flow(scopes=scopes)

    if "error" in flow:
        print(f"Error initiating device code flow: {flow.get('error_description', flow.get('error'))}", file=sys.stderr)
        sys.exit(1)

    # Display the message to user
    print()
    print("=" * 60)
    print(flow["message"])
    print("=" * 60)
    print()

    # Wait for user to complete authentication
    result = app.acquire_token_by_device_flow(flow)

    if "error" in result:
        print(f"Authentication failed: {result.get('error_description', result.get('error'))}", file=sys.stderr)
        sys.exit(1)

    return result


def refresh_token_if_needed() -> Optional[str]:
    """
    Check if token is expired and refresh if possible.
    Returns the current valid token or None if refresh failed.
    """
    token = CONFIG.get("api_token")
    refresh_token = CONFIG.get("refresh_token")
    expires_at = CONFIG.get("token_expires_at")

    if not token:
        return None

    # Check if token is expired (with 5 minute buffer)
    if expires_at:
        import datetime
        try:
            expires = datetime.datetime.fromisoformat(expires_at)
            now = datetime.datetime.now(datetime.timezone.utc)
            if now < expires - datetime.timedelta(minutes=5):
                return token  # Token still valid
        except (ValueError, TypeError):
            pass  # Can't parse, try to use token anyway

    # Token expired or unknown, try to refresh
    if refresh_token and CONFIG.get("azure_client_id") and CONFIG.get("azure_tenant_id"):
        try:
            import msal
            client_id = CONFIG["azure_client_id"]
            tenant_id = CONFIG["azure_tenant_id"]
            authority = f"https://login.microsoftonline.com/{tenant_id}"
            scopes = [f"{client_id}/.default"]

            app = msal.PublicClientApplication(client_id, authority=authority)

            # Try to get cached accounts and acquire token silently
            accounts = app.get_accounts()
            if accounts:
                result = app.acquire_token_silent(scopes, account=accounts[0])
                if result and "access_token" in result:
                    save_token_to_config(result)
                    return result["access_token"]
        except Exception:
            pass  # Refresh failed, return existing token

    return token


def save_token_to_config(result: dict):
    """Save authentication result to config file."""
    import datetime

    CONFIG["api_token"] = result.get("access_token")

    if "refresh_token" in result:
        CONFIG["refresh_token"] = result["refresh_token"]

    # Calculate expiration time
    if "expires_in" in result:
        expires_at = datetime.datetime.now(datetime.timezone.utc) + \
                     datetime.timedelta(seconds=result["expires_in"])
        CONFIG["token_expires_at"] = expires_at.isoformat()

    save_config(CONFIG)


def decode_token_claims(token: str) -> Optional[dict]:
    """Decode JWT token claims without verification (for display only)."""
    try:
        import base64
        import json

        # JWT is header.payload.signature
        parts = token.split(".")
        if len(parts) != 3:
            return None

        # Decode payload (add padding if needed)
        payload = parts[1]
        padding = 4 - len(payload) % 4
        if padding != 4:
            payload += "=" * padding

        decoded = base64.urlsafe_b64decode(payload)
        return json.loads(decoded)
    except Exception:
        return None


# ─────────────────────────────────────────────────────────────────────────────
# API Client
# ─────────────────────────────────────────────────────────────────────────────

class APIError(Exception):
    """API request failed."""
    def __init__(self, message: str, status_code: int = None, response: dict = None):
        super().__init__(message)
        self.status_code = status_code
        self.response = response


class CCP4i2Client:
    """HTTP client for CCP4i2 API."""

    def __init__(self, api_url: str, token: str = None, timeout: int = 30, verify_ssl: bool = True):
        self.api_url = api_url.rstrip("/")
        self.token = token
        self.timeout = timeout
        self.verify_ssl = verify_ssl
        self.session = requests.Session()

        # Set up headers
        self.session.headers["Content-Type"] = "application/json"
        if token:
            self.session.headers["Authorization"] = f"Bearer {token}"

    def _url(self, endpoint: str) -> str:
        """Build full URL for endpoint."""
        # Ensure endpoint starts with /
        if not endpoint.startswith("/"):
            endpoint = "/" + endpoint
        # Ensure endpoint ends with / for Django
        if not endpoint.endswith("/"):
            endpoint = endpoint + "/"
        return self.api_url + endpoint

    def _request(self, method: str, endpoint: str, **kwargs) -> Any:
        """Make HTTP request and handle response."""
        url = self._url(endpoint)
        kwargs.setdefault("timeout", self.timeout)
        kwargs.setdefault("verify", self.verify_ssl)

        try:
            response = self.session.request(method, url, **kwargs)
        except requests.RequestException as e:
            raise APIError(f"Request failed: {e}")

        # Handle errors
        if not response.ok:
            try:
                error_data = response.json()
                message = error_data.get("error") or error_data.get("detail") or str(error_data)
            except (json.JSONDecodeError, ValueError):
                message = response.text or f"HTTP {response.status_code}"
            raise APIError(message, response.status_code, error_data if 'error_data' in dir() else None)

        # Return JSON if available
        if response.content:
            try:
                return response.json()
            except (json.JSONDecodeError, ValueError):
                return response.text
        return None

    def get(self, endpoint: str, **params) -> Any:
        """GET request."""
        return self._request("GET", endpoint, params=params)

    def post(self, endpoint: str, data: dict = None, **kwargs) -> Any:
        """POST request."""
        return self._request("POST", endpoint, json=data, **kwargs)

    def delete(self, endpoint: str) -> Any:
        """DELETE request."""
        return self._request("DELETE", endpoint)

    # ─────────────────────────────────────────────────────────────────────────
    # Projects
    # ─────────────────────────────────────────────────────────────────────────

    def list_projects(self) -> list:
        """List all projects."""
        return self.get("/projects")

    def get_project(self, project_id: str) -> dict:
        """Get project details."""
        return self.get(f"/projects/{project_id}")

    def create_project(self, name: str, description: str = "") -> dict:
        """Create a new project."""
        return self.post("/projects", {"name": name, "description": description})

    def get_project_tree(self, project_id: str) -> dict:
        """Get project job tree."""
        return self.get(f"/projects/{project_id}/job_tree")

    def get_project_jobs(self, project_id: str) -> list:
        """Get all jobs in a project."""
        return self.get(f"/projects/{project_id}/jobs")

    def export_project(self, project_id: str) -> dict:
        """Export project to ZIP."""
        return self.post(f"/projects/{project_id}/export")

    def import_project(self, zip_path: str, progress_callback=None) -> dict:
        """
        Import project from ZIP file. Creates a new project from the archive.

        For large files, uses streaming upload with extended timeout.

        Args:
            zip_path: Path to the ZIP file to upload
            progress_callback: Optional callback(bytes_sent, total_bytes) for progress
        """
        file_path = Path(zip_path)
        file_size = file_path.stat().st_size

        # Calculate appropriate timeout based on file size
        # Allow ~1 MB/s minimum upload speed + 60s buffer
        upload_timeout = max(300, int(file_size / (1024 * 1024)) + 60)

        # Use requests-toolbelt for multipart streaming if available
        try:
            from requests_toolbelt import MultipartEncoder, MultipartEncoderMonitor

            # Create streaming encoder
            encoder = MultipartEncoder(
                fields={
                    "files": (file_path.name, open(zip_path, "rb"), "application/zip")
                }
            )

            # Wrap with monitor for progress if callback provided
            if progress_callback:
                monitor = MultipartEncoderMonitor(
                    encoder,
                    lambda m: progress_callback(m.bytes_read, encoder.len)
                )
                data = monitor
            else:
                data = encoder

            headers = {"Content-Type": encoder.content_type}
            url = self._url("/projects/import_project")

            response = self.session.post(
                url,
                data=data,
                headers=headers,
                timeout=upload_timeout,
                verify=self.verify_ssl,
            )

        except ImportError:
            # Fallback to standard upload without progress
            with open(zip_path, "rb") as f:
                files = {"files": (file_path.name, f, "application/zip")}
                url = self._url("/projects/import_project")
                response = self.session.post(
                    url,
                    files=files,
                    timeout=upload_timeout,
                    verify=self.verify_ssl,
                )

        # Handle response
        if not response.ok:
            try:
                error_data = response.json()
                message = error_data.get("error") or error_data.get("detail") or str(error_data)
            except (json.JSONDecodeError, ValueError):
                message = response.text or f"HTTP {response.status_code}"
            raise APIError(message, response.status_code)

        if response.content:
            try:
                return response.json()
            except (json.JSONDecodeError, ValueError):
                return {"status": "ok", "response": response.text}
        return {"status": "ok"}

    # ─────────────────────────────────────────────────────────────────────────
    # Staged Uploads (for large files)
    # ─────────────────────────────────────────────────────────────────────────

    def request_upload(self, filename: str, upload_type: str, target_job_uuid: Optional[str] = None) -> dict:
        """
        Request a SAS URL for uploading a large file.

        Args:
            filename: Name of the file to upload
            upload_type: 'project_import' or 'unmerged_data'
            target_job_uuid: For unmerged_data, the target job UUID

        Returns:
            Dict with upload_id, sas_url, blob_path, expiry, instructions
        """
        data = {
            "filename": filename,
            "upload_type": upload_type,
        }
        if target_job_uuid:
            data["target_job_uuid"] = target_job_uuid

        # Use longer timeout for upload requests (SAS URL generation can be slow)
        return self.post("/uploads/request", data, timeout=120)

    def complete_upload(self, upload_id: str) -> dict:
        """
        Mark an upload as complete and trigger processing.

        Args:
            upload_id: The upload UUID returned from request_upload

        Returns:
            Dict with processing status and message
        """
        # Use longer timeout - server needs to verify large files and start processing
        return self.post(f"/uploads/{upload_id}/complete", timeout=300)

    def get_upload_status(self, upload_id: str) -> dict:
        """Get the status of a staged upload."""
        return self.get(f"/uploads/{upload_id}")

    def list_uploads(self, status: Optional[str] = None, upload_type: Optional[str] = None) -> list:
        """List staged uploads, optionally filtered."""
        params = {}
        if status:
            params["status"] = status
        if upload_type:
            params["type"] = upload_type
        return self.get("/uploads", **params)

    def cancel_upload(self, upload_id: str) -> dict:
        """Cancel a pending upload."""
        return self.delete(f"/uploads/{upload_id}/cancel")

    def reset_upload(self, upload_id: str) -> dict:
        """Reset a stuck upload back to 'uploaded' status."""
        return self.post(f"/uploads/{upload_id}/reset")

    def force_complete_upload(self, upload_id: str) -> dict:
        """Force completion of an upload, bypassing expiry check."""
        return self.post(f"/uploads/{upload_id}/force-complete", timeout=300)

    # ─────────────────────────────────────────────────────────────────────────
    # Jobs
    # ─────────────────────────────────────────────────────────────────────────

    def list_jobs(self, project_id: str = None) -> list:
        """List jobs, optionally filtered by project."""
        if project_id:
            return self.get_project_jobs(project_id)
        return self.get("/jobs")

    def get_job(self, job_id: str) -> dict:
        """Get job details."""
        return self.get(f"/jobs/{job_id}")

    def get_job_by_number(self, project_id: str, job_number: int) -> dict:
        """Find job by project and job number."""
        jobs = self.get_project_jobs(project_id)
        for job in jobs:
            # API returns 'number' as string
            if str(job.get("number", "")) == str(job_number):
                return job
        raise APIError(f"Job {job_number} not found in project {project_id}")

    def create_job(self, project_id: str, task_name: str) -> dict:
        """Create a new job."""
        return self.post(f"/projects/{project_id}/create_task", {"task_name": task_name})

    def run_job(self, job_id: str) -> dict:
        """Run/submit a job."""
        return self.post(f"/jobs/{job_id}/run")

    def clone_job(self, job_id: str) -> dict:
        """Clone a job."""
        return self.post(f"/jobs/{job_id}/clone")

    def get_job_status(self, job_id: str) -> str:
        """Get job status."""
        job = self.get_job(job_id)
        return job.get("status", "Unknown")

    def get_job_container(self, job_id: str) -> dict:
        """Get full job container."""
        return self.get(f"/jobs/{job_id}/container")

    def get_job_files(self, job_id: str) -> list:
        """Get files associated with a job."""
        return self.get(f"/jobs/{job_id}/files")

    def get_job_report(self, job_id: str) -> dict:
        """Get job report XML."""
        return self.get(f"/jobs/{job_id}/report_xml")

    def get_job_kpi(self, job_id: str) -> dict:
        """Get job KPIs (float and char values)."""
        job = self.get_job(job_id)
        # KPIs are typically embedded in job response or need separate call
        kpis = {
            "float_values": job.get("float_values", []),
            "char_values": job.get("char_values", []),
        }
        return kpis

    def export_job(self, job_id: str) -> bytes:
        """Export job to ZIP."""
        return self.get(f"/jobs/{job_id}/export_job")

    def set_job_parameter(self, job_id: str, param_path: str, value: Any) -> dict:
        """Set a job parameter."""
        return self.post(f"/jobs/{job_id}/set_parameter", {
            "object_path": param_path,
            "value": value
        })

    def get_job_parameter(self, job_id: str, param_path: str) -> Any:
        """Get a job parameter value."""
        return self.get(f"/jobs/{job_id}/get_parameter", object_path=param_path)

    def resolve_fileuse(self, project_id: str, fileuse: str) -> dict:
        """
        Resolve a fileUse reference to file metadata.

        Args:
            project_id: Project ID
            fileuse: FileUse string (e.g., "[-1].XYZOUT[0]")

        Returns:
            Dict with file metadata (project, baseName, dbFileId, relPath, fullPath)
        """
        return self.get(f"/projects/{project_id}/resolve_fileuse", fileuse=fileuse)

    def upload_file_param(self, job_id: str, object_path: str, file_path: str,
                          column_selector: str = None) -> dict:
        """
        Upload a file and set it as a job parameter.

        Args:
            job_id: Job ID
            object_path: Parameter path (e.g., "inputData.XYZIN")
            file_path: Local path to the file to upload
            column_selector: Optional MTZ column selector for reflection files

        Returns:
            Dict with upload result and file metadata
        """
        from pathlib import Path
        file_path = Path(file_path)
        if not file_path.exists():
            raise ValueError(f"File not found: {file_path}")

        url = f"{self.base_url}/jobs/{job_id}/upload_file_param/"
        headers = self._get_headers()
        # Remove Content-Type from headers - requests will set it for multipart
        headers.pop("Content-Type", None)

        data = {"objectPath": object_path}
        if column_selector:
            data["column_selector"] = column_selector

        with open(file_path, "rb") as f:
            files = {"file": (file_path.name, f)}
            response = requests.post(url, headers=headers, data=data, files=files)

        if response.status_code >= 400:
            raise APIError(f"Upload failed: {response.status_code} - {response.text}")

        return response.json()

    # ─────────────────────────────────────────────────────────────────────────
    # Files
    # ─────────────────────────────────────────────────────────────────────────

    def get_file(self, file_id: str) -> dict:
        """Get file metadata."""
        return self.get(f"/files/{file_id}")

    def download_file(self, file_id: str) -> bytes:
        """Download file content."""
        url = self._url(f"/files/{file_id}/download")
        response = self.session.get(url, timeout=self.timeout, verify=self.verify_ssl)
        if not response.ok:
            raise APIError(f"Failed to download file: HTTP {response.status_code}")
        return response.content

    def get_file_digest(self, file_id: str) -> dict:
        """Get file digest/summary."""
        return self.get(f"/files/{file_id}/digest")


# ─────────────────────────────────────────────────────────────────────────────
# FileUse Pattern Detection
# ─────────────────────────────────────────────────────────────────────────────

import re

# Regex patterns for fileUse syntax detection
FILEUSE_PATTERNS = [
    # task_name[jobIndex].jobParamName[paramIndex]
    re.compile(r"^(?:\w+)?\[-?\d+\]\.\w+(?:\[\d+\])?$"),
]


def is_fileuse_pattern(value: str) -> bool:
    """
    Check if a string matches the fileUse pattern.

    FileUse patterns look like:
        [-1].XYZOUT[0]
        prosmart_refmac[-1].XYZOUT
        refmac[-2].HKLOUT[0]

    Args:
        value: String to check

    Returns:
        True if the string matches a fileUse pattern
    """
    if not value or not isinstance(value, str):
        return False

    # Quick pre-check: must contain [ and ] and .
    if '[' not in value or ']' not in value or '.' not in value:
        return False

    # Strip fileUse= prefix if present
    if value.startswith("fileUse="):
        value = value[8:]

    # Strip quotes
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    if value.startswith("'") and value.endswith("'"):
        value = value[1:-1]

    # Try the pattern
    for pattern in FILEUSE_PATTERNS:
        if pattern.match(value):
            return True

    return False


def parse_fileuse_value(value: str) -> str:
    """
    Extract the fileUse string from a value.

    Handles both explicit (fileUse=...) and implicit patterns.

    Args:
        value: Raw value string

    Returns:
        Clean fileUse string
    """
    # Strip fileUse= prefix if present
    if value.startswith("fileUse="):
        value = value[8:]

    # Strip quotes
    if value.startswith('"') and value.endswith('"'):
        value = value[1:-1]
    if value.startswith("'") and value.endswith("'"):
        value = value[1:-1]

    return value


# ─────────────────────────────────────────────────────────────────────────────
# Status Codes (from CCP4i2 CJobStatus)
# ─────────────────────────────────────────────────────────────────────────────

JOB_STATUS = {
    0: "Unknown",
    1: "Pending",
    2: "Running",
    3: "Interrupted",
    4: "Failed",
    5: "Unsatisfactory",
    6: "Finished",
    7: "Exported",
}


def status_name(status_code: int) -> str:
    """Convert numeric status code to name."""
    return JOB_STATUS.get(status_code, f"Unknown({status_code})")


# ─────────────────────────────────────────────────────────────────────────────
# File Validation (using CCP4 libraries via ccp4-python)
# ─────────────────────────────────────────────────────────────────────────────

def validate_mtz(content: bytes, filename: str):
    """Validate MTZ file using gemmi."""
    try:
        import gemmi
        import tempfile
        import os

        # Write to temp file for gemmi to read
        with tempfile.NamedTemporaryFile(suffix=".mtz", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            mtz = gemmi.read_mtz_file(temp_path)
            print(f" OK")
            print(f"      Spacegroup: {mtz.spacegroup.hm}")
            print(f"      Cell: {mtz.cell.a:.2f} {mtz.cell.b:.2f} {mtz.cell.c:.2f} "
                  f"{mtz.cell.alpha:.1f} {mtz.cell.beta:.1f} {mtz.cell.gamma:.1f}")
            print(f"      Columns: {len(mtz.columns)}")
            print(f"      Reflections: {mtz.nreflections}")

            # Show resolution
            print(f"      Resolution: {mtz.resolution_low():.2f} - {mtz.resolution_high():.2f} Å")
        finally:
            os.unlink(temp_path)

    except ImportError:
        print(f" (gemmi not available, {len(content)} bytes)")
    except Exception as e:
        print(f" INVALID: {e}")


def validate_pdb(content: bytes, filename: str):
    """Validate PDB file using gemmi."""
    try:
        import gemmi

        text = content.decode("utf-8", errors="replace")
        structure = gemmi.read_pdb_string(text)

        print(f" OK")
        print(f"      Models: {len(structure)}")
        if structure:
            model = structure[0]
            print(f"      Chains: {len(model)}")
            n_atoms = sum(len(res) for chain in model for res in chain)
            print(f"      Atoms: {n_atoms}")

            # Get spacegroup if available
            if structure.spacegroup_hm:
                print(f"      Spacegroup: {structure.spacegroup_hm}")
            if structure.cell.a > 0:
                print(f"      Cell: {structure.cell.a:.2f} {structure.cell.b:.2f} {structure.cell.c:.2f}")

    except ImportError:
        print(f" (gemmi not available, {len(content)} bytes)")
    except Exception as e:
        print(f" INVALID: {e}")


def validate_cif(content: bytes, filename: str):
    """Validate mmCIF file using gemmi."""
    try:
        import gemmi

        text = content.decode("utf-8", errors="replace")
        structure = gemmi.read_structure_from_string(text)

        print(f" OK")
        print(f"      Models: {len(structure)}")
        if structure:
            model = structure[0]
            print(f"      Chains: {len(model)}")
            n_atoms = sum(len(res) for chain in model for res in chain)
            print(f"      Atoms: {n_atoms}")

    except ImportError:
        print(f" (gemmi not available, {len(content)} bytes)")
    except Exception as e:
        print(f" INVALID: {e}")


# ─────────────────────────────────────────────────────────────────────────────
# Output Formatting
# ─────────────────────────────────────────────────────────────────────────────

def print_table(rows: list, headers: list = None):
    """Print data as a formatted table."""
    if not rows:
        print("(no results)")
        return

    # Get all keys if no headers specified
    if headers is None and rows and isinstance(rows[0], dict):
        headers = list(rows[0].keys())

    # Convert dicts to lists
    if rows and isinstance(rows[0], dict):
        rows = [[str(row.get(h, "")) for h in headers] for row in rows]

    # Calculate column widths
    widths = [len(h) for h in headers]
    for row in rows:
        for i, cell in enumerate(row):
            widths[i] = max(widths[i], len(str(cell)))

    # Print header
    header_line = " | ".join(h.ljust(widths[i]) for i, h in enumerate(headers))
    print(header_line)
    print("-" * len(header_line))

    # Print rows
    for row in rows:
        print(" | ".join(str(cell).ljust(widths[i]) for i, cell in enumerate(row)))


def print_json(data: Any):
    """Print data as formatted JSON."""
    print(json.dumps(data, indent=2, default=str))


def print_tree(node: dict, prefix: str = "", is_last: bool = True):
    """Print a tree structure."""
    connector = "└── " if is_last else "├── "
    name = node.get("task_name") or node.get("name") or str(node.get("id", "?"))
    status = node.get("status", "")
    job_num = node.get("job_number", "")

    line = f"{prefix}{connector}[{job_num}] {name}"
    if status:
        line += f" ({status})"
    print(line)

    children = node.get("children", []) or node.get("jobs", [])
    for i, child in enumerate(children):
        extension = "    " if is_last else "│   "
        print_tree(child, prefix + extension, i == len(children) - 1)


# ─────────────────────────────────────────────────────────────────────────────
# Commands
# ─────────────────────────────────────────────────────────────────────────────

def get_client() -> CCP4i2Client:
    """Get configured API client."""
    if not CONFIG.get("api_url"):
        print("Error: API URL not configured.", file=sys.stderr)
        print("Set CCP4I2_API_URL environment variable or run:", file=sys.stderr)
        print("  i2remote config set api_url <url>", file=sys.stderr)
        sys.exit(1)

    return CCP4i2Client(
        api_url=CONFIG["api_url"],
        token=CONFIG.get("api_token"),
        timeout=CONFIG.get("timeout", 30),
        verify_ssl=CONFIG.get("verify_ssl", True),
    )


def resolve_project(client: CCP4i2Client, identifier: str) -> str:
    """Resolve project identifier to ID (handles name or ID)."""
    # Try as-is first (might be UUID or numeric ID)
    try:
        project = client.get_project(identifier)
        return str(project.get("id") or project.get("uuid") or identifier)
    except APIError:
        pass

    # Try to find by name
    projects = client.list_projects()
    for p in projects:
        if p.get("name", "").lower() == identifier.lower():
            return str(p.get("id") or p.get("uuid"))

    raise APIError(f"Project not found: {identifier}")


def resolve_job(client: CCP4i2Client, project_id: str, job_identifier: str) -> str:
    """Resolve job identifier to ID (handles job number or UUID)."""
    # Try as job number first
    try:
        job_number = int(job_identifier)
        job = client.get_job_by_number(project_id, job_number)
        return str(job.get("id") or job.get("uuid"))
    except ValueError:
        pass
    except APIError:
        pass

    # Try as UUID/ID
    try:
        job = client.get_job(job_identifier)
        return str(job.get("id") or job.get("uuid"))
    except APIError:
        pass

    raise APIError(f"Job not found: {job_identifier}")


# ─────────────────────────────────────────────────────────────────────────────
# Command Handlers
# ─────────────────────────────────────────────────────────────────────────────

def cmd_login(args):
    """Handle login command."""
    client_id, tenant_id = get_azure_config()

    if args.device_code:
        print("Using device code flow (for headless environments)...")
        result = login_device_code(client_id, tenant_id)
    else:
        print("Using interactive browser login...")
        result = login_interactive(client_id, tenant_id)

    # Save tokens
    save_token_to_config(result)

    # Show success
    claims = decode_token_claims(result.get("access_token", ""))
    if claims:
        name = claims.get("name") or claims.get("preferred_username") or claims.get("email") or "Unknown"
        print(f"\n✓ Logged in as: {name}")
    else:
        print("\n✓ Login successful!")

    print(f"Token saved to: {CONFIG_FILE}")


def cmd_logout(args):
    """Handle logout command."""
    # Clear auth-related config
    CONFIG["api_token"] = None
    CONFIG["refresh_token"] = None
    CONFIG["token_expires_at"] = None
    save_config(CONFIG)
    print("Logged out. Credentials cleared.")


def cmd_whoami(args):
    """Handle whoami command - show current authentication status."""
    token = CONFIG.get("api_token")

    if not token:
        print("Not logged in.")
        print("\nTo authenticate, run:")
        print("  i2remote login              # Browser-based login")
        print("  i2remote login --device-code  # For SSH/headless")
        return

    claims = decode_token_claims(token)

    if claims:
        print("Current authentication:")
        print(f"  Name:     {claims.get('name', 'N/A')}")
        print(f"  Email:    {claims.get('preferred_username') or claims.get('email', 'N/A')}")
        print(f"  Subject:  {claims.get('sub', 'N/A')[:20]}...")

        # Check expiration
        import datetime
        exp = claims.get("exp")
        if exp:
            expires = datetime.datetime.fromtimestamp(exp, tz=datetime.timezone.utc)
            now = datetime.datetime.now(datetime.timezone.utc)
            if now > expires:
                print(f"  Status:   EXPIRED (at {expires.isoformat()})")
            else:
                remaining = expires - now
                print(f"  Expires:  {expires.isoformat()} (in {remaining})")
    else:
        print("Token present but could not decode claims.")
        print("Token may be valid - try running a command to test.")

    print(f"\nConfig file: {CONFIG_FILE}")


def cmd_config(args):
    """Handle config command."""
    if args.action == "show" or args.action is None:
        print("Current configuration:")
        for key, value in CONFIG.items():
            # Hide sensitive values
            if key in ("api_token", "refresh_token") and value:
                print(f"  {key}: ****")
            else:
                print(f"  {key}: {value}")
        print(f"\nConfig file: {CONFIG_FILE}")
    elif args.action == "set":
        if not args.key or not args.value:
            print("Usage: i2remote config set <key> <value>", file=sys.stderr)
            sys.exit(1)
        CONFIG[args.key] = args.value
        save_config(CONFIG)
        print(f"Set {args.key} = {args.value}")


def cmd_projects(args):
    """Handle projects command."""
    client = get_client()

    action = args.action or "list"

    if action == "list":
        projects = client.list_projects()
        print_table(projects, ["id", "name", "description", "last_access"])

    elif action == "create":
        if not args.name:
            print("Usage: i2remote projects create <name>", file=sys.stderr)
            sys.exit(1)
        project = client.create_project(args.name, args.description or "")
        print(f"Created project: {project.get('name')} (ID: {project.get('id')})")

    elif action == "show":
        if not args.project:
            print("Usage: i2remote projects show <project>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        project = client.get_project(project_id)
        print_json(project)

    elif action == "tree":
        if not args.project:
            print("Usage: i2remote projects tree <project>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        tree = client.get_project_tree(project_id)
        print(f"Project: {args.project}")
        if isinstance(tree, list):
            for node in tree:
                print_tree(node)
        else:
            print_tree(tree)

    else:
        # Assume it's a project ID for show
        project_id = resolve_project(client, action)
        project = client.get_project(project_id)
        print_json(project)


def cmd_jobs(args):
    """Handle jobs command."""
    client = get_client()

    action = args.action

    # Handle case where action is actually a project ID (i2remote jobs <project>)
    if action and action not in ("list", "create", "run", "clone", "status", "wait", "kpi", "tree", "validate", "set-param", "get-param", "upload-param"):
        # action is actually the project, shift args
        args.project = action
        action = "list"

    if action == "list" or action is None:
        if not args.project:
            print("Usage: i2remote jobs <project> [list]", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        jobs = client.list_jobs(project_id)
        # Transform jobs for display
        for job in jobs:
            job["status_name"] = status_name(job.get("status", 0))
            job["job_number"] = job.get("number", "")
        print_table(jobs, ["job_number", "task_name", "status_name", "title"])

    elif action == "create":
        if not args.project or not args.task:
            print("Usage: i2remote jobs create <project> <task>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job = client.create_job(project_id, args.task)
        print(f"Created job: {job.get('task_name')} (Job #{job.get('job_number')}, ID: {job.get('id')})")

    elif action == "run":
        if not args.project or not args.job:
            print("Usage: i2remote jobs run <project> <job>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        result = client.run_job(job_id)
        print(f"Job submitted: {result}")

    elif action == "clone":
        if not args.project or not args.job:
            print("Usage: i2remote jobs clone <project> <job>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        cloned = client.clone_job(job_id)
        print(f"Cloned job: {cloned.get('task_name')} (Job #{cloned.get('job_number')}, ID: {cloned.get('id')})")

    elif action == "status":
        if not args.project or not args.job:
            print("Usage: i2remote jobs status <project> <job>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        job = client.get_job(job_id)
        print(f"Job #{job.get('number')} - {job.get('task_name')}")
        print(f"  Title: {job.get('title')}")
        print(f"  Status: {status_name(job.get('status', 0))}")
        print(f"  Created: {job.get('creation_time')}")
        print(f"  Finished: {job.get('finish_time') or 'N/A'}")

    elif action == "wait":
        if not args.project or not args.job:
            print("Usage: i2remote jobs wait <project> <job>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)

        print(f"Waiting for job to complete...", end="", flush=True)
        while True:
            job = client.get_job(job_id)
            status = job.get("status", 0)
            status_str = status_name(status)
            if status >= 4:  # Failed, Unsatisfactory, Finished, or Exported
                print(f"\nJob completed with status: {status_str}")
                break
            print(".", end="", flush=True)
            time.sleep(5)

    elif action == "validate":
        # Validate job output files using CCP4 libraries
        if not args.project or not args.job:
            print("Usage: i2remote jobs validate <project> <job>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)

        job = client.get_job(job_id)
        print(f"Validating outputs for Job #{job.get('number')} - {job.get('task_name')}")

        files = client.get_job_files(job_id)
        for f in files:
            filetype = f.get("filetype", "").lower()
            filename = f.get("filename", "")
            role = f.get("role", "")

            # Only validate output files
            if role != "output":
                continue

            print(f"\n  {filename} ({filetype}):", end="")

            try:
                content = client.download_file(f.get("id"))

                # Validate based on file type using CCP4 libraries
                if filetype in ("mtz", "mini_mtz", "observed_mtz", "calculated_mtz"):
                    validate_mtz(content, filename)
                elif filetype in ("pdb", "xyz", "model"):
                    validate_pdb(content, filename)
                elif filetype in ("cif", "mmcif"):
                    validate_cif(content, filename)
                else:
                    print(f" {len(content)} bytes (no validator)")
            except Exception as e:
                print(f" ERROR: {e}")

    elif action == "kpi":
        if not args.project or not args.job:
            print("Usage: i2remote jobs kpi <project> <job>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        job = client.get_job(job_id)

        print(f"KPIs for Job #{job.get('job_number')} - {job.get('task_name')}")
        print("\nFloat Values:")
        for kpi in job.get("float_values", []):
            print(f"  {kpi.get('name')}: {kpi.get('value')}")
        print("\nChar Values:")
        for kpi in job.get("char_values", []):
            print(f"  {kpi.get('name')}: {kpi.get('value')}")

    elif action == "tree":
        if not args.project or not args.job:
            print("Usage: i2remote jobs tree <project> <job>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        files = client.get_job_files(job_id)
        print(f"Files for job {args.job}:")
        print_table(files, ["id", "filename", "filetype", "role"])

    elif action == "set-param":
        # Set a job parameter, with support for fileUse references
        if not args.project or not args.job or not args.param_name:
            print("Usage: i2remote jobs set-param <project> <job> <param_name> <value>...", file=sys.stderr)
            print("\nExamples:", file=sys.stderr)
            print("  i2remote jobs set-param myproject 5 XYZIN fullPath=/path/to/file.pdb", file=sys.stderr)
            print("  i2remote jobs set-param myproject 5 XYZIN '[-1].XYZOUT[0]'", file=sys.stderr)
            print("  i2remote jobs set-param myproject 5 HKLIN 'prosmart_refmac[-1].HKLOUT'", file=sys.stderr)
            sys.exit(1)

        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        param_name = args.param_name
        param_values = args.param_values or []

        if not param_values:
            print("Error: No value provided for parameter", file=sys.stderr)
            sys.exit(1)

        # Join values for processing
        value = " ".join(param_values) if len(param_values) > 1 else param_values[0]

        # Check if this is a fileUse reference (explicit or pattern-matched)
        is_fileuse = value.startswith("fileUse=") or is_fileuse_pattern(value)

        if is_fileuse:
            # Resolve the fileUse reference
            fileuse_str = parse_fileuse_value(value)
            print(f"Resolving fileUse reference: {fileuse_str}")

            try:
                result = client.resolve_fileuse(project_id, fileuse_str)
                if result.get("status") == "Success":
                    file_data = result.get("data", {})
                    print(f"  Resolved to: {file_data.get('baseName')} ({file_data.get('fullPath')})")

                    # Build the value for set_parameter using fullPath
                    value = f"fullPath={file_data.get('fullPath')}"
                else:
                    print(f"Error resolving fileUse: {result.get('error', 'Unknown error')}", file=sys.stderr)
                    sys.exit(1)
            except APIError as e:
                print(f"Error resolving fileUse: {e}", file=sys.stderr)
                sys.exit(1)

        # Set the parameter
        # The param_name should be prefixed with inputData. for most file params
        # But let the user specify the full path if needed
        if "." not in param_name:
            object_path = f"inputData.{param_name}"
        else:
            object_path = param_name

        print(f"Setting {object_path} = {value}")
        result = client.set_job_parameter(job_id, object_path, value)

        if result.get("status") == "Success":
            print(f"Parameter set successfully")
            data = result.get("data", {})
            if data.get("message"):
                print(f"  {data.get('message')}")
        else:
            print(f"Error: {result.get('error', 'Unknown error')}", file=sys.stderr)
            sys.exit(1)

    elif action == "get-param":
        # Get a job parameter value
        if not args.project or not args.job or not args.param_name:
            print("Usage: i2remote jobs get-param <project> <job> <param_name>", file=sys.stderr)
            sys.exit(1)

        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        param_name = args.param_name

        # Prefix with inputData. if no dot in path
        if "." not in param_name:
            object_path = f"inputData.{param_name}"
        else:
            object_path = param_name

        result = client.get_job_parameter(job_id, object_path)

        if result.get("status") == "Success":
            data = result.get("data", {})
            print(f"{object_path}:")
            if data.get("file_path"):
                print(f"  File: {data.get('file_path')}")
                print(f"  Base name: {data.get('base_name')}")
                if data.get("db_file_id"):
                    print(f"  DB ID: {data.get('db_file_id')}")
            elif data.get("value") is not None:
                print(f"  Value: {data.get('value')}")
            else:
                print_json(data)
        else:
            print(f"Error: {result.get('error', 'Unknown error')}", file=sys.stderr)
            sys.exit(1)

    elif action == "upload-param":
        # Upload a file and set it as a job parameter
        if not args.project or not args.job or not args.param_name or not args.param_values:
            print("Usage: i2remote jobs upload-param <project> <job> <param_name> <file_path> [column_selector]", file=sys.stderr)
            print("\nExamples:", file=sys.stderr)
            print("  i2remote jobs upload-param myproject 5 XYZIN /path/to/model.pdb", file=sys.stderr)
            print("  i2remote jobs upload-param myproject 5 F_SIGF /path/to/data.mtz '/*/*/[FP,SIGFP]'", file=sys.stderr)
            sys.exit(1)

        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        param_name = args.param_name
        file_path = args.param_values[0]
        column_selector = args.param_values[1] if len(args.param_values) > 1 else None

        # Build object_path (prefix with inputData. if no dot in path)
        if "." not in param_name:
            object_path = f"inputData.{param_name}"
        else:
            object_path = param_name

        # Check file exists
        if not Path(file_path).exists():
            print(f"Error: File not found: {file_path}", file=sys.stderr)
            sys.exit(1)

        print(f"Uploading {file_path} to {object_path}...")
        if column_selector:
            print(f"  Column selector: {column_selector}")

        try:
            result = client.upload_file_param(job_id, object_path, file_path, column_selector)
            if result.get("status") == "Success":
                data = result.get("data", {})
                updated = data.get("updated_item", {})
                print(f"Upload successful:")
                print(f"  Base name: {updated.get('baseName', 'N/A')}")
                print(f"  Annotation: {updated.get('annotation', 'N/A')}")
                if updated.get('dbFileId'):
                    print(f"  DB ID: {updated.get('dbFileId')}")
            else:
                print(f"Error: {result.get('error', 'Unknown error')}", file=sys.stderr)
                sys.exit(1)
        except APIError as e:
            print(f"Error uploading file: {e}", file=sys.stderr)
            sys.exit(1)
        except ValueError as e:
            print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)

    else:
        print(f"Unknown action: {action}", file=sys.stderr)
        sys.exit(1)


def cmd_files(args):
    """Handle files command."""
    client = get_client()

    action = args.action or "list"

    if action == "list":
        if not args.project or not args.job:
            print("Usage: i2remote files <project> <job> [list]", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        files = client.get_job_files(job_id)
        print_table(files, ["id", "filename", "filetype", "role"])

    elif action == "cat":
        if not args.project or not args.job or not args.filename:
            print("Usage: i2remote files cat <project> <job> <filename>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        files = client.get_job_files(job_id)

        # Find file by name
        target_file = None
        for f in files:
            if f.get("filename") == args.filename or str(f.get("id")) == args.filename:
                target_file = f
                break

        if not target_file:
            print(f"File not found: {args.filename}", file=sys.stderr)
            sys.exit(1)

        content = client.download_file(target_file.get("id"))
        print(content.decode("utf-8", errors="replace"))


def cmd_report(args):
    """Handle report command."""
    client = get_client()

    if not args.project or not args.job:
        print("Usage: i2remote report <project> <job>", file=sys.stderr)
        sys.exit(1)

    project_id = resolve_project(client, args.project)
    job_id = resolve_job(client, project_id, args.job)
    report = client.get_job_report(job_id)

    if isinstance(report, dict) and "xml" in report:
        print(report["xml"])
    else:
        print_json(report)


def cmd_export(args):
    """Handle export command."""
    client = get_client()

    if args.type == "job":
        if not args.project or not args.job:
            print("Usage: i2remote export job <project> <job>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        job_id = resolve_job(client, project_id, args.job)
        result = client.export_job(job_id)
        print(f"Export: {result}")

    elif args.type == "project":
        if not args.project:
            print("Usage: i2remote export project <project>", file=sys.stderr)
            sys.exit(1)
        project_id = resolve_project(client, args.project)
        result = client.export_project(project_id)
        print(f"Export: {result}")


def cmd_import(args):
    """Handle import command."""
    client = get_client()

    if not args.zipfile:
        print("Usage: i2remote import <zipfile>", file=sys.stderr)
        sys.exit(1)

    zip_path = Path(args.zipfile)
    if not zip_path.exists():
        print(f"Error: File not found: {args.zipfile}", file=sys.stderr)
        sys.exit(1)

    if not zip_path.suffix == ".zip":
        print(f"Error: File must be a .zip archive: {args.zipfile}", file=sys.stderr)
        sys.exit(1)

    # Get file size for progress indication
    file_size = zip_path.stat().st_size
    file_size_mb = file_size / (1024 * 1024)
    file_size_gb = file_size_mb / 1024

    # Warn about large files - Azure Container Apps has request size limits
    LARGE_FILE_THRESHOLD_MB = 500  # 500 MB
    if file_size_mb > LARGE_FILE_THRESHOLD_MB:
        print(f"WARNING: Large file detected ({file_size_gb:.1f} GB)", file=sys.stderr)
        print("Azure Container Apps has request size limits (~100MB default).", file=sys.stderr)
        print("For large files, consider uploading directly to Azure Files:", file=sys.stderr)
        print("", file=sys.stderr)
        print("  # Upload to Azure Files share", file=sys.stderr)
        print("  az storage file upload \\", file=sys.stderr)
        print("    --account-name <storage-account> \\", file=sys.stderr)
        print("    --share-name ccp4data \\", file=sys.stderr)
        print(f"    --source '{zip_path}' \\", file=sys.stderr)
        print(f"    --path 'imports/{zip_path.name}'", file=sys.stderr)
        print("", file=sys.stderr)
        print("  # Then run import from server (via maintenance job)", file=sys.stderr)
        print("", file=sys.stderr)
        response = input("Continue with HTTP upload anyway? [y/N]: ")
        if response.lower() != 'y':
            print("Aborted.")
            sys.exit(0)

    # Calculate expected timeout
    upload_timeout = max(300, int(file_size_mb) + 60)

    print(f"Importing project from: {zip_path.name} ({file_size_mb:.1f} MB)")
    print(f"Upload timeout: {upload_timeout}s (adjust if needed)")

    # Progress tracking
    last_percent = [0]
    start_time = [time.time()]

    def progress_callback(bytes_sent, total_bytes):
        percent = int(100 * bytes_sent / total_bytes)
        if percent > last_percent[0]:
            last_percent[0] = percent
            elapsed = time.time() - start_time[0]
            speed_mbps = (bytes_sent / (1024 * 1024)) / elapsed if elapsed > 0 else 0
            remaining_mb = (total_bytes - bytes_sent) / (1024 * 1024)
            eta = remaining_mb / speed_mbps if speed_mbps > 0 else 0
            print(f"\rUploading: {percent}% ({bytes_sent/(1024*1024):.1f}/{total_bytes/(1024*1024):.1f} MB) "
                  f"- {speed_mbps:.1f} MB/s - ETA: {int(eta)}s   ", end="", flush=True)

    print("Uploading to server...")

    try:
        result = client.import_project(str(zip_path), progress_callback=progress_callback)
        print()  # New line after progress
        if result.get("imported"):
            print("Project imported successfully!")
            print("Note: Import runs in background - project will appear shortly")
        else:
            print(f"Import response: {result}")
    except APIError as e:
        print()  # New line after progress
        print(f"Import failed: {e}", file=sys.stderr)
        sys.exit(1)
    except requests.exceptions.SSLError as e:
        print()  # New line after progress
        print(f"SSL/Connection error: {e}", file=sys.stderr)
        print("", file=sys.stderr)
        print("This usually means the file is too large for HTTP upload.", file=sys.stderr)
        print("Use 'i2remote upload-request' for large files instead.", file=sys.stderr)
        sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# Staged Upload Commands (for large files)
# ─────────────────────────────────────────────────────────────────────────────

def cmd_upload_request(args):
    """Request a SAS URL for uploading a large file."""
    client = get_client()

    if not args.file:
        print("Usage: i2remote upload-request <file> [--type project_import|unmerged_data]", file=sys.stderr)
        sys.exit(1)

    file_path = Path(args.file)
    if not file_path.exists():
        print(f"Error: File not found: {args.file}", file=sys.stderr)
        sys.exit(1)

    upload_type = args.type or "project_import"
    if upload_type not in ["project_import", "unmerged_data"]:
        print(f"Error: Invalid upload type: {upload_type}", file=sys.stderr)
        print("Valid types: project_import, unmerged_data", file=sys.stderr)
        sys.exit(1)

    file_size = file_path.stat().st_size
    file_size_mb = file_size / (1024 * 1024)

    print(f"Requesting upload URL for: {file_path.name} ({file_size_mb:.1f} MB)")
    print(f"Upload type: {upload_type}")

    try:
        result = client.request_upload(
            filename=file_path.name,
            upload_type=upload_type,
            target_job_uuid=args.job if hasattr(args, 'job') else None,
        )

        print()
        print("=" * 70)
        print("Upload URL generated successfully!")
        print("=" * 70)
        print()
        print(f"Upload ID: {result['upload_id']}")
        print(f"Expires:   {result['expiry']}")
        print()
        print("To upload, run one of these commands:")
        print()
        print("  # Using azcopy (recommended for large files):")
        print(f"  azcopy copy '{file_path}' '{result['sas_url']}'")
        print()
        print("  # Using curl:")
        print("  curl -X PUT -H 'x-ms-blob-type: BlockBlob' \\")
        print(f"       --data-binary @'{file_path}' \\")
        print(f"       '{result['sas_url']}'")
        print()
        print("After upload completes, run:")
        print(f"  i2remote upload-complete {result['upload_id']}")
        print()

    except APIError as e:
        print(f"Failed to request upload URL: {e}", file=sys.stderr)
        sys.exit(1)


def cmd_upload_complete(args):
    """Complete a staged upload and trigger processing."""
    client = get_client()

    if not args.upload_id:
        print("Usage: i2remote upload-complete <upload_id>", file=sys.stderr)
        sys.exit(1)

    print(f"Completing upload: {args.upload_id}")
    print("Verifying file and triggering processing...")

    try:
        result = client.complete_upload(args.upload_id)

        print()
        print("Upload completed successfully!")
        print(f"Status: {result.get('status', 'unknown')}")
        if result.get('message'):
            print(f"Message: {result['message']}")
        if result.get('note'):
            print(f"Note: {result['note']}")

    except APIError as e:
        print(f"Failed to complete upload: {e}", file=sys.stderr)
        sys.exit(1)


def cmd_upload_reset(args):
    """Reset a stuck upload back to 'uploaded' status."""
    client = get_client()

    if not args.upload_id:
        print("Usage: i2remote upload-reset <upload_id>", file=sys.stderr)
        sys.exit(1)

    print(f"Resetting upload: {args.upload_id}")

    try:
        result = client.reset_upload(args.upload_id)

        print()
        print("Upload reset successfully!")
        print(f"Status: {result.get('status', 'unknown')}")
        print()
        print("You can now run 'upload-complete' to re-trigger processing.")

    except APIError as e:
        print(f"Failed to reset upload: {e}", file=sys.stderr)
        sys.exit(1)


def cmd_upload_force_complete(args):
    """Force complete an expired upload."""
    client = get_client()

    if not args.upload_id:
        print("Usage: i2remote upload-force-complete <upload_id>", file=sys.stderr)
        sys.exit(1)

    print(f"Force completing upload: {args.upload_id}")
    print("Verifying blob exists and triggering processing...")

    try:
        result = client.force_complete_upload(args.upload_id)

        print()
        print("Upload force-completed successfully!")
        print(f"Status: {result.get('status', 'unknown')}")
        if result.get('message'):
            print(f"Message: {result['message']}")
        if result.get('note'):
            print(f"Note: {result['note']}")

    except APIError as e:
        print(f"Failed to force-complete upload: {e}", file=sys.stderr)
        sys.exit(1)


def cmd_upload_status(args):
    """Check status of a staged upload."""
    client = get_client()

    if not args.upload_id:
        print("Usage: i2remote upload-status <upload_id>", file=sys.stderr)
        sys.exit(1)

    try:
        result = client.get_upload_status(args.upload_id)

        print(f"Upload ID:     {result.get('uuid', 'N/A')}")
        print(f"Type:          {result.get('upload_type', 'N/A')}")
        print(f"Status:        {result.get('status', 'N/A')}")
        print(f"Filename:      {result.get('original_filename', 'N/A')}")
        print(f"Created:       {result.get('created_at', 'N/A')}")
        print(f"SAS Expiry:    {result.get('sas_expiry', 'N/A')}")
        print(f"Expired:       {result.get('is_expired', 'N/A')}")
        if result.get('completed_at'):
            print(f"Completed:     {result['completed_at']}")
        if result.get('error_message'):
            print(f"Error:         {result['error_message']}")

    except APIError as e:
        print(f"Failed to get upload status: {e}", file=sys.stderr)
        sys.exit(1)


def cmd_upload_list(args):
    """List staged uploads."""
    client = get_client()

    try:
        uploads = client.list_uploads(
            status=args.status if hasattr(args, 'status') else None,
            upload_type=args.type if hasattr(args, 'type') else None,
        )

        if not uploads:
            print("No staged uploads found.")
            return

        print(f"{'ID':<36}  {'Type':<16}  {'Status':<12}  {'Filename'}")
        print("-" * 100)
        for upload in uploads:
            print(f"{upload.get('uuid', 'N/A'):<36}  "
                  f"{upload.get('upload_type', 'N/A'):<16}  "
                  f"{upload.get('status', 'N/A'):<12}  "
                  f"{upload.get('original_filename', 'N/A')}")

    except APIError as e:
        print(f"Failed to list uploads: {e}", file=sys.stderr)
        sys.exit(1)


# ─────────────────────────────────────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Remote CLI for CCP4i2",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )

    subparsers = parser.add_subparsers(dest="command", help="Command")

    # Login command
    login_parser = subparsers.add_parser("login", help="Authenticate with Azure AD")
    login_parser.add_argument("--device-code", action="store_true",
                              help="Use device code flow (for SSH/headless environments)")

    # Logout command
    subparsers.add_parser("logout", help="Clear saved credentials")

    # Whoami command
    subparsers.add_parser("whoami", help="Show current authentication status")

    # Config command
    config_parser = subparsers.add_parser("config", help="Configuration")
    config_parser.add_argument("action", nargs="?", choices=["show", "set"])
    config_parser.add_argument("key", nargs="?")
    config_parser.add_argument("value", nargs="?")

    # Projects command
    projects_parser = subparsers.add_parser("projects", help="Project operations")
    projects_parser.add_argument("action", nargs="?", help="list|create|show|tree or project ID")
    projects_parser.add_argument("project", nargs="?", help="Project name/ID")
    projects_parser.add_argument("--name", "-n", help="Project name for create")
    projects_parser.add_argument("--description", "-d", help="Project description")

    # Jobs command
    jobs_parser = subparsers.add_parser("jobs", help="Job operations")
    jobs_parser.add_argument("action", nargs="?", help="list|create|run|clone|status|wait|kpi|tree|set-param|get-param|upload-param")
    jobs_parser.add_argument("project", nargs="?", help="Project name/ID")
    jobs_parser.add_argument("job", nargs="?", help="Job number or ID")
    jobs_parser.add_argument("param_name", nargs="?", help="Parameter name for set-param/get-param")
    jobs_parser.add_argument("param_values", nargs="*", help="Parameter value(s) for set-param")
    jobs_parser.add_argument("--task", "-t", help="Task name for create")

    # Files command
    files_parser = subparsers.add_parser("files", help="File operations")
    files_parser.add_argument("action", nargs="?", help="list|cat")
    files_parser.add_argument("project", nargs="?", help="Project name/ID")
    files_parser.add_argument("job", nargs="?", help="Job number or ID")
    files_parser.add_argument("filename", nargs="?", help="Filename for cat")

    # Report command
    report_parser = subparsers.add_parser("report", help="Get job report")
    report_parser.add_argument("project", help="Project name/ID")
    report_parser.add_argument("job", help="Job number or ID")

    # Export command
    export_parser = subparsers.add_parser("export", help="Export job or project")
    export_parser.add_argument("type", choices=["job", "project"])
    export_parser.add_argument("project", help="Project name/ID")
    export_parser.add_argument("job", nargs="?", help="Job number or ID (for job export)")

    # Import command
    import_parser = subparsers.add_parser("import", help="Import project from ZIP")
    import_parser.add_argument("zipfile", help="Path to ZIP file")

    # Upload commands (for large files that exceed HTTP body limits)
    upload_request_parser = subparsers.add_parser(
        "upload-request",
        help="Request a SAS URL for uploading large files"
    )
    upload_request_parser.add_argument("file", help="Path to file to upload")
    upload_request_parser.add_argument(
        "--type", "-t",
        choices=["project_import", "unmerged_data"],
        default="project_import",
        help="Type of upload (default: project_import)"
    )
    upload_request_parser.add_argument(
        "--job", "-j",
        help="Target job UUID (for unmerged_data uploads)"
    )

    upload_complete_parser = subparsers.add_parser(
        "upload-complete",
        help="Complete a staged upload and trigger processing"
    )
    upload_complete_parser.add_argument("upload_id", help="Upload UUID from upload-request")

    upload_status_parser = subparsers.add_parser(
        "upload-status",
        help="Check status of a staged upload"
    )
    upload_status_parser.add_argument("upload_id", help="Upload UUID")

    upload_list_parser = subparsers.add_parser(
        "upload-list",
        help="List staged uploads"
    )
    upload_list_parser.add_argument(
        "--status", "-s",
        choices=["pending", "uploaded", "processing", "completed", "failed", "expired"],
        help="Filter by status"
    )
    upload_list_parser.add_argument(
        "--type", "-t",
        choices=["project_import", "unmerged_data"],
        help="Filter by upload type"
    )

    upload_reset_parser = subparsers.add_parser(
        "upload-reset",
        help="Reset a stuck upload to allow re-processing"
    )
    upload_reset_parser.add_argument("upload_id", help="Upload UUID to reset")

    upload_force_complete_parser = subparsers.add_parser(
        "upload-force-complete",
        help="Force complete an expired upload (blob must still exist)"
    )
    upload_force_complete_parser.add_argument("upload_id", help="Upload UUID to force complete")

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    try:
        if args.command == "login":
            cmd_login(args)
        elif args.command == "logout":
            cmd_logout(args)
        elif args.command == "whoami":
            cmd_whoami(args)
        elif args.command == "config":
            cmd_config(args)
        elif args.command == "projects":
            cmd_projects(args)
        elif args.command == "jobs":
            cmd_jobs(args)
        elif args.command == "files":
            cmd_files(args)
        elif args.command == "report":
            cmd_report(args)
        elif args.command == "export":
            cmd_export(args)
        elif args.command == "import":
            cmd_import(args)
        elif args.command == "upload-request":
            cmd_upload_request(args)
        elif args.command == "upload-complete":
            cmd_upload_complete(args)
        elif args.command == "upload-status":
            cmd_upload_status(args)
        elif args.command == "upload-list":
            cmd_upload_list(args)
        elif args.command == "upload-reset":
            cmd_upload_reset(args)
        elif args.command == "upload-force-complete":
            cmd_upload_force_complete(args)
        else:
            print(f"Unknown command: {args.command}", file=sys.stderr)
            sys.exit(1)

    except APIError as e:
        print(f"API Error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nInterrupted", file=sys.stderr)
        sys.exit(130)


if __name__ == "__main__":
    main()
