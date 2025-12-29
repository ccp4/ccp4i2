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

    i2remote files <project> <job> [list] List files in a job
    i2remote files cat <project> <job> <name>  Display file contents

    i2remote report <project> <job>       Get job report

    i2remote export job <project> <job>   Export a job to zip
    i2remote export project <project>     Export a project to zip
    i2remote import <zipfile>             Import a project from zip

    i2remote config                       Show current configuration
    i2remote config set <key> <value>     Set configuration value

Environment Variables:
    CCP4I2_API_URL      Base URL for the CCP4i2 API (required)
    CCP4I2_API_TOKEN    Authentication token (if required)

Configuration File:
    ~/.ccp4i2remote.json - Persistent configuration

Examples:
    # Set remote server
    export CCP4I2_API_URL=https://myserver.azurecontainerapps.io/api/proxy

    # List projects
    i2remote projects

    # Create and run a job
    i2remote jobs create myproject aimless
    i2remote jobs run myproject 1

    # Wait for completion and check results
    i2remote jobs wait myproject 1
    i2remote jobs kpi myproject 1
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

    def import_project(self, project_id: str, zip_path: str) -> dict:
        """Import project from ZIP."""
        with open(zip_path, "rb") as f:
            files = {"file": (Path(zip_path).name, f, "application/zip")}
            # Need to use form data for file upload
            return self._request("POST", f"/projects/{project_id}/import_project", files=files)

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

def cmd_config(args):
    """Handle config command."""
    if args.action == "show" or args.action is None:
        print("Current configuration:")
        for key, value in CONFIG.items():
            if key == "api_token" and value:
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
    if action and action not in ("list", "create", "run", "clone", "status", "wait", "kpi", "tree", "validate"):
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

    if not Path(args.zipfile).exists():
        print(f"Error: File not found: {args.zipfile}", file=sys.stderr)
        sys.exit(1)

    # Import requires a project - create one or specify
    print("Note: Project import creates a new project from the ZIP")
    # This would need the import endpoint which creates a new project
    print("Import functionality requires server-side implementation")


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
    jobs_parser.add_argument("action", nargs="?", help="list|create|run|clone|status|wait|kpi|tree")
    jobs_parser.add_argument("project", nargs="?", help="Project name/ID")
    jobs_parser.add_argument("job", nargs="?", help="Job number or ID")
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

    args = parser.parse_args()

    if not args.command:
        parser.print_help()
        sys.exit(0)

    try:
        if args.command == "config":
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
