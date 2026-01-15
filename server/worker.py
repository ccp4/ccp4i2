#!/usr/bin/env python3
"""
Standalone worker script to process jobs from Azure Service Bus queue.
This runs independently of Django and monitors the queue for new jobs.

Handles graceful shutdown on SIGTERM/SIGINT by marking in-progress jobs as FAILED.
"""
import os
import json
import time
import signal
import sys
import threading
import logging
import socket
from azure.servicebus import ServiceBusClient
from azure.identity import DefaultAzureCredential

# Get worker identity for logging
WORKER_ID = os.getenv("HOSTNAME", socket.gethostname())

# Track currently processing job for graceful shutdown
_current_job_uuid = None
_current_job_lock = threading.Lock()
_shutdown_requested = False

logging.basicConfig(
    level=logging.INFO,
    format=f"%(asctime)s - [{WORKER_ID}] - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


def get_subprocess_env():
    """
    Get environment dict for subprocesses with azure_packages in PYTHONPATH.

    The worker spawns ccp4-python subprocesses to run Django management commands.
    These need access to packages installed in the container (like django-storages)
    that aren't in the mounted CCP4 py-packages directory.

    The Dockerfile installs these packages to /usr/src/app/azure_packages/
    specifically so they can be added to PYTHONPATH without mixing Python
    site-packages between different Python installations.

    IMPORTANT: azure_packages is APPENDED (not prepended) so that packages from
    py-packages (like Django) take precedence. azure_packages only provides
    packages that aren't in py-packages (like django-storages).
    """
    env = os.environ.copy()

    # Path where Azure packages are installed in the Docker image
    azure_packages_path = "/usr/src/app/azure_packages"

    # Build new PYTHONPATH with azure_packages APPENDED (not prepended)
    # This ensures py-packages Django takes precedence over any Django in azure_packages
    existing_pythonpath = env.get("PYTHONPATH", "")

    if existing_pythonpath:
        env["PYTHONPATH"] = f"{existing_pythonpath}:{azure_packages_path}"
    else:
        env["PYTHONPATH"] = azure_packages_path

    logger.debug("Subprocess PYTHONPATH: %s", env["PYTHONPATH"])
    return env


def handle_shutdown_signal(signum, frame):
    """
    Handle SIGTERM/SIGINT for graceful shutdown.

    When Azure Container Apps scales down or restarts a container, it sends SIGTERM
    followed by SIGKILL after a grace period (default 30s). This handler:
    1. Sets shutdown flag to stop accepting new jobs
    2. Marks any in-progress job as FAILED so it doesn't get stuck in RUNNING state
    """
    global _shutdown_requested, _current_job_uuid

    signal_name = signal.Signals(signum).name
    logger.warning("Received %s signal - initiating graceful shutdown", signal_name)
    _shutdown_requested = True

    # If we have a job in progress, mark it as failed
    with _current_job_lock:
        if _current_job_uuid and _current_job_uuid != "unknown":
            logger.warning(
                "Marking in-progress job %s as FAILED due to worker shutdown",
                _current_job_uuid
            )
            try:
                update_job_status(_current_job_uuid, "FAILED")
                logger.info("Successfully marked job %s as FAILED", _current_job_uuid)
            except Exception as e:
                logger.error(
                    "Failed to mark job %s as FAILED during shutdown: %s",
                    _current_job_uuid, str(e)
                )

    logger.info("Graceful shutdown complete, exiting")
    sys.exit(0)


def set_current_job(job_uuid):
    """Track the currently processing job UUID for graceful shutdown handling."""
    global _current_job_uuid
    with _current_job_lock:
        _current_job_uuid = job_uuid


def clear_current_job():
    """Clear the current job tracker after job completes."""
    global _current_job_uuid
    with _current_job_lock:
        _current_job_uuid = None


def renew_lock_periodically(receiver, msg, stop_event, interval=30):
    """Thread target to renew Service Bus message lock periodically."""
    while not stop_event.is_set():
        try:
            receiver.renew_message_lock(msg)
            logger.debug("Renewed Service Bus message lock")
        except Exception as e:
            logger.error("Failed to renew message lock: %s", str(e))
        stop_event.wait(interval)


def process_job(job_data, receiver=None, msg=None):
    """
    Process a job from the queue.
    If receiver and msg are provided, start a thread to renew the lock during processing.
    """
    global _shutdown_requested

    job_uuid = job_data.get("uuid", job_data.get("job_uuid", "unknown"))
    action = job_data.get("action", "unknown")

    # Check if shutdown was requested before starting
    if _shutdown_requested:
        logger.warning("Shutdown requested, not processing job %s", job_uuid)
        return False

    # Track this job for graceful shutdown handling
    set_current_job(job_uuid)

    logger.info(
        "=== WORKER %s STARTING JOB %s (action: %s) ===", WORKER_ID, job_uuid, action
    )

    lock_stop_event = threading.Event()
    lock_thread = None

    try:
        if receiver is not None and msg is not None:
            # Start lock renewal thread
            lock_thread = threading.Thread(
                target=renew_lock_periodically, args=(receiver, msg, lock_stop_event)
            )
            lock_thread.daemon = True
            lock_thread.start()

        # Route to appropriate handler based on action
        if action == "run_job":
            # Run CCP4 analysis
            result = run_ccp4_analysis(job_data)
            logger.info("CCP4 analysis completed for job %s", job_uuid)

            # Update job status based on analysis result
            if result.get("status") == "completed":
                # update_job_status(job_uuid, "FINISHED")
                # I think that succesfully completed jobs are already marked as FINISHED by manage.py run_job
                return True
            elif result.get("status") == "failed":
                update_job_status(job_uuid, "FAILED")
                return False
            else:
                logger.warning("Unknown result status: %s", result.get("status"))
                update_job_status(job_uuid, "FAILED")
                return False

        elif action == "import_project":
            # Import a project from staged upload
            result = process_project_import(job_data)
            return result.get("status") == "completed"

        elif action == "process_unmerged_data":
            # Process unmerged data from staged upload
            result = process_unmerged_data(job_data)
            return result.get("status") == "completed"

        else:
            logger.warning("Unknown action type: %s", action)
            raise ValueError("Unsupported action: %s" % action)

    except (ValueError, KeyError) as e:
        logger.error("Error processing job %s: %s", job_uuid, str(e))
        # Ensure job status is updated to FAILED on exception
        try:
            update_job_status(job_uuid, "FAILED")
        except Exception as status_error:
            logger.error(
                "Failed to update job status after error: %s", str(status_error)
            )
        return False
    except Exception as e:
        # Catch any unexpected exceptions
        logger.error("Unexpected error processing job %s: %s", job_uuid, str(e))
        # Ensure job status is updated to FAILED
        try:
            update_job_status(job_uuid, "FAILED")
        except Exception as status_error:
            logger.error(
                "Failed to update job status after error: %s", str(status_error)
            )
        return False
    finally:
        # Clear job tracking for graceful shutdown
        clear_current_job()
        if lock_thread is not None:
            lock_stop_event.set()
            lock_thread.join()


def run_ccp4_analysis(parameters):
    """Run CCP4 analysis using the configured Python executable"""
    import subprocess
    import os

    logger.info("Running CCP4 analysis with parameters: %s", parameters)

    # Get the CCP4 Python executable from environment
    ccp4_python = os.getenv("CCP4_PYTHON")
    if not ccp4_python:
        error_msg = "CCP4_PYTHON environment variable not set"
        logger.error(error_msg)
        return {"status": "failed", "error": error_msg}

    job_uuid = parameters.get("job_uuid")
    if not job_uuid:
        error_msg = "job_uuid not provided in parameters"
        logger.error(error_msg)
        return {"status": "failed", "error": error_msg}

    # Prepare command arguments
    cmd = [ccp4_python, "/usr/src/app/manage.py", "run_job", "-ju", job_uuid]

    logger.info("Executing command: %s", " ".join(cmd))

    # Get environment with container site-packages in PYTHONPATH
    env = get_subprocess_env()
    logger.info("PYTHONPATH passed to subprocess: %s", env.get("PYTHONPATH", ""))

    try:
        # Run the command with environment including container packages
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour timeout
            cwd="/usr/src/app",
            check=False,  # We handle return codes manually
            env=env,
        )

        logger.info("Command completed with return code: %s", result.returncode)

        if result.returncode == 0:
            logger.info("CCP4 analysis completed successfully")
            return {
                "status": "completed",
                "result": "analysis_result",
                "stdout": result.stdout,
                "stderr": result.stderr,
            }
        else:
            logger.error("CCP4 analysis failed with return code %s", result.returncode)
            logger.error("STDOUT: %s", result.stdout)
            logger.error("STDERR: %s", result.stderr)
            return {
                "status": "failed",
                "error": "Command failed with return code %s" % result.returncode,
                "stdout": result.stdout,
                "stderr": result.stderr,
            }

    except subprocess.TimeoutExpired:
        error_msg = "CCP4 analysis timed out after 1 hour"
        logger.error(error_msg)
        return {"status": "failed", "error": error_msg}

    except FileNotFoundError:
        error_msg = "CCP4 Python executable not found: %s" % ccp4_python
        logger.error(error_msg)
        return {"status": "failed", "error": error_msg}

    except Exception as e:
        error_msg = "Unexpected error running CCP4 analysis: %s" % str(e)
        logger.error(error_msg)
        return {"status": "failed", "error": error_msg}


def validate_zip_file(file_path):
    """
    Validate that a file is a valid zip archive.

    Checks the zip file header and attempts to read the central directory
    to detect corrupted or incomplete uploads.

    Returns:
        True if file is a valid zip, False otherwise
    """
    import zipfile

    if not os.path.exists(file_path):
        logger.error("File does not exist: %s", file_path)
        return False

    file_size = os.path.getsize(file_path)
    if file_size < 22:  # Minimum size for an empty zip file
        logger.error("File too small to be a valid zip: %d bytes", file_size)
        return False

    try:
        with zipfile.ZipFile(file_path, "r") as zf:
            # Test the zip file integrity by reading the file list
            file_list = zf.namelist()
            logger.info("Zip file validated: %d files found", len(file_list))
            return True
    except zipfile.BadZipFile as e:
        logger.error("Invalid zip file: %s", str(e))
        return False
    except Exception as e:
        logger.error("Error validating zip file: %s", str(e))
        return False


def process_project_import(parameters):
    """
    Process a project import from a staged upload.

    Downloads the blob from Azure Storage and runs the import command.
    """
    import subprocess

    upload_id = parameters.get("upload_id")
    blob_path = parameters.get("blob_path")
    original_filename = parameters.get("original_filename", "unknown.zip")

    logger.info(
        "Processing project import: upload_id=%s, blob_path=%s",
        upload_id, blob_path
    )

    if not blob_path:
        error_msg = "blob_path not provided in parameters"
        logger.error(error_msg)
        update_staged_upload_status(upload_id, "failed", error_msg)
        return {"status": "failed", "error": error_msg}

    try:
        # Download blob to local temp file
        local_path = download_blob_to_local(blob_path)
        logger.info("Downloaded blob to %s", local_path)

        # Validate zip file before attempting import
        if not validate_zip_file(local_path):
            error_msg = "Downloaded file is not a valid zip file (corrupted or incomplete upload)"
            logger.error(error_msg)
            update_staged_upload_status(upload_id, "failed", error_msg)
            # Clean up the invalid file
            try:
                os.remove(local_path)
            except Exception:
                pass
            return {"status": "failed", "error": error_msg}

        # Get CCP4 Python executable
        ccp4_python = os.getenv("CCP4_PYTHON")
        if not ccp4_python:
            error_msg = "CCP4_PYTHON environment variable not set"
            logger.error(error_msg)
            update_staged_upload_status(upload_id, "failed", error_msg)
            return {"status": "failed", "error": error_msg}

        # Run import command
        cmd = [
            ccp4_python,
            "/usr/src/app/manage.py",
            "import_ccp4_project_zip",
            local_path,
        ]

        logger.info("Executing import command: %s", " ".join(cmd))

        # Get environment with azure_packages in PYTHONPATH
        env = get_subprocess_env()
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=7200,  # 2 hour timeout for large imports
            cwd="/usr/src/app",
            check=False,
            env=env,
        )

        if result.returncode == 0:
            logger.info("Project import completed successfully")
            update_staged_upload_status(upload_id, "completed")

            # Clean up: delete blob and local file
            cleanup_after_import(blob_path, local_path)

            return {"status": "completed"}
        else:
            error_msg = "Import command failed with return code %s: %s" % (
                result.returncode, result.stderr
            )
            logger.error(error_msg)
            logger.error("STDOUT: %s", result.stdout)
            update_staged_upload_status(upload_id, "failed", error_msg)
            return {"status": "failed", "error": error_msg}

    except subprocess.TimeoutExpired:
        error_msg = "Project import timed out after 2 hours"
        logger.error(error_msg)
        update_staged_upload_status(upload_id, "failed", error_msg)
        return {"status": "failed", "error": error_msg}

    except Exception as e:
        error_msg = "Error processing project import: %s" % str(e)
        logger.exception(error_msg)
        update_staged_upload_status(upload_id, "failed", error_msg)
        return {"status": "failed", "error": error_msg}


def process_unmerged_data(parameters):
    """
    Process an unmerged data upload.

    Downloads the blob and moves it to the target job's directory.
    """
    upload_id = parameters.get("upload_id")
    blob_path = parameters.get("blob_path")
    target_job_uuid = parameters.get("target_job_uuid")

    logger.info(
        "Processing unmerged data: upload_id=%s, target_job=%s",
        upload_id, target_job_uuid
    )

    if not blob_path or not target_job_uuid:
        error_msg = "Missing required parameters: blob_path or target_job_uuid"
        logger.error(error_msg)
        update_staged_upload_status(upload_id, "failed", error_msg)
        return {"status": "failed", "error": error_msg}

    try:
        # Download blob to local temp file
        local_path = download_blob_to_local(blob_path)
        logger.info("Downloaded blob to %s", local_path)

        # TODO: Move file to job directory and update job parameters
        # This would involve:
        # 1. Finding the job's directory
        # 2. Moving the file there
        # 3. Updating job parameters via manage.py command

        update_staged_upload_status(upload_id, "completed")

        # Clean up blob (keep local file for job processing)
        try:
            delete_blob(blob_path)
        except Exception as e:
            logger.warning("Failed to delete blob %s: %s", blob_path, e)

        return {"status": "completed", "local_path": local_path}

    except Exception as e:
        error_msg = "Error processing unmerged data: %s" % str(e)
        logger.exception(error_msg)
        update_staged_upload_status(upload_id, "failed", error_msg)
        return {"status": "failed", "error": error_msg}


def download_blob_to_local(blob_path, download_dir="/tmp/staged-uploads"):
    """
    Download a blob from Azure Storage to a local path.

    Uses system Python (python3) via subprocess because ccp4-python doesn't
    have azure-storage-blob installed, and the shared /mnt/ccp4data volume
    shouldn't be modified at runtime.

    Args:
        blob_path: Path to the blob within the staging container
        download_dir: Local directory to download to

    Returns:
        Local file path where blob was downloaded
    """
    import subprocess

    storage_account_name = os.environ.get("AZURE_STORAGE_ACCOUNT_NAME")
    if not storage_account_name:
        raise ValueError("AZURE_STORAGE_ACCOUNT_NAME not configured")

    # Create download directory
    os.makedirs(download_dir, exist_ok=True)

    # Extract filename from blob path
    filename = blob_path.split("/")[-1]
    local_path = os.path.join(download_dir, filename)

    logger.info("Downloading blob %s to %s", blob_path, local_path)

    # Use system Python which has azure-storage-blob installed
    # Pass env vars for managed identity
    download_script = '''
import os
import sys
from azure.storage.blob import BlobServiceClient
from azure.identity import DefaultAzureCredential

storage_account_name = os.environ.get("AZURE_STORAGE_ACCOUNT_NAME")
container_name = "staging-uploads"
blob_path = sys.argv[1]
local_path = sys.argv[2]

managed_identity_client_id = os.environ.get("AZURE_CLIENT_ID")
if managed_identity_client_id:
    credential = DefaultAzureCredential(managed_identity_client_id=managed_identity_client_id)
else:
    credential = DefaultAzureCredential()

account_url = f"https://{storage_account_name}.blob.core.windows.net"
blob_service_client = BlobServiceClient(account_url=account_url, credential=credential)

blob_client = blob_service_client.get_blob_client(container=container_name, blob=blob_path)

# Stream download in chunks to avoid loading entire file into memory
download_stream = blob_client.download_blob()
total_bytes = 0
with open(local_path, "wb") as f:
    for chunk in download_stream.chunks():
        f.write(chunk)
        total_bytes += len(chunk)

print(f"Downloaded {total_bytes} bytes")
'''

    result = subprocess.run(
        ["python3", "-c", download_script, blob_path, local_path],
        capture_output=True,
        text=True,
        timeout=7200,  # 2 hour timeout for large files
        env=os.environ.copy(),
    )

    if result.returncode != 0:
        raise RuntimeError(f"Blob download failed: {result.stderr}")

    file_size = os.path.getsize(local_path)
    logger.info("Downloaded %s (%d bytes)", blob_path, file_size)

    return local_path


def delete_blob(blob_path):
    """Delete a blob from the staging container using system Python."""
    import subprocess

    storage_account_name = os.environ.get("AZURE_STORAGE_ACCOUNT_NAME")
    if not storage_account_name:
        raise ValueError("AZURE_STORAGE_ACCOUNT_NAME not configured")

    delete_script = '''
import os
import sys
from azure.storage.blob import BlobServiceClient
from azure.identity import DefaultAzureCredential

storage_account_name = os.environ.get("AZURE_STORAGE_ACCOUNT_NAME")
container_name = "staging-uploads"
blob_path = sys.argv[1]

managed_identity_client_id = os.environ.get("AZURE_CLIENT_ID")
if managed_identity_client_id:
    credential = DefaultAzureCredential(managed_identity_client_id=managed_identity_client_id)
else:
    credential = DefaultAzureCredential()

account_url = f"https://{storage_account_name}.blob.core.windows.net"
blob_service_client = BlobServiceClient(account_url=account_url, credential=credential)

blob_client = blob_service_client.get_blob_client(container=container_name, blob=blob_path)
blob_client.delete_blob()
print(f"Deleted blob: {blob_path}")
'''

    result = subprocess.run(
        ["python3", "-c", delete_script, blob_path],
        capture_output=True,
        text=True,
        timeout=60,
        env=os.environ.copy(),
    )

    if result.returncode != 0:
        raise RuntimeError(f"Blob deletion failed: {result.stderr}")

    logger.info("Deleted blob: %s", blob_path)


def cleanup_after_import(blob_path, local_path):
    """Clean up after successful import."""
    # Delete blob
    try:
        delete_blob(blob_path)
    except Exception as e:
        logger.warning("Failed to delete blob %s: %s", blob_path, e)

    # Delete local file
    try:
        if os.path.exists(local_path):
            os.remove(local_path)
            logger.info("Deleted local file: %s", local_path)
    except Exception as e:
        logger.warning("Failed to delete local file %s: %s", local_path, e)


def update_staged_upload_status(upload_id, status, error_message=None):
    """
    Update the status of a staged upload in the database.

    Uses Django ORM via management command to avoid importing Django in worker.
    """
    import subprocess

    ccp4_python = os.getenv("CCP4_PYTHON")
    if not ccp4_python:
        logger.error("CCP4_PYTHON not set, cannot update upload status")
        return False

    # Use a simple management command to update status
    cmd = [
        ccp4_python,
        "/usr/src/app/manage.py",
        "update_staged_upload",
        "--upload-id", upload_id,
        "--status", status,
    ]
    if error_message:
        cmd.extend(["--error", error_message])

    logger.info("Updating staged upload %s to status %s", upload_id, status)

    # Get environment with container site-packages in PYTHONPATH
    env = get_subprocess_env()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,
            cwd="/usr/src/app",
            check=False,
            env=env,
        )

        if result.returncode == 0:
            logger.info("Updated staged upload status successfully")
            return True
        else:
            logger.error("Failed to update staged upload status: %s", result.stderr)
            return False

    except Exception as e:
        logger.error("Error updating staged upload status: %s", e)
        return False


def update_job_status(job_uuid, status, result=None):
    """Update job status using the set_job_status management command"""

    import subprocess

    logger.info("Updating job %s status to %s", job_uuid, status)

    # Get the CCP4 Python executable from environment
    ccp4_python = os.getenv("CCP4_PYTHON")
    if not ccp4_python:
        logger.error("CCP4_PYTHON environment variable not set")
        return False

    # Prepare command arguments for set_job_status management command
    cmd = [
        ccp4_python,
        "/usr/src/app/manage.py",
        "set_job_status",
        "-ju",
        job_uuid,
        "-s",
        status.upper(),  # Status should be uppercase for the command
    ]

    logger.info("Executing status update command: %s", " ".join(cmd))

    # Get environment with container site-packages in PYTHONPATH
    env = get_subprocess_env()
    logger.info("PYTHONPATH passed to subprocess: %s", env.get("PYTHONPATH", ""))

    try:
        # Run the command with environment including container packages
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=30,  # 30 second timeout for status updates
            cwd="/usr/src/app",
            check=False,
            env=env,
        )

        if result.returncode == 0:
            logger.info("Job status updated successfully for job %s", job_uuid)
            return True
        else:
            logger.error("Failed to update job status for job %s", job_uuid)
            logger.error("STDOUT: %s", result.stdout)
            logger.error("STDERR: %s", result.stderr)
            return False

    except subprocess.TimeoutExpired:
        logger.error("Status update timed out for job %s", job_uuid)
        return False
    except OSError as e:
        logger.error("OS error updating job status for job %s: %s", job_uuid, str(e))
        return False


def create_service_bus_client(connection_string):
    """Create and return a Service Bus client"""
    if connection_string.startswith("https://"):
        # Key Vault reference - use managed identity
        logger.info("Using managed identity for Service Bus authentication")
        credential = DefaultAzureCredential()
        return ServiceBusClient(
            fully_qualified_namespace=connection_string, credential=credential
        )
    else:
        # Direct connection string
        logger.info("Using connection string for Service Bus authentication")
        return ServiceBusClient.from_connection_string(connection_string)


def run_worker_loop(sb_client, queue_name):
    """Run the main worker message processing loop"""
    with sb_client:
        receiver = sb_client.get_queue_receiver(queue_name=queue_name)

        with receiver:
            logger.info("Worker ready to process jobs...")

            while True:
                try:
                    # Receive messages with timeout
                    messages = receiver.receive_messages(
                        max_message_count=1, max_wait_time=30
                    )

                    for msg in messages:
                        try:
                            job_data = json.loads(str(msg))
                            success = process_job(job_data, receiver=receiver, msg=msg)

                            if success:
                                # Mark message as completed to remove it from queue
                                receiver.complete_message(msg)
                                logger.info(
                                    "Job processed and message completed successfully"
                                )
                            else:
                                # Job failed, abandon message so it can be retried
                                receiver.abandon_message(msg)
                                logger.warning(
                                    "Job processing failed, message abandoned for retry"
                                )

                        except (json.JSONDecodeError, ValueError, KeyError) as e:
                            logger.error("Error processing job: %s", str(e))
                            # Send to dead-letter queue for manual inspection
                            try:
                                receiver.dead_letter_message(msg, reason=str(e))
                                logger.info("Invalid message sent to dead-letter queue")
                            except OSError as dlq_error:
                                logger.error(
                                    "Failed to dead-letter message: %s", str(dlq_error)
                                )
                        except Exception as e:
                            # Unexpected error - abandon message for retry
                            logger.error("Unexpected error processing job: %s", str(e))
                            try:
                                receiver.abandon_message(msg)
                            except OSError as abandon_error:
                                logger.error(
                                    "Failed to abandon message: %s", str(abandon_error)
                                )

                except OSError as e:
                    logger.error("Error receiving messages: %s", str(e))
                    time.sleep(5)  # Brief pause before retry


def cleanup_stale_jobs(stale_threshold_hours=2):
    """
    Clean up jobs stuck in RUNNING or RUNNING_REMOTELY state.

    This handles cases where:
    - Worker was OOM-killed without catching signal
    - Container crashed unexpectedly
    - Network partition prevented status update

    Args:
        stale_threshold_hours: Jobs running longer than this are considered stale
    """
    import subprocess

    logger.info("Checking for stale jobs (threshold: %d hours)...", stale_threshold_hours)

    ccp4_python = os.getenv("CCP4_PYTHON")
    if not ccp4_python:
        logger.warning("CCP4_PYTHON not set, skipping stale job cleanup")
        return

    cmd = [
        ccp4_python,
        "/usr/src/app/manage.py",
        "cleanup_stale_jobs",
        "--hours", str(stale_threshold_hours),
    ]

    # Get environment with azure_packages in PYTHONPATH
    env = get_subprocess_env()
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=60,
            cwd="/usr/src/app",
            check=False,
            env=env,
        )

        if result.returncode == 0:
            logger.info("Stale job cleanup completed: %s", result.stdout.strip())
        else:
            logger.warning("Stale job cleanup failed: %s", result.stderr)

    except subprocess.TimeoutExpired:
        logger.warning("Stale job cleanup timed out")
    except FileNotFoundError:
        logger.warning("cleanup_stale_jobs command not found - skipping")
    except Exception as e:
        logger.warning("Error during stale job cleanup: %s", e)


def main():
    """Main worker loop"""
    # Register signal handlers for graceful shutdown
    signal.signal(signal.SIGTERM, handle_shutdown_signal)
    signal.signal(signal.SIGINT, handle_shutdown_signal)
    logger.info("Signal handlers registered for graceful shutdown")

    # Get configuration from environment
    queue_name = os.getenv("SERVICE_BUS_QUEUE_NAME", "ccp4i2-bicep-jobs")
    connection_string = os.getenv("SERVICE_BUS_CONNECTION_STRING")

    if not connection_string:
        logger.error("SERVICE_BUS_CONNECTION_STRING environment variable not set")
        return

    logger.info("Starting worker for queue: %s", queue_name)

    # Clean up any stale jobs from previous worker crashes on startup
    cleanup_stale_jobs()

    # Initialize Service Bus client
    try:
        sb_client = create_service_bus_client(connection_string)
        run_worker_loop(sb_client, queue_name)

    except KeyboardInterrupt:
        logger.info("Worker stopped by user")
    except OSError as e:
        logger.error("Worker error: %s", str(e))
        raise


if __name__ == "__main__":
    main()
