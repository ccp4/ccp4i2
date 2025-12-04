#!/usr/bin/env python3
"""
Standalone worker script to process jobs from Azure Service Bus queue.
This runs independently of Django and monitors the queue for new jobs.
"""
import os
import json
import time
import threading
import logging
import socket
from azure.servicebus import ServiceBusClient
from azure.identity import DefaultAzureCredential

# Get worker identity for logging
WORKER_ID = os.getenv("HOSTNAME", socket.gethostname())

logging.basicConfig(
    level=logging.INFO,
    format=f"%(asctime)s - [{WORKER_ID}] - %(levelname)s - %(message)s",
)
logger = logging.getLogger(__name__)


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
    job_uuid = job_data.get("uuid", job_data.get("job_uuid", "unknown"))
    action = job_data.get("action", "unknown")

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

        # Add your job processing logic here
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

        else:
            logger.warning("Unknown action type: %s", action)
            update_job_status(job_uuid, "FAILED")
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

    # Ensure subprocess inherits PYTHONPATH
    env = os.environ.copy()
    if "PYTHONPATH" in env:
        logger.info("PYTHONPATH passed to subprocess: %s", env["PYTHONPATH"])
    else:
        logger.warning("PYTHONPATH not set in environment")

    try:
        # Run the command with inherited environment
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=3600,  # 1 hour timeout
            cwd="/usr/src/app",
            check=False,  # We handle return codes manually
            env=env,  # <-- ADD THIS LINE
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

    # Ensure subprocess inherits PYTHONPATH
    env = os.environ.copy()
    if "PYTHONPATH" in env:
        logger.info("PYTHONPATH passed to subprocess: %s", env["PYTHONPATH"])
    else:
        logger.warning("PYTHONPATH not set in environment")

    try:
        # Run the command
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


def main():
    """Main worker loop"""
    # Get configuration from environment
    queue_name = os.getenv("SERVICE_BUS_QUEUE_NAME", "ccp4i2-bicep-jobs")
    connection_string = os.getenv("SERVICE_BUS_CONNECTION_STRING")

    if not connection_string:
        logger.error("SERVICE_BUS_CONNECTION_STRING environment variable not set")
        return

    logger.info("Starting worker for queue: %s", queue_name)

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
