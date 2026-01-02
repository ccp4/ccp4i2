"""
Azure Service Bus queue utilities for async task processing.

This module provides utilities for sending messages to Azure Service Bus
for background processing by worker containers.
"""

import json
import logging
import os
from typing import Optional

logger = logging.getLogger(f"azure_extensions:{__name__}")


def send_to_queue(message_body: dict, queue_name: Optional[str] = None) -> bool:
    """
    Send a message to the Azure Service Bus queue.

    Args:
        message_body: Dictionary containing the message payload.
            Must include an 'action' key specifying the handler.
        queue_name: Optional queue name override. Defaults to
            SERVICE_BUS_QUEUE_NAME env var.

    Returns:
        True if message was sent successfully, False otherwise.

    Example:
        >>> send_to_queue({
        ...     "action": "import_project",
        ...     "upload_id": "550e8400-e29b-41d4-a716-446655440000",
        ...     "blob_path": "project_import/abc123.zip",
        ... })
        True
    """
    try:
        from azure.servicebus import ServiceBusClient, ServiceBusMessage
    except ImportError:
        logger.error(
            "azure-servicebus not installed. "
            "Install with: pip install azure-servicebus"
        )
        return False

    connection_string = os.environ.get("SERVICE_BUS_CONNECTION_STRING")
    if not connection_string:
        logger.error("SERVICE_BUS_CONNECTION_STRING not configured")
        return False

    if queue_name is None:
        queue_name = os.environ.get("SERVICE_BUS_QUEUE_NAME", "job-queue")

    try:
        with ServiceBusClient.from_connection_string(connection_string) as client:
            with client.get_queue_sender(queue_name) as sender:
                message = ServiceBusMessage(json.dumps(message_body))
                sender.send_messages(message)

        logger.info(
            "Queued message with action '%s' to queue '%s'",
            message_body.get("action", "unknown"),
            queue_name,
        )
        return True

    except Exception as e:
        logger.exception("Failed to send message to queue: %s", e)
        return False
