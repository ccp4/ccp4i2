"""
Centralized logging configuration for CCP4i2 codebase.

This module provides a consistent logging interface across the entire codebase,
replacing scattered print statements with proper Python logging.

Usage:
    from core.base_object.logging_config import get_logger
    logger = get_logger(__name__)
    logger.debug("Debug message")
    logger.info("Info message")
    logger.warning("Warning message")

Environment Variables:
    CCP4_LOG_LEVEL: Set log level (DEBUG, INFO, WARNING, ERROR, CRITICAL)
                    Default: INFO

Examples:
    # Normal usage (INFO and above)
    python manage.py create_job ...

    # Debug mode (see everything)
    export CCP4_LOG_LEVEL=DEBUG
    python manage.py create_job ...

    # Quiet mode (only errors, clean JSON output)
    export CCP4_LOG_LEVEL=ERROR
    python manage.py create_job --json
"""
import logging
import os
from typing import Optional

# Default log level from environment or INFO
DEFAULT_LOG_LEVEL = os.environ.get('CCP4_LOG_LEVEL', 'INFO').upper()

# Validate log level
VALID_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
if DEFAULT_LOG_LEVEL not in VALID_LEVELS:
    DEFAULT_LOG_LEVEL = 'INFO'


def get_logger(name: str, level: Optional[str] = None) -> logging.Logger:
    """
    Get a configured logger instance.

    Args:
        name: Logger name (typically __name__)
        level: Optional override level (DEBUG, INFO, WARNING, ERROR, CRITICAL)

    Returns:
        Configured logger instance

    Example:
        logger = get_logger(__name__)
        logger.info("Processing job %s", job_id)
    """
    logger = logging.getLogger(name)

    if not logger.handlers:
        # Only configure if not already configured
        handler = logging.StreamHandler()
        formatter = logging.Formatter(
            fmt='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    # Set level
    log_level = level or DEFAULT_LOG_LEVEL
    if log_level not in VALID_LEVELS:
        log_level = 'INFO'
    logger.setLevel(getattr(logging, log_level))

    return logger


def configure_quiet_mode():
    """
    Suppress all output except ERROR and CRITICAL.
    Useful for management commands that need clean JSON output.

    Example:
        from core.base_object.logging_config import configure_quiet_mode
        configure_quiet_mode()
        # Now only errors will be logged
    """
    logging.basicConfig(level=logging.ERROR, force=True)


def configure_debug_mode():
    """
    Enable verbose debug output for troubleshooting.

    Example:
        from core.base_object.logging_config import configure_debug_mode
        configure_debug_mode()
        # Now all debug messages will be shown
    """
    logging.basicConfig(level=logging.DEBUG, force=True)


def configure_silent_mode():
    """
    Suppress all logging output.
    Useful for unit tests or when only stdout output is desired.
    """
    logging.basicConfig(level=logging.CRITICAL + 1, force=True)


# Module-level logger for this configuration module
logger = get_logger(__name__)
logger.debug("Logging configured with level: %s", DEFAULT_LOG_LEVEL)
