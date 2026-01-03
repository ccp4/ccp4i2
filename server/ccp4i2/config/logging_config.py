"""
Django-specific logging configuration for CCP4i2.

This module provides Django LOGGING settings that integrate with the core
logging infrastructure from ccp4i2.core.base_object.logging_config.

Import this in your Django settings.py:
    from ccp4i2.config.logging_config import LOGGING
"""
import os
from pathlib import Path

# Determine log file location
# Default to project root, but allow override
LOG_DIR = os.environ.get('CCP4_LOG_DIR', str(Path(__file__).resolve().parent.parent.parent.parent))
LOG_FILE = os.path.join(LOG_DIR, 'ccp4i2.log')

# Get log level from environment
LOG_LEVEL = os.environ.get('CCP4_LOG_LEVEL', 'INFO').upper()

# Validate log level
VALID_LEVELS = ['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL']
if LOG_LEVEL not in VALID_LEVELS:
    LOG_LEVEL = 'INFO'

# Django LOGGING configuration
LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
            'datefmt': '%Y-%m-%d %H:%M:%S'
        },
        'simple': {
            'format': '%(levelname)s - %(message)s'
        },
        'json': {
            'format': '{"time": "%(asctime)s", "name": "%(name)s", "level": "%(levelname)s", "message": "%(message)s"}',
            'datefmt': '%Y-%m-%d %H:%M:%S'
        },
    },
    'filters': {
        'require_debug_false': {
            '()': 'django.utils.log.RequireDebugFalse',
        },
        'require_debug_true': {
            '()': 'django.utils.log.RequireDebugTrue',
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'verbose',
            'level': LOG_LEVEL,
        },
        'console_simple': {
            'class': 'logging.StreamHandler',
            'formatter': 'simple',
            'level': LOG_LEVEL,
        },
        'file': {
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': LOG_FILE,
            'maxBytes': 10485760,  # 10MB
            'backupCount': 5,
            'formatter': 'verbose',
            'level': 'DEBUG',  # Always log debug to file
        },
        'null': {
            'class': 'logging.NullHandler',
        },
    },
    'loggers': {
        # CCP4i2 application loggers
        'ccp4i2': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        'ccp4i2.api': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        'ccp4i2.db': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        'ccp4i2.lib': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        # Core CData loggers
        'core': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        'core.base_object': {
            'handlers': ['console', 'file'],
            'level': 'DEBUG',
            'propagate': False,
        },
        # Django framework loggers
        'django': {
            'handlers': ['console', 'file'],
            'level': 'INFO',
            'propagate': False,
        },
        'django.server': {
            'handlers': ['console_simple'],
            'level': 'INFO',
            'propagate': False,
        },
        'django.request': {
            'handlers': ['console', 'file'],
            'level': 'ERROR',
            'propagate': False,
        },
        'django.db.backends': {
            'handlers': ['null'],  # Suppress SQL queries by default
            'level': 'DEBUG',
            'propagate': False,
        },
    },
    'root': {
        'handlers': ['console', 'file'],
        'level': LOG_LEVEL,
    },
}


def configure_for_management_command(quiet: bool = False):
    """
    Adjust logging for management command execution.

    Args:
        quiet: If True, suppress all logging except errors
               (useful for commands that output JSON)
    """
    import logging

    if quiet:
        logging.basicConfig(level=logging.ERROR, force=True)
    else:
        logging.basicConfig(level=getattr(logging, LOG_LEVEL), force=True)
