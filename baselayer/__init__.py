"""
CCP4i2 BaseLayer - Environment-Aware Compatibility Layer

This module provides a unified API that works in both:
1. Legacy Qt/PySide2 environment (full CCP4 distribution)
2. Modern Django/Qt-free environment (ccp4i2-django)

Environment detection is automatic based on:
1. CCP4I2_BACKEND environment variable (explicit)
2. DJANGO_SETTINGS_MODULE environment variable
3. Availability of Django vs PySide2 imports (fallback)

Usage:
    from baselayer import QtCore, Signal, Slot, QObject
    # or
    from baselayer import DJANGO, QT
    if DJANGO():
        # Django-specific code
    else:
        # Qt-specific code
"""

import os
import sys

# Cached backend detection result
_backend = None


def _detect_backend():
    """
    Detect which backend to use.

    Priority:
    1. CCP4I2_BACKEND environment variable (explicit override)
    2. DJANGO_SETTINGS_MODULE presence (Django context)
    3. Import availability (try PySide2 first)

    Returns:
        str: 'django' or 'qt'
    """
    global _backend
    if _backend is not None:
        return _backend

    # Priority 1: Explicit environment variable
    explicit = os.environ.get('CCP4I2_BACKEND', '').lower()
    if explicit in ('django', 'modern', 'qt-free'):
        _backend = 'django'
        return _backend
    if explicit in ('qt', 'pyside2', 'legacy'):
        _backend = 'qt'
        return _backend

    # Priority 2: Django settings module present indicates Django context
    if os.environ.get('DJANGO_SETTINGS_MODULE'):
        _backend = 'django'
        return _backend

    # Priority 3: Try imports - prefer PySide2 if available for backward compatibility
    try:
        from PySide2 import QtCore as _QtCore
        _backend = 'qt'
    except ImportError:
        _backend = 'django'

    return _backend


def BACKEND():
    """
    Return current backend identifier.

    Returns:
        str: 'django' or 'qt'
    """
    return _detect_backend()


def DJANGO():
    """
    Check if using Django backend (Qt-free mode).

    Returns:
        bool: True if using Django backend
    """
    return _detect_backend() == 'django'


def QT():
    """
    Check if using Qt backend (legacy PySide2 mode).

    Returns:
        bool: True if using Qt/PySide2 backend
    """
    return _detect_backend() == 'qt'


# Conditional exports based on detected backend
if DJANGO():
    # Modern Qt-free implementations from our stubs
    from .QtCore import (
        Signal,
        Slot,
        QObject,
        QThread,
        QTimer,
        Qt,
        Property,
        QEventLoop,
        QCoreApplication,
        QApplication,
    )

    # Import Qt module stubs for "from baselayer import QtCore" usage
    from . import QtCore
    from . import QtGui
    from . import QtWidgets
    from . import QtSql

else:
    # Real PySide2 - import and re-export
    from PySide2.QtCore import (
        Signal,
        Slot,
        QObject,
        QThread,
        QTimer,
        Qt,
        Property,
        QEventLoop,
        QCoreApplication,
    )
    from PySide2.QtWidgets import QApplication

    # Import real Qt modules
    from PySide2 import QtCore
    from PySide2 import QtGui
    from PySide2 import QtWidgets
    from PySide2 import QtSql


# Version info
__version__ = "1.0.0"
__backend__ = BACKEND()
