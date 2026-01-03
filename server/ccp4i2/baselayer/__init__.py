# Modern Qt-free implementations from our stubs
from .QtCore import (
    Slot,
    QObject,
    Property,
)

# Import Qt module stubs for "from ccp4i2.baselayer import QtCore" usage
from . import QtCore

__all__ = [
    "Slot",
    "QObject",
    "Property",
    "QtCore",
]
