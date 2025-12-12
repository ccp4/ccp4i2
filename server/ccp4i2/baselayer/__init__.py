# Modern Qt-free implementations from our stubs
from .QtCore import (
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

# Import Qt module stubs for "from ccp4i2.baselayer import QtCore" usage
from . import QtCore
