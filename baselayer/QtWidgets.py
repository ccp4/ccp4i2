"""
Stub QtWidgets module for plugin discovery.

Provides minimal Qt widget stubs to allow plugin files to be imported
without requiring PySide2.
"""

from .QtCore import QObject


class QWidget(QObject):
    """Stub QWidget class."""
    pass


class QMainWindow(QWidget):
    """Stub QMainWindow class."""
    pass


class QDialog(QWidget):
    """Stub QDialog class."""
    pass


class QLineEdit(QWidget):
    """Stub QLineEdit class."""
    pass


class QPushButton(QWidget):
    """Stub QPushButton class."""
    pass


class QLabel(QWidget):
    """Stub QLabel class."""
    pass


class QComboBox(QWidget):
    """Stub QComboBox class."""
    pass


class QCheckBox(QWidget):
    """Stub QCheckBox class."""
    pass


class QTableWidget(QWidget):
    """Stub QTableWidget class."""
    pass


class QVBoxLayout:
    """Stub QVBoxLayout class."""
    pass


class QHBoxLayout:
    """Stub QHBoxLayout class."""
    pass


class QApplication:
    """Stub QApplication class."""
    pass
