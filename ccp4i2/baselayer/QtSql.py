"""
Stub PySide2.QtSql module for compatibility.

This provides minimal stubs for Qt SQL classes to allow imports to work
without requiring the full PySide2 installation.
"""

class QSqlDatabase:
    """Stub for QSqlDatabase."""

    @staticmethod
    def addDatabase(driver_type, connection_name=None):
        """Stub for addDatabase."""
        return QSqlDatabase()

    @staticmethod
    def database(connection_name=None):
        """Stub for database."""
        return QSqlDatabase()

    @staticmethod
    def removeDatabase(connection_name):
        """Stub for removeDatabase."""
        pass

    def setDatabaseName(self, name):
        """Stub for setDatabaseName."""
        pass

    def setHostName(self, host):
        """Stub for setHostName."""
        pass

    def setPort(self, port):
        """Stub for setPort."""
        pass

    def setUserName(self, user):
        """Stub for setUserName."""
        pass

    def setPassword(self, password):
        """Stub for setPassword."""
        pass

    def open(self):
        """Stub for open."""
        return True

    def close(self):
        """Stub for close."""
        pass

    def isOpen(self):
        """Stub for isOpen."""
        return False

    def isValid(self):
        """Stub for isValid."""
        return False

    def exec_(self, query):
        """Stub for exec_."""
        return QSqlQuery()

    def tables(self, table_type=None):
        """Stub for tables."""
        return []


class QSqlQuery:
    """Stub for QSqlQuery."""

    def __init__(self, query=None, db=None):
        """Initialize stub query."""
        pass

    def exec_(self, query=None):
        """Stub for exec_."""
        return True

    def next(self):
        """Stub for next."""
        return False

    def value(self, index):
        """Stub for value."""
        return None

    def isActive(self):
        """Stub for isActive."""
        return False

    def isValid(self):
        """Stub for isValid."""
        return False


class QSqlError:
    """Stub for QSqlError."""

    def text(self):
        """Stub for text."""
        return ""

    def type(self):
        """Stub for type."""
        return 0


class QSqlRecord:
    """Stub for QSqlRecord."""

    def count(self):
        """Stub for count."""
        return 0

    def fieldName(self, index):
        """Stub for fieldName."""
        return ""

    def value(self, index):
        """Stub for value."""
        return None


class QSqlTableModel:
    """Stub for QSqlTableModel."""

    def __init__(self, parent=None, db=None):
        """Initialize stub model."""
        pass

    def setTable(self, table_name):
        """Stub for setTable."""
        pass

    def select(self):
        """Stub for select."""
        return True

    def rowCount(self, parent=None):
        """Stub for rowCount."""
        return 0
