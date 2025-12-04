"""
Stub module for legacy dbapi (CCP4DbApi).

This stub allows legacy report classes that import from dbapi to be loaded
for registry discovery. The actual database functionality is provided by
the new Django-based ccp4x.db system.

Note: Calling actual dbapi methods will raise NotImplementedError.
"""

import uuid


class UUIDTYPE:
    """Stub UUID type for CCP4DbApi compatibility."""
    def __init__(self, value):
        if isinstance(value, uuid.UUID):
            self._uuid = value
        elif isinstance(value, str):
            self._uuid = uuid.UUID(value)
        else:
            self._uuid = uuid.UUID(str(value))

    def __str__(self):
        return str(self._uuid)

    def __repr__(self):
        return f"UUIDTYPE('{self._uuid}')"


class CCP4DbApi:
    """
    Stub for the legacy CCP4DbApi class.

    The legacy dbapi provided database access for the Qt-based ccp4i2.
    In the new system, database access is handled by Django models in ccp4x.db.
    """
    UUIDTYPE = UUIDTYPE

    def __init__(self, *args, **kwargs):
        raise NotImplementedError(
            "CCP4DbApi is a legacy component. "
            "Use ccp4x.db.models for database access in the new system."
        )

    @staticmethod
    def getProjectId(*args, **kwargs):
        raise NotImplementedError("Use ccp4x.db.models.Project instead")

    @staticmethod
    def getJobInfo(*args, **kwargs):
        raise NotImplementedError("Use ccp4x.db.models.Job instead")


# For compatibility with imports like "from dbapi import CCP4DbApi"
__all__ = ['CCP4DbApi', 'UUIDTYPE']
