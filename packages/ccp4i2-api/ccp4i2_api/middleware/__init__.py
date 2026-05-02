"""Django middleware shipped by the shared auth package.

Each module in this package implements one auth scheme; consumers select
which to enable via Django ``MIDDLEWARE`` settings. All non-dev schemes
inherit from ``BaseAuthMiddleware`` (or will, after the AzureAD refactor),
which defines the canonical 401 response shape and the ``REQUEST_FLAG_ATTR``
trust signal.
"""

from .azure_ad import AzureADAuthMiddleware
from .base import REQUEST_FLAG_ATTR, BaseAuthMiddleware
from .dev import DevAuthMiddleware
from .dev_admin import DevAdminMiddleware
from .local_session import LocalSessionAuthMiddleware

__all__ = [
    "AzureADAuthMiddleware",
    "BaseAuthMiddleware",
    "DevAdminMiddleware",
    "DevAuthMiddleware",
    "LocalSessionAuthMiddleware",
    "REQUEST_FLAG_ATTR",
]
