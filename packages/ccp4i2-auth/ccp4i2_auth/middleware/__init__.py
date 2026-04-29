"""Django middleware shipped by the shared auth package.

Each module in this package implements one auth scheme; consumers select
which to enable via Django ``MIDDLEWARE`` settings. All non-dev schemes
inherit from ``BaseAuthMiddleware``, which defines the canonical 401
response shape and the ``REQUEST_FLAG_ATTR`` trust signal.
"""

from .base import REQUEST_FLAG_ATTR, BaseAuthMiddleware
from .dev import DevAuthMiddleware
from .local_session import LocalSessionAuthMiddleware

__all__ = [
    "BaseAuthMiddleware",
    "DevAuthMiddleware",
    "LocalSessionAuthMiddleware",
    "REQUEST_FLAG_ATTR",
]
