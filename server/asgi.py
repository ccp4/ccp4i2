"""
Backward-compatible ASGI entry point.

Docker and legacy scripts reference ``asgi:application`` from the server/ directory.
The canonical location is now ``ccp4i2.config.asgi``.
"""
from ccp4i2.config.asgi import application  # noqa: F401
