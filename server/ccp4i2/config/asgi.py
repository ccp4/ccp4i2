"""
ASGI config for ccp4i2 project.

It exposes the ASGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/3.1/howto/deployment/asgi/
"""

import os

# Initialize CCP4 DLL paths for Windows (if needed)
try:
    import ccp4dll
except ImportError:
    pass

from django.core.asgi import get_asgi_application

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "ccp4i2.config.settings")

application = get_asgi_application()

# Under the desktop app, die with it if it dies without a chance to tell us
# (crash, SIGKILL, a debugger stop): an orphaned uvicorn tree bound to a dead
# app's port is the "zombie backends" a Linux tester reported. Started here,
# and only here, so the CLI and job subprocesses are never affected.
from ccp4i2.config.parent_watchdog import start_from_environment as _start_watchdog

_start_watchdog()
