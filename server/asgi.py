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
