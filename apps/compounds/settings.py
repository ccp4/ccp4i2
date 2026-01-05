"""
Compounds app Django settings.

This module extends the core CCP4i2 settings with compound registration
and assay management functionality.

Usage:
    export PYTHONPATH="$PWD:$PWD/apps"
    export DJANGO_SETTINGS_MODULE=compounds.settings

For Azure deployment, this imports azure_extensions.settings which in turn
imports ccp4i2.config.settings. For standalone/desktop deployment, it can
import ccp4i2.config.settings directly.
"""

import os

# Determine which base settings to use
# Azure deployment: use azure_extensions.settings (has staged uploads, etc.)
# Standalone: use ccp4i2.config.settings directly
_use_azure = os.environ.get("COMPOUNDS_USE_AZURE_SETTINGS", "").lower() in ("true", "1", "yes")

if _use_azure:
    try:
        from azure_extensions.settings import *  # noqa: F401, F403
        print("Compounds: Using Azure extensions settings")
    except ImportError:
        print("Warning: COMPOUNDS_USE_AZURE_SETTINGS=true but azure_extensions not found")
        print("Falling back to core CCP4i2 settings")
        from ccp4i2.config.settings import *  # noqa: F401, F403
else:
    from ccp4i2.config.settings import *  # noqa: F401, F403
    print("Compounds: Using core CCP4i2 settings")

# Add compounds apps to INSTALLED_APPS
INSTALLED_APPS = INSTALLED_APPS + [  # noqa: F405
    "compounds.registry",
    "compounds.assays",
    "reversion",
]

# Add compounds URLs to the API
# This is handled via include() in the main urls.py, but we set a flag here
COMPOUNDS_ENABLED = True

# Dev auth configuration
# When DEBUG=True and no other auth is present, auto-login as this user
DEV_USER_EMAIL = os.environ.get("DEV_USER_EMAIL", "dev@localhost")

# Add dev auth middleware for local development
# This must be after AuthenticationMiddleware in the chain
if DEBUG:  # noqa: F405
    # Insert after SessionMiddleware and AuthenticationMiddleware
    _dev_middleware = "ccp4i2.middleware.dev_auth.DevAuthMiddleware"
    if _dev_middleware not in MIDDLEWARE:  # noqa: F405
        # Find AuthenticationMiddleware position and insert after it
        try:
            auth_idx = MIDDLEWARE.index("django.contrib.auth.middleware.AuthenticationMiddleware")  # noqa: F405
            MIDDLEWARE.insert(auth_idx + 1, _dev_middleware)  # noqa: F405
        except ValueError:
            # AuthenticationMiddleware not found, append to end
            MIDDLEWARE.append(_dev_middleware)  # noqa: F405

print(f"Compounds app enabled (registry, assays)")
print(f"  Dev user email: {DEV_USER_EMAIL}")
