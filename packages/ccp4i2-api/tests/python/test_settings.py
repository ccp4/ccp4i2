"""Minimal Django settings for the ccp4i2-api package's pytest suite.

Used by ``pytest --ds=test_settings`` (configured via the
``DJANGO_SETTINGS_MODULE`` key in ``pyproject.toml``'s
``[tool.pytest.ini_options]``). Includes only what the auth-middleware
tests need: a SQLite in-memory database, the auth + contenttypes apps
for the User model, and a placeholder SECRET_KEY.
"""

DEBUG = True
SECRET_KEY = "test-key-not-used-anywhere-real"

INSTALLED_APPS = [
    "django.contrib.auth",
    "django.contrib.contenttypes",
]

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": ":memory:",
    },
}

USE_TZ = True
