"""
Django settings

https://docs.djangoproject.com/en/4.2/topics/settings/
https://docs.djangoproject.com/en/4.2/ref/settings/
"""

import os
from pathlib import Path
from urllib.parse import urlparse, unquote

# BASE_DIR is the directory where your Django project is located (containing manage.py)
BASE_DIR = Path(__file__).resolve().parent.parent

# SECURITY WARNING: keep the secret key used in production secret!
# In production, SECRET_KEY must be set via environment variable.
# For local development, a default is provided.
_secret_key_default = "django-insecure-xq@_ci4r3sl+1!3vt5xz5wurncfvfyq^$k5anjsi3+*wb)(5!v"
SECRET_KEY = os.environ.get("SECRET_KEY")
if not SECRET_KEY:
    if os.environ.get("DEBUG", "True").lower() not in ("true", "1", "yes"):
        raise ValueError(
            "SECRET_KEY environment variable is required in production. "
            "Generate one with: python -c \"import secrets; print(secrets.token_urlsafe(50))\""
        )
    SECRET_KEY = _secret_key_default

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = os.environ.get("DEBUG", "True").lower() in ("true", "1", "yes")

# ALLOWED_HOSTS configuration with environment variable support
ALLOWED_HOSTS_ENV = os.environ.get("ALLOWED_HOSTS")
if ALLOWED_HOSTS_ENV:
    # Parse comma-separated list of hosts
    ALLOWED_HOSTS = [
        host.strip() for host in ALLOWED_HOSTS_ENV.split(",") if host.strip()
    ]
else:
    # Default hosts for development
    ALLOWED_HOSTS = ["localhost", "127.0.0.1"]
    if DEBUG:
        ALLOWED_HOSTS.append("*")  # Allow all hosts in debug mode

INSTALLED_APPS = [
    "corsheaders",
    "django_filters",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.staticfiles",
    "ccp4i2.api.config.ApiConfig",
    "ccp4i2.db.config.DbConfig",
    "rest_framework",
]

MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.middleware.gzip.GZipMiddleware",  # Compress large responses
    "corsheaders.middleware.CorsMiddleware",
    "ccp4i2.middleware.corp.CORPMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

# Authentication middleware selection — fail-closed three-way switch.
#
# Exactly one auth middleware is inserted, picked by deployment shape:
#
# 1. Desktop (Electron):  CCP4I2_LOCAL_SESSION_TOKEN is set
#                         → LocalSessionAuthMiddleware (HMAC against
#                           per-launch secret, request.user is the OS user)
# 2. Cloud (production):  CCP4I2_REQUIRE_AUTH=true
#                         → AzureADAuthMiddleware (JWT validation, group
#                           authorization, no dev fallback)
# 3. Local dev:           neither of the above, AND DEBUG=True
#                         → DevAdminMiddleware (auto-creates dev_admin
#                           superuser; DEBUG-gated inside the middleware
#                           too as defence in depth)
# 4. Otherwise:           NO auth middleware installed. Requests fall
#                         through with AnonymousUser; DRF's IsAuthenticated
#                         rejects them with 401. This is the "production
#                         shaped but missing REQUIRE_AUTH" case — strictly
#                         safer than auto-creating a superuser.
#
# The previous "AzureADAuthMiddleware always inserted, branches internally
# on REQUIRE_AUTH" pattern was a misconfiguration backdoor: a production
# deploy with REQUIRE_AUTH accidentally false would auto-assign dev_admin.
# The new layout keeps each middleware single-responsibility and forces
# misconfiguration to fail closed.
if os.environ.get("CCP4I2_LOCAL_SESSION_TOKEN"):
    MIDDLEWARE.insert(
        0,
        "ccp4i2_api.middleware.local_session.LocalSessionAuthMiddleware",
    )
elif os.environ.get("CCP4I2_REQUIRE_AUTH", "").lower() in ("true", "1", "yes"):
    MIDDLEWARE.insert(0, "ccp4i2_api.middleware.azure_ad.AzureADAuthMiddleware")
elif DEBUG:  # noqa: F821 — DEBUG is defined earlier in this settings file.
    MIDDLEWARE.insert(0, "ccp4i2_api.middleware.dev_admin.DevAdminMiddleware")
# else: no auth middleware. AnonymousUser + DRF IsAuthenticated → 401.

# NOTE: REST_FRAMEWORK is defined once, below. There used to be a second
# definition here that set DEFAULT_FILTER_BACKENDS; because the later assignment
# wins, that one was silently dropping the filter backend and rendering every
# viewset's `filterset_fields` inert. The settings are now merged into the single
# block below.

# CORS configuration with environment variable support
CORS_ALLOWED_ORIGINS_ENV = os.environ.get("CORS_ALLOWED_ORIGINS")
if CORS_ALLOWED_ORIGINS_ENV:
    # Parse comma-separated list of origins
    CORS_ALLOWED_ORIGINS = [
        origin.strip()
        for origin in CORS_ALLOWED_ORIGINS_ENV.split(",")
        if origin.strip()
    ]
    CORS_ALLOW_ALL_ORIGINS = False
else:
    # Default behavior for development
    NEXT_ADDRESS = os.environ.get("NEXT_ADDRESS", "http://localhost:3000")
    if DEBUG:
        CORS_ALLOWED_ORIGINS = [NEXT_ADDRESS]
        CORS_ALLOW_ALL_ORIGINS = True  # Allow all origins in debug mode
    else:
        CORS_ALLOWED_ORIGINS = [NEXT_ADDRESS]
        CORS_ALLOW_ALL_ORIGINS = False

ROOT_URLCONF = "ccp4i2.api.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "APP_DIRS": True,
    },
]

# Static and media URLs
# In production (Azure), static files are served by Next.js web container
# In development, Django serves them directly
STATIC_URL = "/djangostatic/"
MEDIA_URL = "/files/"

from ccp4i2.config import preferences as _preferences

# Shared user preferences (~/.ccp4i2/preferences.json). Precedence for every
# setting resolved below is: environment variable > preferences.json > default.
# Cloud (env-driven) is therefore unaffected; the file is the desktop layer that
# the GUI and the i2/i2run CLI both read, so they stay coherent.
_PREFS = _preferences.load_preferences()

USER_DIR = _preferences.ccp4i2_home()
USER_DIR.mkdir(parents=True, exist_ok=True)
MEDIA_ROOT = USER_DIR / "files"

# Database configuration (env DATABASE_URL > preferences "database" > SQLite default)
DATABASE_URL = _preferences.resolve("database", env="DATABASE_URL", prefs=_PREFS)

if DATABASE_URL:
    # Parse the DATABASE_URL
    url = urlparse(DATABASE_URL)

    # Check if this is a SQLite URL (from startup script when no DB_HOST is set)
    if url.scheme.startswith("sqlite"):
        # Handle SQLite DATABASE_URL (format: sqlite:///path/to/db.sqlite)
        # The path is everything after sqlite:// (url.path contains the full path)
        db_path = url.path
        # Windows drive paths arrive as "/C:/..." (urlparse keeps a leading slash
        # before the drive letter); strip it so SQLite gets "C:/...". POSIX paths
        # ("/data/...") have no drive letter and are left unchanged.
        if len(db_path) >= 3 and db_path[0] == "/" and db_path[1].isalpha() and db_path[2] == ":":
            db_path = db_path[1:]
        DATABASES = {
            "default": {
                "ENGINE": "django.db.backends.sqlite3",
                "NAME": db_path,
            }
        }
    else:
        # PostgreSQL configuration
        # Parse query parameters from the URL
        from urllib.parse import parse_qs

        query_params = parse_qs(url.query)

        # Extract SSL mode from query parameters, default to 'require'
        sslmode = query_params.get("sslmode", ["require"])[0]

        # Extract other common PostgreSQL options from query parameters
        connect_timeout = query_params.get("connect_timeout", ["10"])[0]
        application_name = query_params.get("application_name", ["ccp4i2-django"])[0]

        # Build OPTIONS dictionary with extracted parameters
        db_options = {
            "sslmode": sslmode,
            "connect_timeout": int(connect_timeout),
        }

        # Add application_name if specified
        if "application_name" in query_params:
            db_options["application_name"] = application_name

        # Add any other supported PostgreSQL options from query parameters
        # Common options: sslcert, sslkey, sslrootcert, etc.
        for param in ["sslcert", "sslkey", "sslrootcert", "sslcrl"]:
            if param in query_params:
                db_options[param] = query_params[param][0]

        DATABASES = {
            "default": {
                "ENGINE": "django.db.backends.postgresql",
                "NAME": url.path[1:],  # Remove leading slash
                "USER": unquote(url.username),
                "PASSWORD": unquote(url.password),
                "HOST": unquote(url.hostname),
                "PORT": url.port or 5432,
                "OPTIONS": db_options,
            }
        }

        # Log connection details (hide password for security)
        masked_password = "*" * len(url.password) if url.password else "None"

        # Log other options if they exist
        if len(db_options) > 2:  # More than just sslmode and connect_timeout
            other_options = {
                k: v
                for k, v in db_options.items()
                if k not in ["sslmode", "connect_timeout"]
            }

else:
    # Default SQLite configuration
    DATABASES = {
        "default": {
            "ENGINE": "django.db.backends.sqlite3",
            "NAME": os.environ.get("CCP4I2_DB_FILE", USER_DIR / "db.sqlite3"),
        }
    }

TIME_ZONE = "UTC"
USE_TZ = True
CCP4I2_PROJECTS_DIR = Path(
    _preferences.resolve(
        "projectsDir",
        env="CCP4I2_PROJECTS_DIR",
        default=str(USER_DIR / "CCP4X_PROJECTS"),
        prefs=_PREFS,
    )
)
CCP4I2_PROJECTS_DIR.mkdir(parents=True, exist_ok=True)

REST_FRAMEWORK = {
    "DEFAULT_PARSER_CLASSES": (
        "rest_framework.parsers.JSONParser",
        "rest_framework.parsers.FormParser",
        "rest_framework.parsers.MultiPartParser",
    ),
    # Enables every viewset's `filterset_fields` (e.g. JobViewSet ?project=,
    # FileViewSet ?job=). Previously set in a separate REST_FRAMEWORK block that
    # was clobbered by this one, so filterset_fields was inert API-wide.
    "DEFAULT_FILTER_BACKENDS": [
        "django_filters.rest_framework.DjangoFilterBackend",
    ],
    # Authentication: AzureADAuthentication reads from middleware-validated users
    # In development (no auth), this returns None and AllowAny permits access
    # In production, the middleware validates JWT and sets request.user
    "DEFAULT_AUTHENTICATION_CLASSES": [
        "ccp4i2_api.drf.AzureADAuthentication",
    ],
    "DEFAULT_PERMISSION_CLASSES": [
        "rest_framework.permissions.AllowAny",
    ],
}

# Static files settings
STATIC_URL = "/djangostatic/"
STATIC_ROOT = os.path.join(BASE_DIR, "staticfiles")

# Modified settings for Electron
STATICFILES_STORAGE = (
    "django.contrib.staticfiles.storage.StaticFilesStorage"  # Use default storage
)

# CCP4I2_ROOT: location of the ccp4i2 package (for finding wrappers, pipelines, etc.)
# Defaults to BASE_DIR (the ccp4i2/ package directory) which is correct for all
# deployment modes: pip-installed, Electron, Docker, Azure.
# Environment variable override is available for tests or special configurations.
_ccp4i2_root = _preferences.resolve("ccp4i2Root", env="CCP4I2_ROOT", prefs=_PREFS)
CCP4I2_ROOT = Path(_ccp4i2_root) if _ccp4i2_root else BASE_DIR

# Note: Static files (icons, report assets) are now served by Next.js from public/
# Django staticfiles is only used for Electron desktop app where Next.js handles statics


def parse_size_value(value_str, default):
    """Parse size value that might contain expressions like '104857600*1000'"""
    if not value_str:
        return default
    try:
        # Handle simple multiplication expressions safely
        if "*" in value_str and value_str.count("*") == 1:
            parts = value_str.split("*")
            if len(parts) == 2:
                return int(parts[0].strip()) * int(parts[1].strip())
        # Handle simple numeric values
        return int(value_str)
    except (ValueError, IndexError):
        return default


FILE_UPLOAD_MAX_MEMORY_SIZE = parse_size_value(
    os.environ.get("FILE_UPLOAD_MAX_MEMORY_SIZE"), 104857600
)  # Default: 100MB
DATA_UPLOAD_MAX_MEMORY_SIZE = parse_size_value(
    os.environ.get("DATA_UPLOAD_MAX_MEMORY_SIZE"), 104857600
)  # Default: 100MB
FILE_UPLOAD_MAX_NUMBER_FILES = int(
    os.environ.get("FILE_UPLOAD_MAX_NUMBER_FILES", 10)
)  # Default: 10 files


# =============================================================================
# Logging
# =============================================================================
# Suppress noisy django.request 404 warnings for API endpoints where
# "not found" is a normal response (e.g., diagnostic_xml before a job finishes).
# These appear as WARNING:django.request:Not Found: /api/... in the console
# and confuse non-developer users.

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "loggers": {
        "django.request": {
            "level": "ERROR",
        },
    },
}


# =============================================================================
# Security Settings (Production)
# =============================================================================
# These headers protect against common web vulnerabilities.
# Note: HTTPS/SSL settings are NOT enabled here because the Django server
# runs behind the Next.js proxy (internal HTTP). HTTPS termination happens
# at the Azure Container Apps ingress/web container level.

if not DEBUG:
    # Prevent clickjacking attacks
    X_FRAME_OPTIONS = "DENY"

    # Prevent MIME type sniffing
    SECURE_CONTENT_TYPE_NOSNIFF = True

    # Enable XSS filter in browsers (legacy, but still useful)
    SECURE_BROWSER_XSS_FILTER = True
