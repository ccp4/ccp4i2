"""
Django settings

https://docs.djangoproject.com/en/4.2/topics/settings/
https://docs.djangoproject.com/en/4.2/ref/settings/
"""

import os
from pathlib import Path
from urllib.parse import urlparse, unquote
# DISABLED: Old ccp4i2 import
# from ccp4i2.googlecode import diff_match_patch_py3

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
    print(f"Using ALLOWED_HOSTS from environment: {ALLOWED_HOSTS}")
else:
    # Default hosts for development
    ALLOWED_HOSTS = ["localhost", "127.0.0.1"]
    if DEBUG:
        ALLOWED_HOSTS.append("*")  # Allow all hosts in debug mode
    print(f"Using default ALLOWED_HOSTS: {ALLOWED_HOSTS}")

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
    "corsheaders.middleware.CorsMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

# Azure AD authentication middleware (optional)
# Enable by setting CCP4I2_REQUIRE_AUTH=true along with:
#   AZURE_AD_TENANT_ID - Your Azure AD tenant ID
#   AZURE_AD_CLIENT_ID - Your Azure AD app registration client ID
if os.environ.get("CCP4I2_REQUIRE_AUTH", "").lower() in ("true", "1", "yes"):
    MIDDLEWARE.insert(0, "ccp4i2.middleware.azure_auth.AzureADAuthMiddleware")
    print("Azure AD authentication middleware ENABLED")

REST_FRAMEWORK = {
    "DEFAULT_FILTER_BACKENDS": ["django_filters.rest_framework.DjangoFilterBackend"]
}

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
    print(f"Using CORS_ALLOWED_ORIGINS from environment: {CORS_ALLOWED_ORIGINS}")
else:
    # Default behavior for development
    NEXT_ADDRESS = os.environ.get("NEXT_ADDRESS", "http://localhost:3000")
    if DEBUG:
        CORS_ALLOWED_ORIGINS = [NEXT_ADDRESS]
        CORS_ALLOW_ALL_ORIGINS = True  # Allow all origins in debug mode
        print(
            f"Debug mode: CORS_ALLOWED_ORIGINS={CORS_ALLOWED_ORIGINS}, CORS_ALLOW_ALL_ORIGINS=True"
        )
    else:
        CORS_ALLOWED_ORIGINS = [NEXT_ADDRESS]
        CORS_ALLOW_ALL_ORIGINS = False
        print(f"Production mode: CORS_ALLOWED_ORIGINS={CORS_ALLOWED_ORIGINS}")

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

USER_DIR = Path.home().resolve() / ".ccp4i2"
USER_DIR.mkdir(exist_ok=True)
MEDIA_ROOT = USER_DIR / "files"

# Database configuration with DATABASE_URL support
DATABASE_URL = os.environ.get("DATABASE_URL")

if DATABASE_URL:
    # Parse the DATABASE_URL
    url = urlparse(DATABASE_URL)

    # Check if this is a SQLite URL (from startup script when no DB_HOST is set)
    if url.scheme.startswith("sqlite"):
        # Handle SQLite DATABASE_URL (format: sqlite:///path/to/db.sqlite)
        # The path is everything after sqlite:// (url.path contains the full path)
        db_path = url.path
        DATABASES = {
            "default": {
                "ENGINE": "django.db.backends.sqlite3",
                "NAME": db_path,
            }
        }
        print(f"Using SQLite database from DATABASE_URL: {db_path}")
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
        print(
            f"Using PostgreSQL database: {url.username}:{masked_password}@{url.hostname}:{url.port}/{url.path[1:]} "
            f"with SSL mode: {sslmode}"
        )

        # Log other options if they exist
        if len(db_options) > 2:  # More than just sslmode and connect_timeout
            other_options = {
                k: v
                for k, v in db_options.items()
                if k not in ["sslmode", "connect_timeout"]
            }
            print(f"Additional database options: {other_options}")

else:
    # Default SQLite configuration
    DATABASES = {
        "default": {
            "ENGINE": "django.db.backends.sqlite3",
            "NAME": os.environ.get("CCP4I2_DB_FILE", USER_DIR / "db.sqlite3"),
        }
    }
    print(f"Using SQLite database: {DATABASES['default']['NAME']}")

TIME_ZONE = "UTC"
USE_TZ = True
CCP4I2_PROJECTS_DIR = Path(
    os.environ.get(
        "CCP4I2_PROJECTS_DIR", Path.home().resolve() / ".ccp4i2" / "CCP4X_PROJECTS"
    )
)
CCP4I2_PROJECTS_DIR.mkdir(exist_ok=True)

REST_FRAMEWORK = {
    "DEFAULT_PARSER_CLASSES": (
        "rest_framework.parsers.JSONParser",
        "rest_framework.parsers.FormParser",
        "rest_framework.parsers.MultiPartParser",
    ),
    # Allow unauthenticated API access for development
    # The Next.js frontend proxies requests and doesn't have access to CSRF tokens
    # In production, Azure AD authentication is handled at the container app level
    "DEFAULT_AUTHENTICATION_CLASSES": [],
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

# CCP4I2_ROOT is set by Electron app (packaged: Resources/ccp4i2, dev: project root)
# Fall back to calculating from BASE_DIR for standalone Django usage
CCP4I2_ROOT_ENV = os.environ.get("CCP4I2_ROOT")
if CCP4I2_ROOT_ENV:
    CCP4I2_ROOT = Path(CCP4I2_ROOT_ENV)
else:
    # Standalone Django / Docker: BASE_DIR is ccp4i2/ directory containing qticons/svgicons
    # BASE_DIR = Path(__file__).parent.parent = server/ccp4i2 (where icons are located)
    CCP4I2_ROOT = BASE_DIR

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
        print(f"Warning: Invalid size value '{value_str}', using default {default}")
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

print(
    f"File upload settings: MAX_MEMORY_SIZE={FILE_UPLOAD_MAX_MEMORY_SIZE}, DATA_MAX_MEMORY_SIZE={DATA_UPLOAD_MAX_MEMORY_SIZE}, MAX_FILES={FILE_UPLOAD_MAX_NUMBER_FILES}"
)

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
