# CCP4i2 Django Server Dockerfile
#
# This Dockerfile builds the Django REST API backend for CCP4i2.
# It requires CCP4 to be mounted at runtime (not bundled in image).
#
# Build:
#   docker build -t ccp4i2-server -f Docker/server/Dockerfile .
#
# Run (local development):
#   docker run -p 8000:8000 \
#     -v /path/to/ccp4:/mnt/ccp4data/ccp4-9:ro \
#     -v /path/to/projects:/mnt/ccp4data/ccp4i2-projects \
#     -e SECRET_KEY=your-secret-key \
#     ccp4i2-server

FROM python:3.11-slim

# Set environment variables
# Use azure_extensions.settings to enable Azure-specific features
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/usr/src/app \
    DJANGO_SETTINGS_MODULE=azure_extensions.settings

# Install system dependencies
# libgl1, libfreetype6, libxcb-*, etc. are required for ccp4mg/PySide2 imports used by various wrappers
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libpq-dev \
    curl \
    libgl1 \
    libfreetype6 \
    libfontconfig1 \
    libxcb1 \
    libxcb-render0 \
    libxcb-shm0 \
    libxcb-icccm4 \
    libxcb-image0 \
    libxcb-keysyms1 \
    libxcb-randr0 \
    libxcb-render-util0 \
    libxcb-shape0 \
    libxcb-xfixes0 \
    libxcb-xinerama0 \
    libxcb-xkb1 \
    libxkbcommon0 \
    libxkbcommon-x11-0 \
    libx11-6 \
    libxext6 \
    libxrender1 \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /usr/src/app

# Copy server code
# Note: Build context should be repository root
COPY server/ .

# Copy Azure extensions Django app (only used in Azure deployment)
# This keeps Azure-specific code separate from the core ccp4i2 distribution
COPY Docker/azure/azure_extensions /usr/src/app/azure_extensions

# =============================================================================
# Overlay users Django app (platform admin system)
# =============================================================================
COPY apps/users /usr/src/app/users

# =============================================================================
# Overlay compounds Django apps (registry, assays & constructs)
# =============================================================================
COPY apps/compounds/registry /usr/src/app/compounds/registry
COPY apps/compounds/assays /usr/src/app/compounds/assays
COPY apps/compounds/constructs /usr/src/app/compounds/constructs
COPY apps/compounds/__init__.py /usr/src/app/compounds/__init__.py
COPY apps/compounds/urls.py /usr/src/app/compounds/urls.py
COPY apps/compounds/admin_views.py /usr/src/app/compounds/admin_views.py

# Install Python dependencies
# Core Django dependencies (quoted to prevent shell interpretation of < >)
RUN pip install --no-cache-dir \
    "django>=4.2,<5.0" \
    "djangorestframework>=3.14" \
    "django-cors-headers>=4.0" \
    "django-filter>=23.0" \
    "uvicorn>=0.20.0" \
    "gunicorn>=21.0" \
    "psycopg2-binary>=2.9" \
    "dj-database-url>=2.0" \
    "python-dotenv>=1.0" \
    "django-reversion>=5.0"

# Install Azure integration dependencies
RUN pip install --no-cache-dir -r requirements-azure.txt

# Install from pyproject.toml if additional deps needed
RUN pip install --no-cache-dir -e . 2>/dev/null || true

# Create directories for runtime
RUN mkdir -p /mnt/ccp4data/ccp4i2-projects

# Copy startup scripts
COPY Docker/scripts/startup-server.sh /usr/src/app/startup.sh
RUN chmod +x /usr/src/app/startup.sh

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8000/api/ccp4i2/health/ || exit 1

# Default command
ENTRYPOINT ["/usr/src/app/startup.sh"]
