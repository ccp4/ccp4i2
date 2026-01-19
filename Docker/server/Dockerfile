# CCP4i2 Django Server Dockerfile (with bundled CCP4)
#
# This Dockerfile builds the Django REST API backend for CCP4i2 using the
# CCP4 base image. CCP4 is bundled in the image rather than mounted at runtime.
#
# Prerequisites:
#   Build the base images first:
#   ./Docker/azure-uksouth/scripts/build-base-image.sh
#   ./Docker/azure-uksouth/scripts/build-arpwarp-image.sh
#
# Build:
#   docker build -t ccp4i2-server -f Docker/server/Dockerfile .
#
# Run:
#   docker run -p 8000:8000 \
#     -v /path/to/projects:/mnt/projects \
#     -e SECRET_KEY=your-secret-key \
#     ccp4i2-server

ARG ACR_LOGIN_SERVER=ccp4acrukbwmx.azurecr.io
ARG CCP4_VERSION=ccp4-20251105
ARG BASE_IMAGE_NAME=ccp4i2/base-arpwarp

FROM ${ACR_LOGIN_SERVER}/${BASE_IMAGE_NAME}:${CCP4_VERSION}

# Set environment variables
# Use azure_extensions.settings to enable Azure-specific features
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/usr/src/app \
    DJANGO_SETTINGS_MODULE=azure_extensions.settings \
    CCP4_PYTHON=/opt/ccp4/bin/ccp4-python

# Create app directory
WORKDIR /usr/src/app

# Copy server code
# Note: Build context should be repository root
COPY server/ .

# Copy Azure extensions Django app (only used in Azure deployment)
# This keeps Azure-specific code separate from the core ccp4i2 distribution
COPY Docker/azure-uksouth/azure_extensions /usr/src/app/azure_extensions

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
COPY apps/compounds/media_views.py /usr/src/app/compounds/media_views.py
COPY apps/compounds/validators.py /usr/src/app/compounds/validators.py

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
    "django-reversion>=5.0" \
    "pytz"

# Install Azure integration dependencies
# Install to main site-packages for container's Python (uvicorn/gunicorn)
RUN pip install --no-cache-dir -r requirements-azure.txt

# Also install packages to a dedicated directory for ccp4-python subprocesses.
# These are packages that ccp4-python needs but aren't in the CCP4 py-packages.
# Use --no-deps for Django apps to avoid duplicating Django itself (which is in py-packages).
# The subprocess PYTHONPATH puts azure_packages AFTER py-packages, so Django comes from py-packages.
RUN pip install --no-cache-dir --target=/usr/src/app/azure_packages \
    --no-deps "django-storages==1.14.4" \
    --no-deps "django-reversion>=5.0" && \
    pip install --no-cache-dir --target=/usr/src/app/azure_packages \
    "azure-storage-blob>=12.26.0"

# Install from pyproject.toml if additional deps needed
RUN pip install --no-cache-dir -e . 2>/dev/null || true

# Create directories for runtime
RUN mkdir -p /mnt/projects

# Copy startup scripts
COPY Docker/scripts/startup-server.sh /usr/src/app/startup.sh
COPY Docker/scripts/startup-worker.sh /usr/src/app/startup-worker.sh
RUN chmod +x /usr/src/app/startup.sh /usr/src/app/startup-worker.sh

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8000/api/ccp4i2/health/ || exit 1

# Default command
ENTRYPOINT ["/usr/src/app/startup.sh"]
