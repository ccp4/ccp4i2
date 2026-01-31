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
ARG BUILD_TIMESTAMP
ARG GIT_COMMIT

FROM ${ACR_LOGIN_SERVER}/${BASE_IMAGE_NAME}:${CCP4_VERSION}

# Set environment variables
# Use azure_extensions.settings to enable Azure-specific features
# CCP4 env var is inherited from base image (e.g., /mnt/ccp4data/ccp4-20251105)
# CCP4_PYTHON points to the ccp4-python binary in that installation
# BUILD_TIMESTAMP and GIT_COMMIT are for version tracking
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/usr/src/app \
    DJANGO_SETTINGS_MODULE=azure_extensions.settings \
    CCP4_PYTHON=${CCP4}/bin/ccp4-python \
    BUILD_TIMESTAMP=${BUILD_TIMESTAMP} \
    GIT_COMMIT=${GIT_COMMIT}

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
COPY apps/compounds/settings.py /usr/src/app/compounds/settings.py
COPY apps/compounds/formatting.py /usr/src/app/compounds/formatting.py
COPY apps/compounds/admin_views.py /usr/src/app/compounds/admin_views.py
COPY apps/compounds/config_views.py /usr/src/app/compounds/config_views.py
COPY apps/compounds/media_views.py /usr/src/app/compounds/media_views.py
COPY apps/compounds/validators.py /usr/src/app/compounds/validators.py
COPY apps/compounds/utils.py /usr/src/app/compounds/utils.py

# Install ALL Python dependencies to azure_packages for ccp4-python
# This is the ONLY location for Django packages - ccp4-python uses PYTHONPATH to find them
# entrypoint.sh sets PYTHONPATH to include azure_packages before CCP4's paths
RUN pip install --no-cache-dir --target=/usr/src/app/azure_packages \
    "django>=4.2,<5.0" \
    "djangorestframework>=3.14" \
    "django-cors-headers>=4.0" \
    "django-filter>=23.0" \
    "django-storages[azure]>=1.14" \
    "django-reversion>=5.0" \
    "uvicorn>=0.20.0" \
    "gunicorn>=21.0" \
    "psycopg2-binary>=2.9" \
    "dj-database-url>=2.0" \
    "python-dotenv>=1.0" \
    "openpyxl>=3.1.0" \
    "pytz" \
    "azure-storage-blob>=12.26.0" \
    "azure-identity>=1.25.0" \
    "azure-servicebus>=7.14" \
    "azure-keyvault-secrets>=4.10" \
    "azure-monitor-opentelemetry>=1.8" \
    "PyJWT[crypto]>=2.10" \
    "certifi>=2024.0.0" \
    "Pillow>=10.0.0"

# Install from pyproject.toml if additional deps needed
RUN pip install --no-cache-dir -e . 2>/dev/null || true

# Create directories for runtime
RUN mkdir -p /mnt/projects

# Copy startup scripts
# entrypoint.sh: Sets up environment for ALL commands (including interactive shells)
# startup.sh: Runs migrations and starts the Django server
# startup-worker.sh: Starts the background worker
COPY Docker/scripts/entrypoint.sh /usr/src/app/entrypoint.sh
COPY Docker/scripts/startup-server.sh /usr/src/app/startup.sh
COPY Docker/scripts/startup-worker.sh /usr/src/app/startup-worker.sh
RUN chmod +x /usr/src/app/entrypoint.sh /usr/src/app/startup.sh /usr/src/app/startup-worker.sh

# Expose port
EXPOSE 8000

# Health check
HEALTHCHECK --interval=30s --timeout=10s --start-period=60s --retries=3 \
    CMD curl -f http://localhost:8000/api/ccp4i2/health/ || exit 1

# Entrypoint sets up environment for ALL commands (including interactive shells)
# This ensures 'az containerapp exec' gets the same environment as the server
ENTRYPOINT ["/usr/src/app/entrypoint.sh"]
CMD ["/usr/src/app/startup.sh"]
