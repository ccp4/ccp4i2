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
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PYTHONPATH=/usr/src/app \
    DJANGO_SETTINGS_MODULE=ccp4i2.config.settings

# Install system dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libpq-dev \
    curl \
    && rm -rf /var/lib/apt/lists/*

# Create app directory
WORKDIR /usr/src/app

# Copy server code
# Note: Build context should be repository root
COPY server/ .

# Install Python dependencies
# Core Django dependencies
RUN pip install --no-cache-dir \
    django>=4.2,<5.0 \
    djangorestframework>=3.14 \
    django-cors-headers>=4.0 \
    django-filter>=23.0 \
    whitenoise>=6.0 \
    uvicorn>=0.20.0 \
    gunicorn>=21.0 \
    psycopg2-binary>=2.9 \
    dj-database-url>=2.0 \
    python-dotenv>=1.0

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
    CMD curl -f http://localhost:8000/api/health/ || exit 1

# Default command
ENTRYPOINT ["/usr/src/app/startup.sh"]
