# CCP4i2 Docker Deployment

This directory contains Docker configurations for deploying CCP4i2 in containerized environments.

## Overview

CCP4i2 consists of two main components:
- **Server**: Django REST API backend (`server/`)
- **Client**: Next.js/React web frontend (`client/`)

Both can be deployed as Docker containers for local development, testing, or cloud deployment.

## Directory Structure

```
Docker/
├── server/
│   └── Dockerfile          # Django API server image
├── client/
│   └── Dockerfile          # Next.js web client image
├── azure/
│   ├── container-app.bicep # Azure Container Apps deployment
│   ├── parameters.example.json
│   └── azure_ad_auth.py    # Azure AD authentication middleware
├── scripts/
│   └── startup-server.sh   # Server container entrypoint
├── docker-compose.yml      # Local development compose file
├── .env.example            # Example environment configuration
└── README.md               # This file
```

## Prerequisites

1. **Docker** and **Docker Compose** installed
2. **CCP4 Suite** installed locally (for job execution)
3. For Azure deployment: Azure CLI (`az`) installed and configured

## Quick Start (Local Development)

### 1. Configure Environment

```bash
cd Docker
cp .env.example .env

# Edit .env to set your CCP4 installation path
# CCP4_PATH=/path/to/your/ccp4-9
```

### 2. Start Services

```bash
# From repository root
docker compose -f Docker/docker-compose.yml up

# Or build fresh images
docker compose -f Docker/docker-compose.yml up --build
```

### 3. Access the Application

- **Web Client**: http://localhost:3000
- **API Server**: http://localhost:8000
- **API Docs**: http://localhost:8000/api/

## Shell Container for CCP4 Setup

The shell container provides an interactive environment matching the server,
with **read-write access** to the CCP4 data volume. Use this to:

- Install/build CCP4 for Docker deployment
- Configure CCP4 paths for containerized execution
- Test CCP4 commands in the Docker environment
- Prepare CCP4 on an Azure Files share

### Start Shell Container

```bash
# Interactive session (recommended)
docker compose -f Docker/docker-compose.yml --profile shell run --rm shell

# Or start in background and attach
docker compose -f Docker/docker-compose.yml --profile shell up -d shell
docker compose -f Docker/docker-compose.yml exec shell bash
```

### Build CCP4 for Docker

Once in the shell container:

```bash
# You're now at /usr/src/app with /mnt/ccp4data mounted read-write

# Example: Install CCP4 from a downloaded tarball
cd /mnt/ccp4data
tar xzf /path/to/ccp4-9.0-linux64.tar.gz

# Configure paths if needed
cd ccp4-9
./bin/ccp4.setup-sh

# Test it works
source ./bin/ccp4.setup-sh
ccp4-python -c "import gemmi; print('CCP4 OK')"
```

## Azure Files for CCP4 Data

Two approaches for setting up CCP4 on Azure Files:

1. **AzCopy (Recommended)**: Untar locally, then use azcopy to upload - much faster
2. **SMB Mount**: Mount the share and work directly - slower but interactive

### Option A: Upload with AzCopy (Recommended)

This is much faster than working over an SMB mount.

#### 1. Configure Azure settings in .env

```bash
# In Docker/.env, add your Azure configuration:
AZURE_STORAGE_ACCOUNT=yourstorageaccount
AZURE_RESOURCE_GROUP=yourresourcegroup
AZURE_FILE_SHARE=ccp4data
```

#### 2. Untar CCP4 locally

```bash
# Extract CCP4 to a local directory
cd /tmp
tar xzf ~/Downloads/ccp4-9.0-linux64.tar.gz
# This creates /tmp/ccp4-9/
```

#### 3. Upload to Azure Files

```bash
# From the Docker directory
cd Docker/azure

# Upload the entire CCP4 directory
./azcopy-files.sh /tmp/ccp4-9 ccp4-9

# Or use sync mode (faster for subsequent updates)
./azcopy-files.sh --sync /tmp/ccp4-9 ccp4-9
```

The script will:
- Read Azure config from `.env`
- Log you into Azure if needed
- Generate a time-limited SAS token
- Upload files using azcopy (parallel, resumable)

### Option B: Mount Azure Files Locally

For interactive work or small changes, you can mount the share directly.

#### Mount on macOS

```bash
# Get storage account key
STORAGE_KEY=$(az storage account keys list \
  --resource-group YOUR_RG \
  --account-name YOUR_STORAGE_ACCOUNT \
  --query '[0].value' -o tsv)

# Create mount point
mkdir -p /Volumes/ccp4data

# Mount the share
mount_smbfs "//YOUR_STORAGE_ACCOUNT:${STORAGE_KEY}@YOUR_STORAGE_ACCOUNT.file.core.windows.net/ccp4data" /Volumes/ccp4data
```

#### Mount on Linux

```bash
# Install cifs-utils
sudo apt-get install cifs-utils

# Create mount point
sudo mkdir -p /mnt/azure/ccp4data

# Create credentials file
cat << EOF | sudo tee /etc/smbcredentials/ccp4data.cred
username=YOUR_STORAGE_ACCOUNT
password=YOUR_STORAGE_KEY
EOF
sudo chmod 600 /etc/smbcredentials/ccp4data.cred

# Mount
sudo mount -t cifs //YOUR_STORAGE_ACCOUNT.file.core.windows.net/ccp4data /mnt/azure/ccp4data \
  -o credentials=/etc/smbcredentials/ccp4data.cred,serverino,nosharesock,actimeo=30
```

#### Configure .env for mounted share

```bash
# In Docker/.env
CCP4_DATA_PATH=/Volumes/ccp4data  # or /mnt/azure/ccp4data on Linux
CCP4_PATH=/Volumes/ccp4data/ccp4-9
```

#### Use Shell Container

```bash
# Start shell with Azure Files mounted
docker compose -f Docker/docker-compose.yml --profile shell run --rm shell

# Inside container, CCP4 data is at /mnt/ccp4data
# Any changes persist to Azure Files
```

## Build Images Individually

### Server Image

```bash
# From repository root
docker build -t ccp4i2-server -f Docker/server/Dockerfile .

# Run with CCP4 mounted
docker run -p 8000:8000 \
  -v /path/to/ccp4-9:/mnt/ccp4data/ccp4-9:ro \
  -v $(pwd)/data/projects:/mnt/ccp4data/ccp4i2-projects \
  -e SECRET_KEY=your-secret-key \
  ccp4i2-server
```

### Client Image

```bash
# From repository root
docker build -t ccp4i2-client -f Docker/client/Dockerfile .

# Run
docker run -p 3000:3000 \
  -e NEXT_PUBLIC_API_URL=http://localhost:8000 \
  ccp4i2-client
```

## Environment Variables

### Server

| Variable | Description | Default |
|----------|-------------|---------|
| `SECRET_KEY` | Django secret key (required) | - |
| `DEBUG` | Enable debug mode | `false` |
| `DATABASE_URL` | Database connection URL | SQLite |
| `CCP4_DATA_PATH` | Path to CCP4 data mount | `/mnt/ccp4data` |
| `CCP4I2_PROJECTS_DIR` | Projects storage path | `/mnt/ccp4data/ccp4i2-projects` |
| `ALLOWED_HOSTS` | Comma-separated allowed hosts | `localhost` |
| `CORS_ALLOWED_ORIGINS` | CORS origins | `http://localhost:3000` |
| `USE_GUNICORN` | Use gunicorn (production) | `true` |
| `WORKERS` | Number of gunicorn workers | `2` |

### Client

| Variable | Description | Default |
|----------|-------------|---------|
| `NEXT_PUBLIC_API_URL` | Backend API URL | `http://localhost:8000` |
| `NODE_ENV` | Node environment | `production` |

## Database Options

### SQLite (Default)

For development and single-user deployments:
```bash
DATABASE_URL=sqlite:////mnt/ccp4data/ccp4i2-projects/ccp4i2.sqlite
```

### PostgreSQL

For production deployments:
```bash
# Start with PostgreSQL profile
docker compose -f Docker/docker-compose.yml --profile postgres up

# Or use external PostgreSQL
DATABASE_URL=postgresql://user:password@host:5432/dbname
```

## Azure Deployment

### Prerequisites

1. Azure subscription
2. Azure CLI installed and logged in
3. Resource group created

### Deploy with Bicep

```bash
# Create resource group
az group create --name ccp4i2-rg --location uksouth

# Copy and edit parameters
cp Docker/azure/parameters.example.json Docker/azure/parameters.json
# Edit parameters.json with your values

# Deploy
az deployment group create \
  --resource-group ccp4i2-rg \
  --template-file Docker/azure/container-app.bicep \
  --parameters @Docker/azure/parameters.json
```

### Upload CCP4 to Azure Files

After deployment, upload CCP4 to the created file share:

```bash
# Get storage account key
STORAGE_KEY=$(az storage account keys list \
  --resource-group ccp4i2-rg \
  --account-name ccp4i2storage \
  --query '[0].value' -o tsv)

# Upload CCP4 (use azcopy for large directories)
azcopy copy '/path/to/ccp4-9/*' \
  "https://ccp4i2storage.file.core.windows.net/ccp4data/ccp4-9?$STORAGE_KEY" \
  --recursive
```

### Azure AD Authentication

For enterprise deployments, enable Azure AD authentication:

1. Register an application in Azure AD
2. Configure the middleware (see `azure/azure_ad_auth.py`)
3. Set environment variables:
   - `AZURE_CLIENT_ID`
   - `AZURE_CLIENT_SECRET`
   - `AZURE_TENANT_ID`
   - `AZURE_ADMIN_GROUP_ID`
   - `AZURE_USER_GROUP_ID`

## CCP4 Requirement

CCP4i2 requires the CCP4 crystallographic suite for job execution. The CCP4 installation must be:

1. **Mounted at runtime** - CCP4 is not bundled in Docker images
2. **Read-only access** - Container only reads CCP4 binaries and libraries
3. **Compatible version** - CCP4 8.0+ recommended

### Why CCP4 is Not Bundled

- CCP4 is large (~5GB) and would create unwieldy images
- CCP4 licensing requires separate installation
- Allows using different CCP4 versions without rebuilding images
- File share mounting works well for cloud deployments

## Development Workflow

### Hot Reload (Development)

For development with hot-reload:

```bash
# Server only (for Electron client development)
docker compose -f Docker/docker-compose.yml up server

# Then run Electron client locally
cd client
npm install
npm run start:electron
```

### Rebuilding After Changes

```bash
# Rebuild specific service
docker compose -f Docker/docker-compose.yml build server

# Rebuild and restart
docker compose -f Docker/docker-compose.yml up --build server
```

## Troubleshooting

### CCP4 Not Found

```
WARNING: CCP4 not found at /mnt/ccp4data/ccp4-9/bin/ccp4.setup-sh
```

Ensure `CCP4_PATH` in `.env` points to your CCP4 installation.

### Database Connection Issues

```
ERROR: Failed to construct DATABASE_URL
```

Check that all database environment variables are set correctly.

### Permission Denied on Projects Directory

```bash
# Fix permissions on the mounted projects directory
sudo chown -R 1000:1000 ./data/projects
```

## Security Notes

1. **Never commit `.env` files** with real credentials
2. **Use secrets management** in production (Azure Key Vault, etc.)
3. **Enable HTTPS** for production deployments
4. **Restrict CORS origins** to your actual frontend domain
5. **Use strong passwords** for database and Django secret key
