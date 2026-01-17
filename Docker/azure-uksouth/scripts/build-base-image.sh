#!/bin/bash
# Build and push the CCP4 base image
#
# This script builds the base image containing CCP4, which is used as the
# foundation for the server image. Because CCP4 is large (~10GB) and rarely
# changes, separating it into a base image makes subsequent builds much faster.
#
# The CCP4 installation is NOT stored in the git repo. Instead, provide the path
# to a dereferenced tarball which will be extracted to a temporary build context.
#
# Prerequisites:
#   Create a dereferenced CCP4 tarball (symlinks converted to copies):
#   cd /path/to/ccp4data
#   tar --dereference -czf ccp4-20251105-dereferenced.tgz ccp4-20251105
#
# Usage:
#   ./Docker/azure/scripts/build-base-image.sh <tarball-path> [ccp4-version]
#
# Examples:
#   ./Docker/azure/scripts/build-base-image.sh ~/Developer/ccp4data/ccp4-20251105-dereferenced.tgz
#   ./Docker/azure/scripts/build-base-image.sh /path/to/ccp4.tgz ccp4-20251105

set -e

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"

# Load environment
if [ -f "$SCRIPT_DIR/../.env.deployment" ]; then
    source "$SCRIPT_DIR/../.env.deployment"
fi

# Parse arguments
TARBALL_PATH="$1"
CCP4_VERSION="${2:-${CCP4_VERSION:-ccp4-20251105}}"

if [ -z "$TARBALL_PATH" ]; then
    echo "Usage: $0 <tarball-path> [ccp4-version]"
    echo ""
    echo "Arguments:"
    echo "  tarball-path   Path to dereferenced CCP4 tarball"
    echo "  ccp4-version   CCP4 version name (default: ccp4-20251105)"
    echo ""
    echo "Example:"
    echo "  $0 ~/Developer/ccp4data/ccp4-20251105-dereferenced.tgz"
    echo ""
    echo "To create the tarball:"
    echo "  cd ~/Developer/ccp4data"
    echo "  tar --dereference -czf ccp4-20251105-dereferenced.tgz ccp4-20251105"
    exit 1
fi

if [ ! -f "$TARBALL_PATH" ]; then
    echo "ERROR: Tarball not found at $TARBALL_PATH"
    exit 1
fi

# Configuration - ACR_NAME should be set from .env.deployment after infrastructure deployment
ACR_NAME="${ACR_NAME:-}"
ACR_LOGIN_SERVER="${ACR_LOGIN_SERVER:-${ACR_NAME}.azurecr.io}"
IMAGE_NAME="ccp4i2/base"
IMAGE_TAG="${CCP4_VERSION}"

echo "=========================================="
echo "Building CCP4 Base Image"
echo "=========================================="
echo "Tarball: $TARBALL_PATH"
echo "CCP4 Version: $CCP4_VERSION"
echo "Image: $ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG"
echo ""

# Create temporary build context
BUILD_CONTEXT=$(mktemp -d)
echo "Build context: $BUILD_CONTEXT"

# Cleanup on exit
cleanup() {
    echo "Cleaning up build context..."
    rm -rf "$BUILD_CONTEXT"
}
trap cleanup EXIT

# Copy Dockerfile to build context
cp "$REPO_ROOT/Docker/base/Dockerfile.ccp4-base" "$BUILD_CONTEXT/"

# Extract tarball to build context
echo ""
echo "Extracting CCP4 tarball (this may take a few minutes)..."
tar -xzf "$TARBALL_PATH" -C "$BUILD_CONTEXT"

# Verify extraction
if [ ! -d "$BUILD_CONTEXT/$CCP4_VERSION" ]; then
    echo "ERROR: Expected directory $CCP4_VERSION not found after extraction"
    echo "Contents of build context:"
    ls -la "$BUILD_CONTEXT"
    exit 1
fi

# Check CCP4 directory size
CCP4_SIZE=$(du -sh "$BUILD_CONTEXT/$CCP4_VERSION" | cut -f1)
echo "CCP4 directory size: $CCP4_SIZE"
echo ""

# Login to ACR
echo "Logging in to Azure Container Registry..."
az acr login --name "$ACR_NAME"

# Build the image locally
# For large images like this, local build + push is faster than ACR Tasks
echo ""
echo "Building image locally (this will take a while for the first build)..."
echo "Subsequent builds will use layer cache."
echo ""

cd "$BUILD_CONTEXT"

# Build for linux/amd64 since ACR builds run on amd64
# On Apple Silicon Macs, this uses emulation (slower but required)
docker build \
    --platform linux/amd64 \
    --build-arg CCP4_VERSION="$CCP4_VERSION" \
    -t "$ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG" \
    -t "$ACR_LOGIN_SERVER/$IMAGE_NAME:latest" \
    -f Dockerfile.ccp4-base \
    .

echo ""
echo "Build complete. Image size:"
docker images "$ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG" --format "{{.Size}}"

# Push to ACR
echo ""
echo "Pushing to Azure Container Registry..."
echo "This may take a while for the first push (~10GB)..."
echo ""

docker push "$ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG"
docker push "$ACR_LOGIN_SERVER/$IMAGE_NAME:latest"

echo ""
echo "=========================================="
echo "Base image build complete!"
echo "=========================================="
echo ""
echo "Image: $ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG"
echo ""
echo "Next steps:"
echo "  1. Rebuild and deploy the server image:"
echo "     ./Docker/azure-uksouth/scripts/build-and-push.sh server"
echo "     ./Docker/azure-uksouth/scripts/deploy-applications.sh server"
echo ""
