#!/bin/bash
# Build and push the ARP/wARP layer image
#
# This script builds a layer on top of the CCP4 base image, adding ARP/wARP.
# Much faster than rebuilding the entire base image.
#
# Prerequisites:
#   1. CCP4 base image must already exist in ACR
#   2. Create a dereferenced ARP/wARP tarball:
#      tar --dereference -czf arp-warp-8.0-dereferenced.tgz arp_warp_8.0
#
# Usage:
#   ./Docker/azure-uksouth/scripts/build-arpwarp-image.sh <tarball-path> [arpwarp-version]
#
# Examples:
#   ./Docker/azure-uksouth/scripts/build-arpwarp-image.sh ~/Developer/ccp4data/arp-warp-8.0-dereferenced.tgz
#   ./Docker/azure-uksouth/scripts/build-arpwarp-image.sh /path/to/arpwarp.tgz arp_warp_8.0

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
ARPWARP_VERSION="${2:-arp_warp_8.0}"

if [ -z "$TARBALL_PATH" ]; then
    echo "Usage: $0 <tarball-path> [arpwarp-version]"
    echo ""
    echo "Arguments:"
    echo "  tarball-path     Path to dereferenced ARP/wARP tarball"
    echo "  arpwarp-version  ARP/wARP directory name (default: arp_warp_8.0)"
    echo ""
    echo "Example:"
    echo "  $0 ~/Developer/ccp4data/arp-warp-8.0-dereferenced.tgz"
    echo ""
    echo "To create the tarball:"
    echo "  cd /path/to/arpwarp/parent"
    echo "  tar --dereference -czf arp-warp-8.0-dereferenced.tgz arp_warp_8.0"
    exit 1
fi

if [ ! -f "$TARBALL_PATH" ]; then
    echo "ERROR: Tarball not found at $TARBALL_PATH"
    exit 1
fi

# Configuration
ACR_NAME="${ACR_NAME:-}"
ACR_LOGIN_SERVER="${ACR_LOGIN_SERVER:-${ACR_NAME}.azurecr.io}"
CCP4_VERSION="${CCP4_VERSION:-ccp4-20251105}"
IMAGE_NAME="ccp4i2/base-arpwarp"
IMAGE_TAG="${CCP4_VERSION}"

echo "=========================================="
echo "Building ARP/wARP Layer Image"
echo "=========================================="
echo "Tarball: $TARBALL_PATH"
echo "ARP/wARP Version: $ARPWARP_VERSION"
echo "Base Image: $ACR_LOGIN_SERVER/ccp4i2/base:$CCP4_VERSION"
echo "Output Image: $ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG"
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
cp "$REPO_ROOT/Docker/base/Dockerfile.arpwarp" "$BUILD_CONTEXT/"

# Extract tarball to build context
echo ""
echo "Extracting ARP/wARP tarball..."
tar -xzf "$TARBALL_PATH" -C "$BUILD_CONTEXT"

# Verify extraction
if [ ! -d "$BUILD_CONTEXT/$ARPWARP_VERSION" ]; then
    echo "ERROR: Expected directory $ARPWARP_VERSION not found after extraction"
    echo "Contents of build context:"
    ls -la "$BUILD_CONTEXT"
    exit 1
fi

# Check directory size
ARPWARP_SIZE=$(du -sh "$BUILD_CONTEXT/$ARPWARP_VERSION" | cut -f1)
echo "ARP/wARP directory size: $ARPWARP_SIZE"
echo ""

# Login to ACR
echo "Logging in to Azure Container Registry..."
az acr login --name "$ACR_NAME"

# Build the image locally
echo ""
echo "Building image locally..."
echo ""

cd "$BUILD_CONTEXT"

# Build for linux/amd64 since ACR builds run on amd64
docker build \
    --platform linux/amd64 \
    --build-arg ACR_LOGIN_SERVER="$ACR_LOGIN_SERVER" \
    --build-arg CCP4_VERSION="$CCP4_VERSION" \
    --build-arg ARPWARP_VERSION="$ARPWARP_VERSION" \
    -t "$ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG" \
    -t "$ACR_LOGIN_SERVER/$IMAGE_NAME:latest" \
    -f Dockerfile.arpwarp \
    .

echo ""
echo "Build complete. Image size:"
docker images "$ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG" --format "{{.Size}}"

# Push to ACR
echo ""
echo "Pushing to Azure Container Registry..."
echo ""

docker push "$ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG"
docker push "$ACR_LOGIN_SERVER/$IMAGE_NAME:latest"

echo ""
echo "=========================================="
echo "ARP/wARP layer build complete!"
echo "=========================================="
echo ""
echo "Image: $ACR_LOGIN_SERVER/$IMAGE_NAME:$IMAGE_TAG"
echo ""
echo "Next steps:"
echo "  1. Update build-and-push.sh to use base-arpwarp as the base for server"
echo "  2. Build and deploy the server image:"
echo "     ./Docker/azure-uksouth/scripts/build-and-push.sh server"
echo "     ./Docker/azure-uksouth/scripts/deploy-applications.sh"
echo ""
