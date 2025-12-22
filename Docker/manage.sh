#!/bin/bash
# CCP4i2 Docker Management Script
#
# Usage:
#   ./manage.sh <command> [options]
#
# Commands:
#   up        Start all services (server + client)
#   down      Stop all services
#   server    Start only the server (for Electron development)
#   client    Start only the client
#   shell     Start interactive shell container for CCP4 setup
#   logs      View logs from services
#   build     Build/rebuild containers
#   db        Start with PostgreSQL database
#   clean     Remove all containers and volumes
#   mount     Mount Azure Files share (macOS)
#   unmount   Unmount Azure Files share (macOS)

set -e

# Change to Docker directory
cd "$(dirname "$0")"

# Load environment variables
if [ -f .env ]; then
    export $(grep -v '^#' .env | xargs)
fi

COMPOSE_FILE="docker-compose.yml"

# Azure Files configuration (set in .env or override here)
AZURE_STORAGE_ACCOUNT="${AZURE_STORAGE_ACCOUNT:-}"
AZURE_SHARE_NAME="${AZURE_SHARE_NAME:-ccp4data}"
AZURE_MOUNT_POINT="${AZURE_MOUNT_POINT:-/Volumes/ccp4data}"

show_help() {
    echo "CCP4i2 Docker Management"
    echo ""
    echo "Usage: ./manage.sh <command> [options]"
    echo ""
    echo "Commands:"
    echo "  up          Start all services (server + client)"
    echo "  down        Stop all services"
    echo "  server      Start only the Django server"
    echo "  client      Start only the Next.js client"
    echo "  shell       Start interactive shell container"
    echo "  logs [svc]  View logs (optionally for specific service)"
    echo "  build       Build/rebuild containers"
    echo "  db          Start with PostgreSQL database"
    echo "  clean       Remove all containers and volumes"
    echo "  status      Show status of services"
    echo "  mount       Mount Azure Files share (macOS)"
    echo "  unmount     Unmount Azure Files share (macOS)"
    echo ""
    echo "Examples:"
    echo "  ./manage.sh up              # Start full stack"
    echo "  ./manage.sh server          # Start server only"
    echo "  ./manage.sh shell           # Interactive CCP4 shell"
    echo "  ./manage.sh logs server     # View server logs"
    echo "  ./manage.sh mount           # Mount Azure Files"
}

case "$1" in
    up)
        echo "Starting CCP4i2 services..."
        docker compose -f "$COMPOSE_FILE" up "${@:2}"
        ;;
    down)
        echo "Stopping CCP4i2 services..."
        docker compose -f "$COMPOSE_FILE" down "${@:2}"
        ;;
    server)
        echo "Starting Django server..."
        docker compose -f "$COMPOSE_FILE" up server "${@:2}"
        ;;
    client)
        echo "Starting Next.js client..."
        docker compose -f "$COMPOSE_FILE" up client "${@:2}"
        ;;
    shell)
        echo "Starting CCP4i2 shell container..."
        echo "CCP4 data will be available at /mnt/ccp4data"
        docker compose -f "$COMPOSE_FILE" run --rm shell "${@:2}"
        ;;
    logs)
        docker compose -f "$COMPOSE_FILE" logs -f "${@:2}"
        ;;
    build)
        echo "Building containers..."
        docker compose -f "$COMPOSE_FILE" build "${@:2}"
        ;;
    db)
        echo "Starting with PostgreSQL database..."
        docker compose -f "$COMPOSE_FILE" --profile postgres up "${@:2}"
        ;;
    clean)
        echo "Removing all containers and volumes..."
        docker compose -f "$COMPOSE_FILE" down -v --remove-orphans
        ;;
    status)
        docker compose -f "$COMPOSE_FILE" ps
        ;;
    mount)
        # Mount Azure Files share on macOS using Finder
        if [ -z "$AZURE_STORAGE_ACCOUNT" ]; then
            echo "Error: AZURE_STORAGE_ACCOUNT not set in .env"
            echo "Add: AZURE_STORAGE_ACCOUNT=your_storage_account_name"
            exit 1
        fi
        echo "Mounting Azure Files share..."
        echo "Storage Account: $AZURE_STORAGE_ACCOUNT"
        echo "Share: $AZURE_SHARE_NAME"
        echo "Mount Point: $AZURE_MOUNT_POINT"
        echo ""
        echo "Opening Finder to mount SMB share..."
        echo "You will be prompted for credentials:"
        echo "  Username: $AZURE_STORAGE_ACCOUNT"
        echo "  Password: Your storage account key"
        echo ""
        open "smb://${AZURE_STORAGE_ACCOUNT}.file.core.windows.net/${AZURE_SHARE_NAME}"
        echo ""
        echo "After mounting, the share will be available at: $AZURE_MOUNT_POINT"
        echo "Update CCP4_DATA_PATH in .env to: $AZURE_MOUNT_POINT"
        ;;
    unmount)
        if [ -d "$AZURE_MOUNT_POINT" ]; then
            echo "Unmounting $AZURE_MOUNT_POINT..."
            diskutil unmount "$AZURE_MOUNT_POINT" || umount "$AZURE_MOUNT_POINT"
            echo "Unmounted successfully"
        else
            echo "Mount point $AZURE_MOUNT_POINT not found"
        fi
        ;;
    *)
        show_help
        ;;
esac
