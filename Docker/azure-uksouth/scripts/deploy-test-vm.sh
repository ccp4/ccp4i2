#!/bin/bash
# Deploy or manage the Windows test VM
#
# Usage:
#   ./deploy-test-vm.sh deploy    # Create the VM
#   ./deploy-test-vm.sh start     # Start a deallocated VM
#   ./deploy-test-vm.sh stop      # Deallocate (stop billing)
#   ./deploy-test-vm.sh rdp       # Show RDP connection details
#   ./deploy-test-vm.sh delete    # Delete all VM resources

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/../.env.deployment"

VM_NAME="ccp4i2-bicep-win-test"
BICEP_FILE="$SCRIPT_DIR/../infrastructure/windows-test-vm.bicep"

case "${1:-}" in
  deploy)
    # Get your public IP for RDP access (try IPv4 first, fall back to IPv6)
    echo "Detecting your public IP..."
    MY_IP=$(curl -4 -s --max-time 5 https://ifconfig.me 2>/dev/null || curl -s https://ifconfig.me)
    echo "Your IP: $MY_IP"

    # Determine CIDR suffix: /32 for IPv4, /128 for IPv6
    if [[ "$MY_IP" == *:* ]]; then
      IP_CIDR="${MY_IP}/128"
    else
      IP_CIDR="${MY_IP}/32"
    fi

    echo ""
    echo "Enter a password for the VM admin account (ccp4admin):"
    echo "  (must be 12+ chars with uppercase, lowercase, number, and special char)"
    read -s ADMIN_PASSWORD
    echo ""

    echo "Deploying Windows test VM to $RESOURCE_GROUP..."
    az deployment group create \
      --resource-group "$RESOURCE_GROUP" \
      --template-file "$BICEP_FILE" \
      --parameters \
        adminPassword="$ADMIN_PASSWORD" \
        allowedRdpSourceIp="$IP_CIDR" \
      --query "properties.outputs" \
      --output table

    echo ""
    echo "VM deployed. Connect via Microsoft Remote Desktop:"
    PUBLIC_IP=$(az vm show -d --name "$VM_NAME" -g "$RESOURCE_GROUP" --query publicIps -o tsv)
    echo "  Address: $PUBLIC_IP"
    echo "  Username: ccp4admin"
    echo ""
    echo "Auto-shutdown is set to 7pm GMT. Run './deploy-test-vm.sh stop' to deallocate manually."
    ;;

  start)
    echo "Starting VM $VM_NAME..."
    az vm start --name "$VM_NAME" --resource-group "$RESOURCE_GROUP"
    PUBLIC_IP=$(az vm show -d --name "$VM_NAME" -g "$RESOURCE_GROUP" --query publicIps -o tsv)
    echo "VM started. RDP to: $PUBLIC_IP"
    ;;

  stop)
    echo "Deallocating VM $VM_NAME (this stops billing)..."
    az vm deallocate --name "$VM_NAME" --resource-group "$RESOURCE_GROUP"
    echo "VM deallocated. No compute charges while stopped."
    ;;

  rdp)
    PUBLIC_IP=$(az vm show -d --name "$VM_NAME" -g "$RESOURCE_GROUP" --query publicIps -o tsv)
    echo "RDP connection details:"
    echo "  Address:  $PUBLIC_IP"
    echo "  Username: ccp4admin"
    echo ""
    echo "On macOS, open Microsoft Remote Desktop and add a new PC with the address above."
    ;;

  delete)
    echo "This will delete the VM and all associated resources (disk, NIC, public IP, NSG, subnet)."
    read -p "Are you sure? (y/N) " confirm
    if [[ "$confirm" =~ ^[Yy]$ ]]; then
      echo "Deleting VM..."
      az vm delete --name "$VM_NAME" -g "$RESOURCE_GROUP" --yes
      echo "Deleting NIC..."
      az network nic delete --name "${VM_NAME}-nic" -g "$RESOURCE_GROUP" 2>/dev/null || true
      echo "Deleting public IP..."
      az network public-ip delete --name "${VM_NAME}-pip" -g "$RESOURCE_GROUP" 2>/dev/null || true
      echo "Deleting NSG..."
      az network nsg delete --name "${VM_NAME}-nsg" -g "$RESOURCE_GROUP" 2>/dev/null || true
      echo "Deleting OS disk..."
      az disk delete --name "${VM_NAME}-osdisk" -g "$RESOURCE_GROUP" --yes 2>/dev/null || true
      echo "Done. All VM resources deleted."
    fi
    ;;

  *)
    echo "Usage: $0 {deploy|start|stop|rdp|delete}"
    exit 1
    ;;
esac
