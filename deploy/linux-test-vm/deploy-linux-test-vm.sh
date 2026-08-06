#!/bin/bash
# Deploy or manage the CCP4i2 Ubuntu 24.04 Linux test VM.
#
# A self-contained, throwaway rig for reproducing the Electron app's behaviour on
# a modern locked-down kernel (24.04 restricts unprivileged user namespaces by
# default -- the thing that breaks the AppImage sandbox and that 22.04 can't
# show). It owns its own resource group + VNet, so it has zero dependency on any
# shared/production infrastructure and can be created and deleted freely.
#
# Usage:
#   CCP4_VM_PASSWORD='...' ./deploy-linux-test-vm.sh deploy   # Create RG + VM
#   ./deploy-linux-test-vm.sh start                            # Start (deallocated -> running)
#   ./deploy-linux-test-vm.sh stop                             # Deallocate (stop billing)
#   ./deploy-linux-test-vm.sh ip                               # Show public IP
#   ./deploy-linux-test-vm.sh ssh [cmd...]                     # SSH in (or run a command)
#   ./deploy-linux-test-vm.sh rdp                              # Show RDP details
#   ./deploy-linux-test-vm.sh allow-ip                         # Re-point NSG at your current IP
#   ./deploy-linux-test-vm.sh wait-ready                       # Poll until cloud-init has finished
#   ./deploy-linux-test-vm.sh setup <CCP4_URL>                 # Install CCP4 + test ccp4i2 install
#   ./deploy-linux-test-vm.sh delete                           # Delete the whole resource group
#
# Override any default with an env var (RESOURCE_GROUP, LOCATION, PREFIX,
# ADMIN_USER, SSH_PUBKEY_FILE, VM_SIZE).

set -e -o pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RESOURCE_GROUP="${RESOURCE_GROUP:-ccp4i2-linux-test-rg}"
LOCATION="${LOCATION:-uksouth}"
PREFIX="${PREFIX:-ccp4i2-linux-test}"
ADMIN_USER="${ADMIN_USER:-ccp4admin}"
SSH_PUBKEY_FILE="${SSH_PUBKEY_FILE:-$HOME/.ssh/id_ed25519.pub}"
VM_SIZE="${VM_SIZE:-Standard_D4s_v5}"
VM_NAME="${PREFIX}-vm"
NSG_NAME="${VM_NAME}-nsg"
BICEP_FILE="$SCRIPT_DIR/linux-test-vm.bicep"

detect_ip() {
  local ip
  ip=$(curl -4 -s --max-time 6 https://ifconfig.me 2>/dev/null || curl -s https://ifconfig.me)
  if [[ "$ip" == *:* ]]; then echo "${ip}/128"; else echo "${ip}/32"; fi
}

public_ip() {
  az vm show -d --name "$VM_NAME" -g "$RESOURCE_GROUP" --query publicIps -o tsv
}

case "${1:-}" in
  deploy)
    IP_CIDR="$(detect_ip)"
    echo "Locking access to your current IP: $IP_CIDR"

    if [[ -z "${CCP4_VM_PASSWORD:-}" ]]; then
      echo "Enter a password for the VM admin account ($ADMIN_USER):"
      echo "  (12+ chars: upper, lower, number, special) -- used for RDP + password SSH"
      read -r -s CCP4_VM_PASSWORD
      echo ""
    fi

    if [[ ! -f "$SSH_PUBKEY_FILE" ]]; then
      echo "ERROR: SSH public key not found at $SSH_PUBKEY_FILE" >&2
      echo "       Set SSH_PUBKEY_FILE=/path/to/key.pub or generate one with ssh-keygen." >&2
      exit 1
    fi
    SSH_PUBKEY="$(cat "$SSH_PUBKEY_FILE")"

    echo "Ensuring resource group $RESOURCE_GROUP ($LOCATION) exists..."
    az group create --name "$RESOURCE_GROUP" --location "$LOCATION" -o none

    echo "Deploying Linux test VM to $RESOURCE_GROUP..."
    az deployment group create \
      --resource-group "$RESOURCE_GROUP" \
      --template-file "$BICEP_FILE" \
      --parameters \
        prefix="$PREFIX" \
        vmSize="$VM_SIZE" \
        adminPassword="$CCP4_VM_PASSWORD" \
        adminSshPublicKey="$SSH_PUBKEY" \
        allowedSourceIp="$IP_CIDR" \
      --query "properties.outputs" \
      --output table

    PUBLIC_IP="$(public_ip)"
    echo ""
    echo "VM deployed at $PUBLIC_IP"
    echo "  SSH: ssh $ADMIN_USER@$PUBLIC_IP"
    echo "  RDP: Microsoft Remote Desktop -> $PUBLIC_IP (user $ADMIN_USER)"
    echo ""
    echo "cloud-init installs the desktop in the background (~3-5 min)."
    echo "Poll with:  $0 wait-ready   then:  $0 setup <CCP4_BUILD_URL>"
    echo "Auto-shutdown is 7pm GMT. Run '$0 stop' to deallocate manually."
    ;;

  start)
    az vm start --name "$VM_NAME" --resource-group "$RESOURCE_GROUP"
    echo "Started. Public IP: $(public_ip)"
    echo "Reminder: run '$0 allow-ip' if your IP changed since last time."
    ;;

  stop)
    az vm deallocate --name "$VM_NAME" --resource-group "$RESOURCE_GROUP"
    echo "Deallocated. No compute charges while stopped."
    ;;

  ip)
    public_ip
    ;;

  ssh)
    shift
    exec ssh -o StrictHostKeyChecking=accept-new "$ADMIN_USER@$(public_ip)" "$@"
    ;;

  rdp)
    echo "RDP: Microsoft Remote Desktop -> $(public_ip) (user $ADMIN_USER)"
    ;;

  allow-ip)
    IP_CIDR="$(detect_ip)"
    echo "Re-pointing NSG rules AllowSSH + AllowRDP at $IP_CIDR ..."
    az network nsg rule update -g "$RESOURCE_GROUP" --nsg-name "$NSG_NAME" \
      -n AllowSSH --source-address-prefixes "$IP_CIDR" -o none
    az network nsg rule update -g "$RESOURCE_GROUP" --nsg-name "$NSG_NAME" \
      -n AllowRDP --source-address-prefixes "$IP_CIDR" -o none
    echo "Done."
    ;;

  wait-ready)
    IP="$(public_ip)"
    echo "Polling $VM_NAME for the cloud-init breadcrumb (up to ~10 min)..."
    for _ in $(seq 1 60); do
      if ssh -o StrictHostKeyChecking=accept-new -o ConnectTimeout=8 \
           "$ADMIN_USER@$IP" 'test -f /var/lib/ccp4i2-cloudinit-done' 2>/dev/null; then
        echo "[OK] cloud-init finished; desktop is ready."
        exit 0
      fi
      sleep 10
    done
    echo "[WARN] breadcrumb not seen yet -- cloud-init may still be running." >&2
    exit 1
    ;;

  setup)
    CCP4_URL="${2:-}"
    if [[ -z "$CCP4_URL" ]]; then
      echo "Usage: $0 setup <CCP4_BUILD_URL>" >&2
      exit 2
    fi
    IP="$(public_ip)"
    echo "Copying setup-linux-test-vm.sh to $VM_NAME ..."
    scp -o StrictHostKeyChecking=accept-new \
      "$SCRIPT_DIR/setup-linux-test-vm.sh" "$ADMIN_USER@$IP:/tmp/setup-linux-test-vm.sh"
    echo "Running setup on the VM (CCP4 download + ccp4i2 install)..."
    ssh -o StrictHostKeyChecking=accept-new "$ADMIN_USER@$IP" \
      "bash /tmp/setup-linux-test-vm.sh '$CCP4_URL'"
    ;;

  delete)
    echo "This deletes the ENTIRE resource group '$RESOURCE_GROUP' (VM, disk, VNet,"
    echo "NIC, public IP, NSG -- everything). The rig owns its own RG, so this is clean."
    read -r -p "Are you sure? (y/N) " confirm
    if [[ "$confirm" =~ ^[Yy]$ ]]; then
      az group delete --name "$RESOURCE_GROUP" --yes
      echo "Done."
    fi
    ;;

  *)
    echo "Usage: $0 {deploy|start|stop|ip|ssh|rdp|allow-ip|wait-ready|setup <url>|delete}"
    exit 1
    ;;
esac
