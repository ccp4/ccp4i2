# CCP4i2 Linux test VM

A self-contained, throwaway **Ubuntu 24.04** desktop VM in Azure for reproducing
how the CCP4i2 Electron app behaves on a **modern, locked-down kernel** — i.e.
one where `kernel.apparmor_restrict_unprivileged_userns=1` by default. That
restriction is what breaks the AppImage's Chromium sandbox (the SUID-sandbox
crash / blank window seen from a12+), and it does **not** exist on Ubuntu 22.04 —
so a 22.04 rig reproduces nothing. This one does.

It owns its **own resource group and VNet**, with no dependency on any shared or
production network, so it can be created and deleted freely.

## Prerequisites

- `az` logged in to the usual subscription.
- An SSH public key at `~/.ssh/id_ed25519.pub` (or set `SSH_PUBKEY_FILE`).
- The [Microsoft Remote Desktop] client for the GUI (RDP).

## Quick start

```bash
cd deploy/linux-test-vm

# Create the RG + VM (prompts for an admin password; locks SSH/RDP to your IP)
CCP4_VM_PASSWORD='Str0ng-Pass!' ./deploy-linux-test-vm.sh deploy

# Wait for the desktop to finish installing, then install CCP4 + fetch artifacts
./deploy-linux-test-vm.sh wait-ready
./deploy-linux-test-vm.sh setup 'https://.../ccp4-XXXX-linux.tar.gz'

# RDP in to click around; or drive it over SSH
./deploy-linux-test-vm.sh rdp
./deploy-linux-test-vm.sh ssh 'cat /proc/sys/kernel/apparmor_restrict_unprivileged_userns'
```

## Day-to-day

| Command | What it does |
|---|---|
| `deploy` | Create the RG + VM (idempotent-ish; re-running redeploys the template) |
| `start` / `stop` | Boot / deallocate (stop billing) |
| `allow-ip` | Re-point the firewall at your current IP (run after your IP changes) |
| `wait-ready` | Poll until cloud-init has finished provisioning the desktop |
| `ssh [cmd]` | SSH in, or run a one-off command on the VM |
| `rdp` | Print the RDP connection details |
| `setup <CCP4_URL>` | On the VM: install CCP4, test `pip install ccp4i2`, fetch AppImage + `.deb`, install the `.deb` |
| `delete` | Delete the **entire** resource group |

Auto-shutdown is 19:00 GMT to avoid surprise bills; `start` brings it back.

## What to test

On this kernel the `.deb` and the AppImage should diverge:

- **`.deb`** — its `postinst` sets `chrome-sandbox` to `root:root` mode `4755`
  and installs an AppArmor `userns` profile, so the Chromium sandbox starts
  normally. **Expected: launches with a full sandbox, no flags.**
- **AppImage** — can't ship a setuid `chrome-sandbox`, so it falls back to
  `--no-sandbox` (+ `--disable-dev-shm-usage`). **Expected: the portable
  fallback**, with the caveats documented in the app.

## Files

| File | Purpose |
|---|---|
| `linux-test-vm.bicep` | Self-contained infra: RG-local VNet/subnet, NSG, public IP, NIC, VM (Ubuntu 24.04), auto-shutdown |
| `linux-test-vm.cloud-init.yaml` | xfce4 + xrdp desktop, Electron/AppImage runtime libs (handles the 24.04 `t64` renames), software GL |
| `deploy-linux-test-vm.sh` | Lifecycle wrapper around `az` |
| `setup-linux-test-vm.sh` | Runs **on** the VM: CCP4 + ccp4i2 install + artifact fetch/install |

[Microsoft Remote Desktop]: https://apps.microsoft.com/detail/9wzdncrfj3ps
