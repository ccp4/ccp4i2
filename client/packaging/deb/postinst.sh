#!/bin/bash
# Runs as root during `apt install` / `dpkg -i`. Does the two things an AppImage
# cannot, which together let the Chromium sandbox start on locked-down kernels
# (Ubuntu 24.04+ restricts unprivileged user namespaces by default) — so the
# .deb needs NO --no-sandbox and keeps full sandboxing.
#
# NOTE: electron-builder runs its own ${macro} substitution over this file, so it
# must contain NO ${...} shell syntax — use brace-less $VAR only. A bare ${VAR}
# here fails the build with "Macro VAR is not defined".
set -e

APP_DIR=/opt/ccp4i2-django
SANDBOX="$APP_DIR/chrome-sandbox"

# 1. The Chromium SUID sandbox helper must be owned by root with mode 4755.
#    Without this the app aborts with "The SUID sandbox helper binary was found,
#    but is not configured correctly". This alone makes the (privileged) SUID
#    sandbox work regardless of the unprivileged-userns restriction.
if [ -f "$SANDBOX" ]; then
  chown root:root "$SANDBOX" || true
  chmod 4755 "$SANDBOX" || true
fi

# 2. Also install an AppArmor profile granting user namespaces, so the
#    (preferred) namespace sandbox works too on 24.04+. Belt-and-suspenders and
#    a no-op where AppArmor isn't enforcing the restriction.
PROFILE=/etc/apparmor.d/ccp4i2-django
if [ -d /etc/apparmor.d ] && command -v apparmor_parser >/dev/null 2>&1; then
  cat > "$PROFILE" <<'PROF'
abi <abi/4.0>,
include <tunables/global>

profile ccp4i2-django "/opt/ccp4i2-django/ccp4i2-django" flags=(unconfined) {
  userns,
  include if exists <local/ccp4i2-django>
}
PROF
  # Non-fatal: older AppArmor without abi/4.0 will reject this; that's fine
  # because those kernels don't enforce the restriction in the first place.
  apparmor_parser -r -W "$PROFILE" 2>/dev/null || true
fi

exit 0
